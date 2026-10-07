#!/usr/bin/env python3
import os
import glob
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scvi
import torch

# -------------------------
# Args
# -------------------------
# 1) out_path (base output folder)
# 2) label_col
# 3) unlabeled_category
# 4) batch_col

out_path = sys.argv[1]
label_col = sys.argv[2]
unlabeled_category = sys.argv[3]
batch_col = sys.argv[4]

in_dir = os.path.join(out_path, "hvg_selection")
out_dir = os.path.join(out_path, "trained")
os.makedirs(out_dir, exist_ok=True)

# -------------------------
# Settings
# -------------------------
scvi.settings.seed = 42
torch.set_float32_matmul_precision("medium")
scvi.settings.dl_num_workers = int(os.environ.get("SCVI_DL_NUM_WORKERS", "8"))

if not torch.cuda.is_available():
    raise RuntimeError("No GPU visible (torch.cuda.is_available() is False).")

subset_files = sorted(glob.glob(os.path.join(in_dir, "subset_*.h5ad")))
if not subset_files:
    raise FileNotFoundError(f"No subset_*.h5ad found in {in_dir}")

print(f"Found {len(subset_files)} subset files in {in_dir}")

for subset_path in subset_files:
    base = os.path.basename(subset_path)  # subset_T_NK_ILC.h5ad
    annotation = base.replace("subset_", "").replace(".h5ad", "")
    safe_name = annotation.replace(" ", "_").replace("/", "_")

    print(f"\n=== Training {annotation} ===")
    subset_adata = ad.read_h5ad(subset_path)

    # checks
    if "counts" not in subset_adata.layers:
        raise ValueError(f"{base}: missing layers['counts'] (expected raw counts).")
    if batch_col not in subset_adata.obs.columns:
        raise ValueError(f"{base}: missing batch_col '{batch_col}' in obs.")
    if label_col not in subset_adata.obs.columns:
        raise ValueError(f"{base}: missing label_col '{label_col}' in obs.")
    if "modality" not in subset_adata.obs.columns:
        raise ValueError(f"{base}: missing 'modality' column in obs.")
    
    subset_adata.obs["modality"] = subset_adata.obs["modality"].astype("category")

    # label handling
    labels = subset_adata.obs[label_col].copy().astype("object")
    labels = pd.Series(labels, index=subset_adata.obs_names)

    labels = labels.fillna(unlabeled_category)
    labels = labels.replace({
        "unknown": unlabeled_category,
        "Unknown": unlabeled_category,
    })

    subset_adata.obs[label_col] = pd.Categorical(labels)
    labels = subset_adata.obs[label_col]

    labeled_only = labels[labels != unlabeled_category]
    n_labeled_classes = pd.unique(labeled_only).size

    print(f"Labels: {list(labels.cat.categories)}")
    print(f"n_labeled_classes (excl '{unlabeled_category}'): {n_labeled_classes}")

    if n_labeled_classes < 2:
        raise ValueError(
            f"{base}: SCANVI requires at least 2 labeled classes excluding "
            f"'{unlabeled_category}', but found {n_labeled_classes}."
        )

    # setup + train SCANVI directly
    scvi.model.SCANVI.setup_anndata(
        subset_adata,
        layer="counts",
        batch_key=batch_col,
        labels_key=label_col,
        unlabeled_category=unlabeled_category,
        categorical_covariate_keys=["modality"],
        
    )

    model = scvi.model.SCANVI(
        subset_adata,
        n_latent=30,
        n_hidden=128,
        n_layers=2,
        dispersion="gene",
        gene_likelihood="zinb",
    )

    plan_kwargs = {"weight_decay": 0.0}

    model.train(
        max_epochs=100,
        early_stopping=True,
        early_stopping_patience=15,
        accelerator="gpu",
        devices=1,
        plan_kwargs=plan_kwargs,
    )
    
    # number of epochs actually trained
    if "elbo_train" in model.history:
        n_epochs_trained = len(model.history["elbo_train"])
    else:
        train_keys = [k for k in model.history.keys() if k.endswith("_train")]
        n_epochs_trained = len(model.history[train_keys[0]]) if train_keys else -1

    print(f"✓ {annotation}: trained for {n_epochs_trained} epochs")
    subset_adata.uns["n_epochs_trained"] = int(n_epochs_trained)

    subset_adata.obs["label:scANVI"] = model.predict()
    subset_adata.uns["used_scanvi"] = True

    latent = model.get_latent_representation().astype(np.float32)
    embedding_key = f"X_scanvi-{safe_name}"
    subset_adata.obsm[embedding_key] = latent
    subset_adata.uns["embedding_key"] = embedding_key
    subset_adata.uns["annotation"] = annotation

    out_file = os.path.join(out_dir, f"subset_{safe_name}.h5ad")
    subset_adata.write_h5ad(out_file)
    print(f"✓ wrote {out_file}  (embedding_key={embedding_key}, shape={latent.shape})")

print("\nAll subsets trained.")
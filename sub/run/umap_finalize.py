#!/usr/bin/env python3
import os
import sys
import glob

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csc_matrix

sc.settings.verbosity = 3

# -------------------------
# Args
# -------------------------
# 1) adata_path (global merged_annotated.h5ad)
# 2) out_path   (base output folder containing trained/)
# 3) resolutions (comma separated, e.g. "0.5,1.0,2.0")
# 4) out_path_base (final global h5ad path)
# 5) out_path_cellxgene (final cellxgene h5ad path)
adata_path = sys.argv[1]
out_path = sys.argv[2]
resolutions = [float(x) for x in sys.argv[3].split(",")]
out_path_base = sys.argv[4]
out_path_cellxgene = sys.argv[5]

trained_dir = os.path.join(out_path, "trained")
trained_files = sorted(glob.glob(os.path.join(trained_dir, "subset_*.h5ad")))
if not trained_files:
    raise FileNotFoundError(f"No trained subset files found in: {trained_dir}")

print(f"Reading global adata: {adata_path}")
adata = ad.read_h5ad(adata_path)

n_cells_total = adata.n_obs
all_cells = adata.obs_names
cell_index_map = {cell: idx for idx, cell in enumerate(all_cells)}

print(f"Global cells: {n_cells_total}")
print(f"Found {len(trained_files)} trained subsets")

for subset_path in trained_files:
    base = os.path.basename(subset_path)  # subset_T_NK_ILC.h5ad
    annotation = base.replace("subset_", "").replace(".h5ad", "")
    safe_name = annotation.replace(" ", "_").replace("/", "_")

    print(f"\n=== Postprocess {annotation} ===")
    subset_adata = ad.read_h5ad(subset_path)

    # Figure out which embedding to use
    embedding_key = subset_adata.uns.get("embedding_key", None)
    if embedding_key is None:
        # fallback: try to find an obsm key that looks like X_scvi-... or X_scanvi-...
        candidates = [k for k in subset_adata.obsm_keys() if k.startswith("X_scvi-") or k.startswith("X_scanvi-")]
        if not candidates:
            raise ValueError(f"{base}: Could not find embedding in .obsm and no uns['embedding_key']")
        embedding_key = candidates[0]

    if embedding_key not in subset_adata.obsm:
        raise ValueError(f"{base}: embedding_key '{embedding_key}' not found in subset_adata.obsm")

    # Determine prefix ("scanvi" vs "scvi") from embedding key
    if embedding_key.startswith("X_scanvi-"):
        prefix = "scanvi"
    elif embedding_key.startswith("X_scvi-"):
        prefix = "scvi"
    else:
        prefix = "latent"

    # ---- Merge embedding back into global adata (NaN-filled) ----
    subset_cells = subset_adata.obs_names
    subset_idx = np.array([cell_index_map[c] for c in subset_cells], dtype=int)

    latent_rep = subset_adata.obsm[embedding_key]
    n_dim = latent_rep.shape[1]

    full_embedding = np.full((n_cells_total, n_dim), np.nan, dtype=np.float32)
    full_embedding[subset_idx] = latent_rep.astype(np.float32)
    adata.obsm[embedding_key] = full_embedding

    # ---- Compute UMAP on subset using the stored embedding ----
    sc.pp.neighbors(subset_adata, use_rep=embedding_key)
    sc.tl.umap(subset_adata, random_state=42)

    umap_key = f"{embedding_key}_umap"  # e.g. X_scanvi-T_NK_ILC_umap
    umap_coords = subset_adata.obsm["X_umap"].astype(np.float32)

    full_umap = np.full((n_cells_total, 2), np.nan, dtype=np.float32)
    full_umap[subset_idx] = umap_coords
    adata.obsm[umap_key] = full_umap

    # ---- Leiden clustering on subset ----
    for res in resolutions:
        leiden_key = f"{prefix}-{safe_name}-{res}_leiden"

        # If igraph not installed, Scanpy will fall back.
        sc.tl.leiden(
            subset_adata,
            resolution=res,
            key_added=leiden_key,
            random_state=42,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )

        all_categories = ["Not in subset"] + list(subset_adata.obs[leiden_key].cat.categories)
        full_leiden = pd.Categorical(["Not in subset"] * n_cells_total, categories=all_categories)

        # Fill only subset cells
        for i, cell_name in enumerate(subset_cells):
            full_leiden[cell_index_map[cell_name]] = subset_adata.obs[leiden_key].iloc[i]

        adata.obs[leiden_key] = full_leiden

    print(f"✓ Added embedding: {embedding_key}")
    print(f"✓ Added UMAP: {umap_key}")
    print(f"✓ Added Leiden clusters for resolutions {resolutions}")

# ---- Normalize Unknowns across categorical columns ----
for c in adata.obs.columns:
    if isinstance(adata.obs[c].dtype, pd.CategoricalDtype):
        if "Unknown" not in adata.obs[c].cat.categories:
            adata.obs[c] = adata.obs[c].cat.add_categories(["Unknown"])
    adata.obs[c] = adata.obs[c].fillna("Unknown")
    adata.obs[c] = adata.obs[c].replace("unknown", "Unknown")

# ---- Verification ----
print(f"\n{'='*60}")
print("VERIFICATION")
print(f"{'='*60}")
print(f"Total cells: {adata.n_obs}")

print("\nNew embeddings in obsm:")
for key in adata.obsm.keys():
    if ("scanvi" in key or "scvi" in key):
        shape = adata.obsm[key].shape
        # count non-nan rows if matrix has at least 1 column
        if shape[1] > 0:
            non_nan = (~np.isnan(adata.obsm[key][:, 0])).sum()
        else:
            non_nan = 0
        print(f"  {key}: {shape}, non-NaN cells: {non_nan}")

print("\nNew Leiden clusters in obs:")
for col in adata.obs.columns:
    if (("scanvi" in col or "scvi" in col) and "leiden" in col):
        unique_counts = adata.obs[col].value_counts()
        print(f"\n  {col}:")
        print(f"    Categories: {list(adata.obs[col].cat.categories)}")
        print(f"    'Not in subset' count: {unique_counts.get('Not in subset', 0)}")


# ---- Store raw + normalized ----
adata.layers["raw_counts"] = csc_matrix(adata.X).astype(np.float32)

adata_tmp = adata.copy()
sc.pp.normalize_total(adata_tmp, target_sum=1e4)
sc.pp.log1p(adata_tmp)
adata.layers["normalized"] = csc_matrix(adata_tmp.X).astype(np.float32)
del adata_tmp

# ---- Save global output ----
os.makedirs(os.path.dirname(out_path_base), exist_ok=True)
adata.write(out_path_base)
print(f"\nWrote: {out_path_base}")

# ---- Prepare cellxgene instance ----
adata_cxg = adata.copy()
adata_cxg.X = adata_cxg.layers["normalized"].copy()

integration_methods = ["harmony", "scvi", "scanvi", "scimilarity", "seurat", "bbknn", "combat"]
for integration in integration_methods:
    embedding_key = f"X_{integration}"
    if embedding_key in adata_cxg.obsm:
        adata_cxg.obsm[integration] = adata_cxg.obsm.pop(embedding_key)

for layer in list(adata_cxg.layers.keys()):
    adata_cxg.layers[layer] = csc_matrix(adata_cxg.layers[layer]).astype(np.float32)

adata_cxg.X = csc_matrix(adata_cxg.X).astype(np.float32)

keep = {"scanvi", "scvi"} | {k for k in adata_cxg.obsm_keys() if k.endswith("_umap")}
for k in list(adata_cxg.obsm_keys()):
    if k not in keep:
        del adata_cxg.obsm[k]

os.makedirs(os.path.dirname(out_path_cellxgene), exist_ok=True)
adata_cxg.write(out_path_cellxgene)
print(f"Wrote: {out_path_cellxgene}")
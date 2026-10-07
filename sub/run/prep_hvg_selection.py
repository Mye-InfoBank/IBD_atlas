import sys

from matplotlib import use
import umap
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
#from scdownstream.modules.local.scanpy.leiden.templates import leiden
import scvi
from scvi.model import SCANVI, SCVI
import torch
import platform
from scipy.sparse import csc_matrix
import numpy as np
import scipy as sp
import yaml
from sympy import N
import os

sc.settings.verbosity = 3
scvi.settings.seed = 42

adata_path = sys.argv[1]
out_path = sys.argv[2]
column_split_cell_type = sys.argv[3]
column_label_cell_type = sys.argv[4]
column_label_cell_type_unlabeled = sys.argv[5]
column_batch = sys.argv[6]
resolutions = sys.argv[7].split(",")
resolutions = [float(i) for i in resolutions]
n_hvgs = int(sys.argv[8])

# read in adata
adata = sc.read_h5ad(adata_path)
# out_path_base is the directory out path + 'base_sub_atlas_manual.h5ad'
out_path_base = f"{out_path}/sub_atlas_manual.h5ad"
out_path_cellxgene = f"{out_path}/sub_atlas_manual_cellxgene.h5ad"

# generate adata inner and outer: intersected genes to train on
if "intersection" in adata.var.columns:
    adata_int = adata[:, adata.var["intersection"].values].copy()
    print(f"Using intersection in adata vars intersected genes for scanvi training.")
else:
    adata_int = adata.copy()
    print(f"No intersection in adata vars, using all genes for scanvi training.")

print("n_vars union:", adata.n_vars)
print("n_vars intersection, here the input genes:", adata_int.n_vars)

assert adata_int.n_obs == adata.n_obs
assert (adata_int.obs_names == adata.obs_names).all()

# Get coarse annotations
coarse_annotations = adata_int.obs[column_split_cell_type].unique()
print(f"Coarse annotations: {list(coarse_annotations)}")

# Prepare cell index mapping for NaN filling
all_cells = adata_int.obs_names
cell_index_map = {cell: idx for idx, cell in enumerate(all_cells)}

os.makedirs(f"{out_path}/hvg_selection", exist_ok=True)

for annotation in coarse_annotations:
    print(f"\n=== Processing {annotation} ===")
    
    mask = adata_int.obs[column_split_cell_type] == annotation
    subset_adata = adata_int[mask].copy()
    print(f"Subset shape: {subset_adata.shape}")
    
    subset_adata.layers["counts"] = subset_adata.X.copy()

    if n_hvgs > 0:
        # exact number of HVGs
        sc.pp.highly_variable_genes(
            subset_adata,
            n_top_genes=n_hvgs,
            subset=False,
            layer="counts",
            flavor="seurat_v3",
        )
        subset_adata.X = subset_adata.layers["counts"].copy()
        subset_adata = subset_adata[:, subset_adata.var["highly_variable"]].copy()

    elif n_hvgs == 0:
        # automatic HVG detection by thresholds
        sc.pp.normalize_total(subset_adata, target_sum=1e4)
        sc.pp.log1p(subset_adata)
        sc.pp.highly_variable_genes(
            subset_adata,
            subset=False,
            layer=None,
            flavor="seurat",
        )
        subset_adata.X = subset_adata.layers["counts"].copy()
        subset_adata = subset_adata[:, subset_adata.var["highly_variable"]].copy()

    else:
        # n_hvgs < 0 -> use all genes
        subset_adata.X = subset_adata.layers["counts"].copy()

    print(f"After HVG filtering: {subset_adata.shape}")
    print(f"Number of HVGs: {subset_adata.n_vars}")
    
    if subset_adata.n_vars == 0:
        print(f"WARNING: No HVGs found for {annotation}! Skipping...")
        continue
    
    subset_adata.layers["counts"] = subset_adata.X.copy()
    subset_out_path = f"{out_path}/hvg_selection/subset_{annotation}.h5ad"
    subset_adata.write(subset_out_path)
nextflow.enable.dsl = 2

params.adata_path     = "/nfs/data/COST_IBD/merged_annotated_cell_types.h5ad"
params.out_path       = "/nfs/data/COST_IBD/sub/results"
params.column_split   = "cell_type:coarse"
params.label_col      = "cell_type:tier_2"
params.unlabeled_cat  = "Unknown"
params.batch_col      = "dataset"
params.resolutions    = "0.5,1.0,2.0"
params.n_hvgs         = 0

params.script_prep    = file("${baseDir}/prep_hvg_selection.py")
params.script_train   = file("${baseDir}/scanvi.py")
params.script_post    = file("${baseDir}/umap_finalize.py")

params.scvi_sif       = "/nfs/scratch/nf-core_work/${System.getenv('USER')}/manual_subatlas/apptainer_library/scvi-tools-py3.13-cu12-1.4.1-runtime.sif"

workflow {
  adata_ch     = Channel.value(file(params.adata_path))
  prep_script  = Channel.value(params.script_prep)
  train_script = Channel.value(params.script_train)
  post_script  = Channel.value(params.script_post)
  start_ch     = Channel.value(true)

  prep_done  = PREP_HVG_SELECTION(start_ch, adata_ch, prep_script)
  train_done = SCANVI_GPU(prep_done, train_script)
  UMAP_FINALIZE(train_done, adata_ch, post_script)
}

process PREP_HVG_SELECTION {
  tag "prep_hvg_selection"
  publishDir "${params.out_path}", mode: 'copy', overwrite: true

  input:
    val dummy
    path adata
    path prep_script

  output:
    val true

  script:
  """
  set -euo pipefail
  mkdir -p "${params.out_path}/hvg_selection"

  python3 "${prep_script}" \
    "${adata}" \
    "${params.out_path}" \
    "${params.column_split}" \
    "${params.label_col}" \
    "${params.unlabeled_cat}" \
    "${params.batch_col}" \
    "${params.resolutions}" \
    "${params.n_hvgs}"
  """
}

process SCANVI_GPU {
  tag "scanvi"
  label "process_gpu"
  container "${params.scvi_sif}"
  publishDir "${params.out_path}", mode: 'copy', overwrite: true

  input:
    val dummy
    path train_script

  output:
    val true

  script:
  """
  set -euo pipefail

  mkdir -p "${params.out_path}/trained"

  nvidia-smi || true
  python -c "import torch; print('torch', torch.__version__); print('cuda:', torch.cuda.is_available(), 'torch_cuda:', torch.version.cuda)" || true

  python "${train_script}" \
    "${params.out_path}" \
    "${params.label_col}" \
    "${params.unlabeled_cat}" \
    "${params.batch_col}"
  """
}

process UMAP_FINALIZE {
  tag "umap_finalize"
  publishDir "${params.out_path}", mode: 'copy', overwrite: true

  input:
    val dummy
    path adata
    path post_script

  output:
    path "sub_atlas_manual.h5ad"
    path "sub_atlas_manual_cellxgene.h5ad"

  script:
  """
  set -euo pipefail

  python3 "${post_script}" \
    "${adata}" \
    "${params.out_path}" \
    "${params.resolutions}" \
    "${params.out_path}/sub_atlas_manual.h5ad" \
    "${params.out_path}/sub_atlas_manual_cellxgene.h5ad"

  cp "${params.out_path}/sub_atlas_manual.h5ad" .
  cp "${params.out_path}/sub_atlas_manual_cellxgene.h5ad" .
  """
}
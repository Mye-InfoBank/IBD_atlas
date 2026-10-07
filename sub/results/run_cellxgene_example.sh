#!/usr/bin/env bash
conda activate cellxgene_ibd

cellxgene launch \
/nfs/data/COST_IBD/versions/IBD/12_00_00/build/results/finalized/merged_cellxgene.h5ad \
--host 0.0.0.0 \
--port 5005 \
--backed \
--disable-annotations \
--disable-diffexp \
--title "IBD Atlas"
#!/usr/bin/env bash

#export NXF_SYNTAX_PARSER=v1

# to start well  without cache and newest version (-r dev)
# nextflow run nf-core/scdownstream -profile apptainer,daisybio,gpu -r dev -c nextflow.config

# to start with cache
nextflow run nf-core/scdownstream -profile apptainer,daisybio,gpu -r dev -c nextflow.config -resume 

# to resume specific runs
#nextflow run nf-core/scdownstream -profile apptainer,daisybio,gpu -r dev -c nextflow.config -resume wise_lamarr

# to start from this version on the server !NOT RECOMMENDED!
#nextflow run /nfs/data/COST_IBD/scdownstream -profile apptainer,daisybio,gpu -c nextflow.config -resume
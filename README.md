# IBD Atlas

Scripts and workflows for the construction, integration, and generation of the core IBD atlas and cell-type-specific sub-atlases of the Mye-InfoBank IBD project.

## Core atlas construction

The core atlas is generated using the nf-core/scdownstream workflow. The required input data and metadata are prepared as follows.

### 1. Prepare input datasets

Place all preprocessed `.h5ad` datasets used for atlas construction in the corresponding `data/` directory.

Each independent `.h5ad` file represents one dataset that will subsequently be integrated into the core atlas.

### 2. Prepare sample metadata

Annotate the samples from all input datasets in:

`core/metadata/sample_table.csv`

The table contains the sample-level metadata used for the atlas, including dataset, patient, condition, tissue, disease status, and other available clinical and technical information.

Samples are matched between the individual `.h5ad` datasets and the metadata using `Sample ID`.

### 3. Generate atlas barcodes

For every cell, generate a unique `atlas_barcodes` identifier by combining the dataset identifier with the original cell barcode:

`atlas_barcodes = dataset + "_" + barcode`

This identifier is used to consistently match cells across the different processing and annotation steps.

### 4. Add manual annotations

Cell-type annotations that were manually curated are stored in:

`core/metadata/manual_annotations.csv`

Merge these annotations into the corresponding cells using `atlas_barcodes` and label all cells that do not have a manual annotation as `Unknown`.

### 5. Remove excluded cells

Cells that should not be included in the final atlas are listed in:

`core/metadata/removed_cells_from_atlas.csv`

Remove these cells from the input data by matching their `atlas_barcodes`.

### 6. Configure nf-core/scdownstream

The files required to run the core atlas integration are located in:

`core/run/`

The datasets to be integrated are specified in:

`core/run/samplesheet.csv`

Pipeline parameters and resource settings are specified in:

`core/run/nextflow.config`

An example command for running the pipeline is provided in:

`core/run/run_example.sh`

The prepared datasets are then integrated with nf-core/scdownstream.

## Post-processing of the core atlas

After completion of the nf-core/scdownstream core integration, restore the original names of the columns that were temporarily supplied to the pipeline as `label` and `batch`:

- `label` → rename back to the manually specified annotation column
- `batch` → rename back to the original column used as the batch variable during integration

The resulting AnnData object represents the integrated core IBD atlas and serves as input for the cell-type-specific sub-atlas generation.

## Sub-atlas generation

Cell-type-specific sub-atlases are generated from the integrated core atlas.

The corresponding Nextflow configurations and scripts can be accessed in the `sub/run` folder.
The sub-atlas run is performed using the 

`sub/run/main.nf` Nextflow script, which takes the integrated core atlas as input and generates sub-atlases for each cell type defined in the `cell_type:coarse` column of the core atlas.

The results of the sub-atlas generation are stored in manually specified output folders, which are defined in the Nextflow configuration.

The resulting sub-atlases can be visualized and explored using the provided cellxgene instances, which are generated for each sub-atlas. An example command for running the cellxgene instance is provided in:

 `sub/results/run_cellxgene_example.sh`.
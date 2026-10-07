# versions IBD
## version naming:
Versions are named in 3 steps: xx_yy_zz: 
- x: if the core atlas changed (new dataset)
- y: if the extended atlas changed (not used for IBD)
- z: if the sub atlas changed
if the core changed the numbering of extended and sub will start by 0 again

## versions overview:

### v00–v03: Initial atlas construction and annotation setup
- Initial IBD atlas assembled and manually coarse-annotated.
- Transition from earlier integration approaches to nf-core/scdownstream.
- Inclusion/exclusion of datasets refined; pediatric and datasets without suitable raw counts were removed.
- Metadata harmonization expanded.
- UMAP generation was adjusted to avoid package-related issues.
- Coarse cell-type annotations were manually revised by experts.
- Higher Leiden clustering resolutions were tested for improved subpopulation resolution.
### v04: Integration strategy optimization
- Integration switched to the patient level to better identify problematic patients and reduce inappropriate mixing.
- integration_hvgs = 0 introduced to use automatic HVG selection.
- Three unpublished patients with fewer than 50 cells were removed.
- UMAP calculation was already performed locally within the core-atlas workflow.
### v05: Gene harmonization and technical cleanup
- Gene symbols were harmonized using the HUGO gene symbol unifier.
- Duplicated genes were handled systematically.
- Cell barcodes were made unique before integration.
- Sparse matrices were standardized to memory-efficient numeric types.
- One unusually large/problematic sample was removed.
### v06–v07: Annotation and metadata refinement
- One dataset was removed (Devlin).
- Expert-derived Hamburg annotations were transferred back to the original datasets and then used as input for scANVI integration, replacing heterogeneous publication-level labels.
Extensive sample metadata were added and harmonized.
Missing labels were completed using majority voting within Leiden clusters.
Cell-type annotation columns were reorganized into harmonized annotation levels.
Medication metadata were refined.
Ambient RNA correction with DecontX was re-enabled.
### v08: Standardized scANVI-based atlas and sub-atlas construction
- Latest scdownstream version.
- Batch definition changed from smaller units to dataset-level batches to improve integration and clustering.
Automatic HVG selection was fixed and used correctly.
scANVI became the standard integration/annotation approach for the sub-atlases.
Large-cell-type annotations were manually re-evaluated and updated.
### v09: Core atlas cleanup
- PBMC-derived samples were removed to keep the atlas focused on intestinal tissue.
- Core and sub-atlases were rebuilt.
### v10: Systematic testing of sub-atlas integration strategies
- Tested use of the union/all available genes versus more restricted gene sets.
- Tested scVI pretraining followed by scANVI versus alternative configurations.
- Compared several model/input strategies to improve cell-state separation.
These versions were mainly exploratory and used to optimize the final modeling strategy!
### v11: Final sub-atlas model optimization
- Compared scANVI alone versus scVI + scANVI.
- Compared fixed training of 100 epochs with automatic/unrestricted training.
- Tested inclusion versus exclusion of the Martin dataset.
Best-performing setup:
    - scANVI only
    - 100 epochs
    - Martin dataset excluded
This gave better cell separation and lower annotation entropy.
### v12: Final atlas expansion and biological refinement
- Martin dataset removed because it strongly reduced the intersecting gene set.
- Manually identified artifact clusters were removed.
- One new high-quality 10x scRNA-seq dataset was added:
    - increased total cell number,
    - improved neutrophil representation,
    - improved coverage of epithelial compartments
    - New datasets also introduced fibroblast states not previously represented in the atlas.
- Two single-nucleus RNA-seq datasets were added to broaden the atlas across sequencing modalities and conditions.
- A categorical integration covariate was added to distinguish scRNA-seq versus snRNA-seq.
- For scANVI, the finest available expert-curated annotation, cell_type:tier_2, was used as the label column.


# atlases IBD
## resulting atlases and structure:
The resulting atlases (.h5ad files) are divided in core and sub atlas for each version; these are in seperate folders.

For the core atlas, the 'results/finalized' folder 2 .h5ad files can be accessed:
- 'merged.h5ad': the atlas containing the raw counts
- 'merged_cellxgene.h5ad': the cellxgene instance which can be and is deployed
- (merged.rds: for working with it as R object)
- (merged_metadata.csv: metadata of the atlas as csv)

For the sub atlas, in the specified results folder 2 .h5ad files can be accessed:
- 'sub_atlas_manual.h5ad': the atlas containing the raw counts
- 'sub_atlas_manual_cellxgene.h5ad': the cellxgene instance which can be and is deployed

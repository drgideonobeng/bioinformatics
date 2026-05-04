# Single-Cell Transcriptomic Profiling of the E9.5 Tbx1+ Cardiopharyngeal Microenvironment

## Executive Summary
This report details the computational processing and biological annotation of a single-cell RNA-sequencing (scRNA-seq) dataset derived from an E9.5 Wild-Type mouse embryo (Tbx1-wt-E95).

The experimental design targeted the cardiopharyngeal network. To isolate this specific spatial domain, embryos underwent precise micro-dissection: the upper head (forebrain and midbrain) and the lower trunk were removed prior to FACS isolation of Tbx1+ lineages and their intimately associated microenvironment. Using a custom, parallelized 10-step Seurat pipeline, raw matrices were processed into a fully annotated dataset.

Reflecting the advanced developmental stage at day 9.5, the analysis captured increased biological complexity compared to earlier timepoint(E8.5), revealing 11 distinct cellular populations spanning cardiovascular, ectodermal/neural, and endodermal lineages.

Repository Structure
Plaintext
.
├── config.sh
├── data/
│   └── filtered_feature_bc_matrix/
├── README.md
├── results/
│   ├── objects/
│   ├── plots/
│   └── tables/
└── scripts/
    ├── 01_data_download.R
    ├── 02_create_seurat_obj.R
    ├── 03_qc_visualize.R
    ├── 04_filter_cells.R
    ├── 05_normalize_data.R
    ├── 06_run_dim_reduction.R
    ├── 07_visualize_clusters.R
    ├── 07b_run_clustree.R
    ├── 08_find_markers.R
    ├── 09_plot_marker_genes.R
    └── 10_annotate_clusters.R
    
###  Quality Control and Preprocessing
The initial phase of the pipeline (Scripts 01-04) established robust data integrity for the E9.5 dataset.

Quality Metrics: Violin plots (results/plots/03_qc_violins.pdf) demonstrated high-quality transcript capture. Cells exhibited optimal gene detection rates (nFeature_RNA) and robust transcript counts (nCount_RNA), indicating healthy, viable cells.

Filtering: Cells with elevated mitochondrial read percentages (percent.mt) or outlier feature counts were strictly filtered out to remove dying cells, empty droplets, and multiplet artifacts.

### Dimensionality Reduction and Clustering Optimization
Post-normalization and scaling (Script 05), Principal Component Analysis (PCA) was performed.

Variance Analysis: The Elbow Plot (results/plots/02_elbow_plot.pdf) exhibited a clear inflection point, capturing the high transcriptomic complexity of the maturing E9.5 tissue.

Resolution Validation: To ensure optimal community detection, the clustree algorithm was implemented (Script 07b). The resulting stability network (results/plots/07b_clustree_resolutions.pdf) empirically validated our chosen clustering resolution, identifying the exact "Goldilocks" zone where 11 stable biological populations emerged before the dataset could over-cluster into statistical noise.

Clustering: Graph-based clustering partitioned the cells into 11 distinct transcriptional states, visualized via UMAP (results/plots/07_umap_clusters.pdf).

###  Cell Type Annotation and Lineage Mapping
Differential expression testing (Script 08), accelerated via multicore parallel processing, generated a comprehensive marker gene catalog (results/tables/05_top_cluster_markers_cheatsheet.csv). These markers allowed for the unambiguous mapping of all 11 clusters to their E9.5 anatomical identities (Script 10).

The Mesodermal & Cardiovascular Core
Cluster 0 (Pharyngeal Mesenchyme): Expanding mesenchymal populations in the arches, marked by Ptx3 and Crym.

Cluster 1 (Epicardial Progenitors): Developing epicardium and associated splanchnic mesoderm, indicated by Tbx18 and Lum.

Cluster 3 (Second Heart Field): Multipotent splanchnic mesoderm actively contributing to the heart tube, defined by Hand2 and Osr1.

Cluster 4 (Pharyngeal Mesoderm): Core branchiomeric muscle progenitors within the arches, marked by Tcf21 and Msc.

Cluster 5 (Pharyngeal Arch Endothelium): Maturing vascular networks and arch arteries, defined by definitive endothelial markers Cdh5 (VE-Cadherin) and Cldn5.

Cluster 7 (Cardiomyocytes): Differentiated, contracting heart muscle defined by robust sarcomeric gene expression (Myl7, Myl3, Tnnc1).

The Ectodermal & Neural Lineages
Cluster 8 (Cranial Neural Crest): Migratory neural crest-derived mesenchyme intimately intertwined with the Tbx1 domain, identified by Sox10 and Dlx5.

Cluster 6 (Dorsal Neural Tube): Early neural progenitor populations, defined by Zic1 and Msx1.

Cluster 9 (Hindbrain): Specific retained rhombomere tissue (rhombomeres 3 and 5) from the roof of the pharyngeal arches, precisely marked by Egr2 (Krox20).

Cluster 10 (Cranial Ganglia / Neuroblasts): Actively differentiating neuroblasts migrating into the arches to form cranial nerves, marked by Neurod1.

The Endodermal Hub
Cluster 2 (Pharyngeal Endoderm): The epithelial lining of the pharyngeal pouches, marked by Pax1 and Krt7.

### Key Biological Insights
The annotated UMAP (results/plots/10_umap_annotated.pdf) perfectly captures the spatial and developmental realities of the E9.5 cardiopharyngeal region:

Developmental Progression: Compared to the E8.5 timepoint, the E9.5 dataset shows advanced differentiation. New, mature populations such as Epicardial Progenitors (Cluster 1) and Cranial Ganglia (Cluster 10) have clearly emerged.

Micro-dissection Validation: The presence of specific Hindbrain (Cluster 9) and migrating Cranial Ganglia (Cluster 10) populations, rather than forebrain tissue, perfectly validates the physical micro-dissection. It confirms the successful isolation of the precise junction where the rhombencephalon, cranial nerves, and cardiopharyngeal mesoderm intersect.

Cardiovascular Continuum: The spatial orientation of the UMAP maintains a clear biological trajectory, showing the progression from Second Heart Field and Pharyngeal Mesoderm progenitors toward mature Cardiomyocytes.

### Conclusion and Next Steps
The computational pipeline successfully transformed raw sequencing reads into a highly pure, mathematically validated, and anatomically accurate reference map of the E9.5 Tbx1+ microenvironment.

### Future Directions:
With this highly refined Wild-Type E9.5 object (10_seurat_annotated.rds) finalized, the project is fully prepped for the integration phase. The next analytical step will be to project the Tbx1KO (Knockout) dataset onto this baseline to identify specific lineage disruptions, altered developmental trajectories, and transcriptomic shifts resulting from the loss of Tbx1.

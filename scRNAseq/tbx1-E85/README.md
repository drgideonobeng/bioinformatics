Final Project Report: Single-Cell Transcriptomic Profiling of the E8.5 Tbx1+ Cardiopharyngeal Microenvironment
## Executive Summary
This report details the end-to-end computational processing and biological annotation of a single-cell RNA-sequencing (scRNA-seq) dataset derived from an E8.5 Wild-Type mouse embryo (Bioject : PRJNA704805;Sample GSM5169142). The experimental design specifically targeted the cardiopharyngeal network by micro-dissecting the embryo (removing the head and lower trunk) and performing FACS isolation for Tbx1+ cells and their intimately associated microenvironment.

Using a custom, modular 10-step Seurat pipeline (v4/v5), raw matrices were successfully processed into a fully annotated, biologically interpretable dataset. The analysis revealed 9 distinct cellular populations spanning all three embryonic germ layers, capturing the dynamic transition of multipotent progenitors into definitive cardiac and vascular lineages.

# Quality Control and Preprocessing
The initial phase of the pipeline (Scripts 01-04) focused on data integrity.

Quality Metrics: The violin plots (03_qc_violins.pdf) demonstrated robust transcript capture. The majority of cells exhibited high gene detection rates (nFeature_RNA between 2000 and 6000) and substantial transcript counts (nCount_RNA).

Filtering: Cells with aberrant mitochondrial read percentages (percent.mt) or abnormally low/high feature counts were filtered out to remove dying cells, empty droplets, and potential doublets, establishing a high-confidence dataset for downstream analysis.

# Dimensionality Reduction and Clustering
Post-normalization and scaling (Script 05), Principal Component Analysis (PCA) was performed.

Elbow Plot Analysis: The standard deviation of principal components (02_elbow_plot.pdf) exhibited a clear inflection point around PC 15-20. This indicates a high degree of transcriptomic complexity, appropriate for a developing embryonic region.

Clustering: Using graph-based clustering (Script 07), the cells partitioned into 9 distinct transcriptional states (Clusters 0 through 8), visualized via Uniform Manifold Approximation and Projection (UMAP) (07_umap_clusters.pdf).

# Cell Type Annotation and Lineage Mapping
Differential expression testing (Script 08) generated a comprehensive marker gene catalog (05_top_cluster_markers_cheatsheet.csv). By correlating these statistical markers with known developmental biology, the 9 clusters were unambiguously mapped to their anatomical identities (Script 10).

The Mesodermal & Cardiovascular Core
Cluster 4 (Pharyngeal Mesoderm): The core Tbx1+ multipotent progenitor pool, defined by Tcf21, Nkx2-6, and Edn1.

Cluster 2 (Cardiomyocytes): The maturing Second Heart Field derivatives actively forming the outflow tract and right ventricle, marked by heavy expression of Tnnt2, Actc1, and Hand2.

Cluster 5 (Pharyngeal Arch Endothelium): Early vascular cells forming the arch arteries, identified by Cdh5 (VE-Cadherin), Cldn5, and Sox17.

Cluster 1 (Cranial Paraxial Mesoderm): Somitic/paraxial tissue contributing to head musculature, marked by Meox1 and Tcf15.

Clusters 0 & 3 (Pharyngeal / Transitional Mesenchyme): Fibroblast-like and transitioning support cells expressing matrix components like Col3a1 and Lum.

The Ectodermal Lineages
Cluster 6 (Hindbrain / Neural Tube): Captured due to its spatial proximity dorsal to the arches, expressing classic neural master regulators Pax6, Sox2, and neuronal Tubb3.

Cluster 7 (Cranial Neural Crest): A highly distinct population of migratory cells pouring into the pharyngeal arches to form facial structures and the cardiac septum, beautifully delineated by Sox10, Foxd3, and Tfap2a.

The Endodermal Signaling Hub
Cluster 8 (Pharyngeal Pouch Endoderm): The epithelial lining of the pharyngeal arches, acting as a crucial signaling center. Defined by epithelial markers Epcam, Cdh1 (E-cadherin), and endodermal Foxa1.

# Key Biological Insights
The annotated UMAP (10_umap_annotated.pdf) perfectly mirrors the spatial and developmental realities of the E8.5 cardiopharyngeal region:

Developmental Trajectory: There is a clear transcriptional continuum on the UMAP spanning from the Pharyngeal Mesoderm (progenitors) through the Transitional Mesenchyme, terminating at the Cardiomyocytes. This visually captures the biological differentiation of the Second Heart Field.

Distinct Ectodermal Origins: The Hindbrain and Cranial Neural Crest populations cluster far away from the mesodermal core, reflecting their completely distinct lineage histories prior to migration into the Tbx1 domain.

Epithelial Isolation: The Pharyngeal Pouch Endoderm clusters as a tight, isolated island, consistent with its stable epithelial state compared to the highly dynamic, transitioning mesenchyme around it.

# Conclusion and Future Directions
The pipeline successfully transformed raw sequencing reads into a highly pure, anatomically verified map of the Wild-Type E8.5 Tbx1+ microenvironment. The physical micro-dissection strategy employed by the researchers was highly effective, eliminating extraneous tissue noise and restricting the data to the precise biological niche of interest.

Next Steps:
This object (10_seurat_annotated.rds) serves as an ideal baseline reference. The immediate next phase of the project will involve integrating this Wild-Type reference with the Tbx1-wt-E9.5 dataset. By performing global and local differential expression across these matched clusters, we can identify which genes and pathways are involved in developmental trajectories.

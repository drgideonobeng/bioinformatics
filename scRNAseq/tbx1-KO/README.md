# Single-Cell Transcriptomic Profiling of the E9.5 Tbx1-Knockout Cardiopharyngeal Micrenvironment
## Executive Summary
This report details the end-to-end computational processing and biological annotation of a single-cell RNA-sequencing (scRNA-seq) dataset derived from an E9.5 Tbx1-Knockout mouse embryo (Sample GSM5169145). The experimental design targeted the cardiopharyngeal region following the loss of Tbx1, a master regulatory gene whose deletion results in severe hypoplasia of the pharyngeal arches and cardiac outflow tract anomalies.

Using the custom built Seurat pipeline, raw matrices were successfully processed into a fully annotated dataset. Notably, the analysis revealed 13 distinct cellular populations, a higher degree of cluster fragmentation compared to the 11 clusters observed in the healthy E9.5 Wild-Type baseline. This fractured transcriptomic landscape perfectly captures the perturbed developmental trajectories, stalled progenitor states, and aberrant cellular behavior characteristic of the Tbx1 mutant phenotype.

## Quality Control and Preprocessing
The initial phase of the pipeline rigorously established data integrity for the knockout tissue.

Quality Metrics: The violin plots (03_qc_violins.pdf) demonstrated robust transcript capture across the sample. The majority of cells exhibited optimal gene detection rates (nFeature_RNA largely clustered between 2000 and 5000) and substantial transcript counts (nCount_RNA), indicating a highly viable cellular suspension despite the genetic perturbation.

Filtering: Cells with aberrant mitochondrial read percentages (percent.mt) or abnormal feature counts were strictly filtered out to remove dying cells and empty droplets, establishing a clean, high-confidence dataset for downstream clustering.

## Dimensionality Reduction and Clustering
Post-normalization and scaling, Principal Component Analysis (PCA) and graph-based clustering were performed.

Elbow Plot Analysis: The standard deviation of principal components (02_elbow_plot.pdf) exhibited a clear inflection point around PC 15-20, successfully capturing the high variance of the tissue's underlying biology.

Resolution Validation: The clustree visualization (07b_clustree_resolutions.pdf) tracked the stability of the cellular groupings. The tree validated that the higher degree of fracturing in the KO dataset was mathematically stable, leading to the optimal resolution that yielded 13 clusters.

Clustering: The cells partitioned into 13 distinct transcriptional states (Clusters 0 through 12), visualized via UMAP (07a_umap_clusters.pdf).

## Cell Type Annotation and Lineage Mapping
Differential expression testing yielded a definitive marker gene catalog (05_top_cluster_markers_cheatsheet.csv), allowing the 13 clusters to be mapped to their biological identities (10_umap_annotated.pdf).

## The Mesodermal & Cardiovascular Core

Cluster 12 (Cardiomyocytes): Differentiated, beating heart muscle cells defined by structural sarcomere genes like Myh7 and Mhrt.

Cluster 10 (Endothelium): Definitive vascular endothelial cells lining the perturbed arch arteries and heart tube, marked by Cdh5 and Icam2.

Cluster 7 (Anterior SHF / Arch Mesenchyme): Anterior Second Heart Field progenitors and arch mesenchyme heavily involved in right ventricle and arch artery development, defined by Hand2 and Tnc.

Cluster 9 (Posterior SHF / Epicardial Progenitors): The posterior domain of the SHF and proepicardial organ, heavily marked by Aldh1a2 and Tbx5.

Cluster 6 (Pharyngeal Mesoderm): Core branchiomeric muscle progenitors defined by Msc and Lhx2.

The Ectodermal & Neural Lineages

Cluster 1 (Cranial Neural Crest): Migrating neural crest cells heavily marked by Msx1 and Zic1.

Cluster 3 (Hindbrain): Retained rhombomere tissue from the roof of the pharyngeal arches, precisely patterned by Hoxa2 and Mafb.

Cluster 11 (Cranial Ganglia / Neuroblasts): Actively differentiating neuroblasts forming the cranial nerves, marked by Neurod4 and Nhlh1.

Stroma, Endoderm, & Sub-States

Cluster 0 (Pharyngeal Mesenchyme): A massive pool of undifferentiated mesenchymal fibroblasts expressing structural matrix genes like Lum and Cthrc1.

Cluster 4 (Proliferating / Cycling Cells): A highly distinct group defined entirely by mitotic spindle and cell-cycle genes (Nusap1, Prc1, Aspm).

Cluster 2 (Posterior Endoderm / Mesoderm): Marked by Pax8 and Hoxb1.

Cluster 8 (Pharyngeal Epithelium): Marked by tight junction genes Cldn4 and Crb3.

Cluster 5 (Transitional Mesenchyme): A specialized or stalled intermediate state marked by Crym and Cygb.

## Key Biological Insights
The annotated UMAP (10_umap_annotated.pdf) reveals the severe transcriptomic consequences of the Tbx1 knockout:

Lineage Fragmentation: The emergence of 13 clusters (compared to 11 in the WT) highlights how Tbx1 deletion fractures normal development. Progenitor pools that normally group together smoothly are now splitting into disjointed sub-populations.

Cell Cycle Uncoupling: The emergence of Cluster 4 (Proliferating / Cycling Cells) as an isolated island is a classic hallmark of perturbed datasets. Instead of integrating smoothly into their respective parent lineages, these cells have decoupled, grouping together purely based on their hyper-proliferative stress state.

SHF Bifurcation: The Second Heart Field has clearly split into spatially distinct Anterior (Cluster 7) and Posterior (Cluster 9) domains on the UMAP, suggesting that the normal continuous flow of multi-potent progenitors into the heart tube has been severely disrupted by the loss of Tbx1.

# Conclusion and Future Directions
The computational pipeline successfully mapped the highly perturbed microenvironment of the E9.5 Tbx1 Knockout. By clustering this dataset independently, we have identified the specific progenitor pools that stall, fracture, or hyper-proliferate in the absence of Tbx1.

Next Steps:
While this independent KO map is highly informative, the ultimate goal is to directly quantify the lineage disruptions. The immediate next phase of the project will involve using the RPCI (RISC) algorithm to mathematically project these Knockout cells onto the fully validated Wild-Type E9.5 reference map. This asymmetrical integration will perfectly highlight exactly where the knockout trajectories derail and which definitive cell types (like mature arch endothelium or specific branchiomeric muscles) fail to form.

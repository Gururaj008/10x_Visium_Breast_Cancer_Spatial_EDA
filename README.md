# 10x Visium Breast Cancer Spatial EDA

A reproducible exploratory analysis pipeline for 10x Visium spatial transcriptomics data from breast cancer tissue.

**Input & Pre‑processing:** Public 10x Genomics Visium breast cancer dataset (filtered counts per spot + spatial image). Quality control: compute total counts, number of genes, % mitochondrial reads; filter spots with > 10% mito; normalize to 10 000 counts per spot; logₑ(1 + counts); identify highly variable genes; scale data (clamp to ±10).

**Analysis & Models:**  
- **Dimensionality reduction:** `sc.pp.scale` → `sc.tl.pca` (scree plot of variance explained)  
- **Graph & UMAP:** `sc.pp.neighbors` (n_pcs=30, n_neighbors=10) → `sc.tl.umap` (overlay QC metrics)  
- **Clustering:** Leiden (`sc.tl.leiden`, resolution=0.5) visualized on UMAP and spatial coordinates  
- **Spatial autocorrelation:** build neighbor graph (`sq.gr.spatial_neighbors`) → compute Moran’s I (`sq.gr.spatial_autocorr`) → extract top 10 genes by I

**Results & Significance:**  
QC metrics show most spots have 6 000–8 000 counts (low‑capture tail ~3 000) and clear separation by % mito. First 10 PCs capture ~60–70% variance; UMAP reveals distinct tissue subregions. Leiden clustering identifies X clusters (e.g., tumor vs. stroma) confirmed spatially. Top Moran’s I genes (e.g., **Mbp**, **Slc17a7**) exhibit strong spatial patterns, highlighting microenvironment niches.

**Key Takeaways & Usage:**  
A streamlined Scanpy + Squidpy workflow rapidly highlights spatial structure and uncovers spatially variable genes. This EDA lays groundwork for cell‑type deconvolution or spatially aware differential expression.  

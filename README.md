This project follows a step-by-step reproducible workflow integrating Python, R, and Cytoscape.
## Step 1: Data Preprocessing (Python)</strong></summary>
Prepare the RNA-seq count matrix and sample metadata for downstream DESeq2 analysis.
python scripts/preprocessing_data_colab.py

This script:
Cleans column and index names
Matches metadata and expression samples
Adds Condition labels (Tumor vs. Normal)
Outputs ready-to-use CSVs:
* counts_matrix_deseq2_ready.csv
* sample_metadata_deseq2_ready.csv

## Step 2: Differential Gene Expression Analysis (R)
# Run in RStudio or R terminal
source("scripts/deseq2_analysis.R")

This script:
Loads data and creates DESeq2 object
Filters low-expressed genes
Runs DESeq2 pipeline
Applies FDR < 0.05 and |log2FC| > 1 filter
Exports:
* DEGs_filtered_FDR0.05_Log2FC1.csv
* volcano_plot.png
* heatmap.png (top 50 DEGs)

## Step 3: Functional Enrichment Analysis (Python)</strong></summary>
# Perform GO and KEGG enrichment analysis using g:Profiler and generate a barplot of top enriched terms.
python scripts/gprofiler_enrichment.py

This script:
Uses top DEGs as input
Queries g:Profiler for enriched pathways
Outputs:
* enrichment_results.csv
* top10_enriched_terms.png

## Step 4: Pathway Network Visualization (Cytoscape)</strong></summary>
Visualize functional enrichment results as a network using EnrichmentMap in Cytoscape.

# Manual steps:
Open Cytoscape
Install the EnrichmentMap plugin from App Manager
Go to: Apps → EnrichmentMap → Create Enrichment Map
Upload:
Enrichments File = .gem file from g:Profiler
Gene Sets File = .gmt file from g:Profiler

Set:
* FDR q-value cutoff = 0.05
* Click Build
* Apply a layout (e.g., yFiles Organic)

 Scripts for each phase are inside /scripts/ and ready to run independently.

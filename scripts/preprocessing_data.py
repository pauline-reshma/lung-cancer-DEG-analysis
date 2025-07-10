# Preprocessing RNA-seq data for DESeq2 analysis
import pandas as pd

# Load raw count matrix and metadata
count_data = pd.read_csv("GSE81089_counts.csv", index_col=0)
metadata = pd.read_csv("GSE81089_sample_metadata.csv", index_col=0)

# Fix column naming (if needed)
count_data.columns = count_data.columns.str.replace(".", "-", regex=False)
metadata.index = metadata.index.str.replace(".", "-", regex=False)

# Match and filter samples
common_samples = list(set(count_data.columns).intersection(metadata.index))
count_data = count_data[common_samples]
metadata = metadata.loc[common_samples]

# Add Condition label
metadata["Condition"] = ["Tumor" if "tumor" in s.lower() else "Normal" for s in metadata["Sample_Type"]]

# Save for DESeq2
count_data.to_csv("counts_matrix_deseq2_ready.csv")
metadata.to_csv("sample_metadata_deseq2_ready.csv")

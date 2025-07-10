# Functional Enrichment Analysis using g:Profiler
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from gprofiler import GProfiler

# Load DEG list
deg_list = pd.read_csv("../data/DEG_gene_list.txt", header=None)[0].tolist()

# Run g:Profiler query
gp = GProfiler(return_dataframe=True)
results = gp.profile(organism='hsapiens', query=deg_list)

# Save enrichment results
results.to_csv("../results/enrichment_results.csv", index=False)

# Plot top 10 enriched terms
top_results = results.head(10)
top_results["-log10(p_value)"] = -np.log10(top_results["p_value"])

plt.figure(figsize=(10, 6))
sns.barplot(data=top_results, y='name', x='-log10(p_value)', hue='source')
plt.xlabel("-log10(p-value)")
plt.ylabel("Enriched Terms")
plt.title("Top 10 Enriched GO/KEGG Terms")
plt.legend(title="Source")
plt.tight_layout()
plt.savefig("../results/top10_enriched_terms.png")
plt.show()

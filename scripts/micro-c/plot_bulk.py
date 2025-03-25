#################################
## Identidy marker genes in bulk_atlas ##
##################################

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

 
adata = sc.read_h5ad("/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad")
adata.obs['tissue'] = [name[:-8] for name in adata.obs_names]

sc.tl.rank_genes_groups(adata, 'tissue', method='wilcoxon') #Find marker genes by tissue instead of by leiden clustering (done earlier)
marker_genes = adata.uns['rank_genes_groups']['names']

sc.tl.dendrogram(adata, groupby='tissue')

# Get the top N genes for each cluster
n_top_genes = 5
top_genes = pd.DataFrame(marker_genes).iloc[:n_top_genes]

# Dot plot of top marker genes across clusters

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='tissue',  # Use tissue labels for grouping instead of 'leiden'
    n_genes=4,
    values_to_plot="logfoldchanges", cmap='viridis',
    vmin=-4,
    vmax=4,
    min_logfoldchange=3,
    colorbar_title='log fold change'
)
plt.savefig('/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/figures/17MAR2025_marker_genes_by_tissue.pdf')


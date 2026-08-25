import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.stats import hypergeom
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import multipletests
from gprofiler import GProfiler
from matplotlib.colors import LinearSegmentedColormap
from tqdm import trange
from warnings import filterwarnings


CSV_PATH = r"Data\genes_top5_bottom_custom_annotated_dendro_age_crr.csv" # top, bottom

df = pd.read_csv(CSV_PATH)

og_df = df.copy()
pivot_column = 'p_fdr'
df = df[['source', 'name', 'query', 'p_fdr']]

og_df = og_df.set_index('name')


sub_data = df[df['source'] == "GO:CC"] #"GO:BP"
sub_data.drop(columns=['source'], inplace=True)

row_order = pd.Index(sub_data['name']).drop_duplicates()
col_order = sorted(pd.Index(sub_data['query']).drop_duplicates(),
                   key=lambda x: int(x[6:]))

heatmap = sub_data.fillna(50)
heatmap = heatmap.pivot(index='name', columns='query', values='p_fdr')
heatmap = heatmap.dropna(axis=0, how='any')

heatmap = heatmap.reindex(index=row_order, columns=col_order)

Z = linkage(squareform(pdist(heatmap), force='tomatrix'), 'ward')
dendro = dendrogram(Z, labels=heatmap.index, color_threshold=30, orientation="right")


t = 30
c = 2.0
Zt = Z.copy()
Zt[:,2] = np.arcsinh(Zt[:,2] / c) * c
t_t = np.arcsinh(t / c) * c
fig, ax = plt.subplots(figsize=(6, 20))
dendrogram(Zt, labels=heatmap.index, orientation="right", color_threshold=t_t, ax=ax)
ax.axvline(t_t, ls="--", lw=1)
plt.tight_layout()

fig.savefig(r"F:\gene-MRI_submission\dendrogram_CC_bottom_arcsin_long.pdf", bbox_inches="tight")  # top, bottom

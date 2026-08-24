
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.pyplot as plt
import kaleido

INPUT_XLSX  = 'Data\disgenet_results_celltype.xlsx'   
OUTPUT_PATH = "Data\similarity_celltypes.pdf"

SIG_P_THRESHOLD = 1e-3          
TOP_K = 9                
USE_WEIGHTED_JACCARD = False  

AXIS_RANGE  = [-1.2, 1.2]
BASE_WIDTH  = 10.0               
SIM_THRESHOLD = 0.0            

VALID_DISEASES = {"ASD","BD","SCZ","MDD"}
DISEASE_ORDER  = {"SCZ":0, "BD":1, "MDD":2, "ASD":3}


CUSTOM_COLORS = ["#FFE8AB", "#F8AB58", "#EB5A24","#FF0000"]
cmap = LinearSegmentedColormap.from_list("custom", CUSTOM_COLORS)
norm = Normalize(vmin=0, vmax=1)


df = pd.read_excel(INPUT_XLSX)
expected = {"subtype1", "subtype2", "p_corr"}
if not expected.issubset(df.columns):
    raise ValueError(f"Expected columns {expected}, found {set(df.columns)}")

df["peak"]     = df["subtype1"].astype(str)
df["celltype"] = df["subtype2"].astype(str)
df["disease"]  = df["peak"].str.split().str[0]


df_sig = df[df["p_corr"].astype(float) < SIG_P_THRESHOLD].copy()


if df_sig.empty:
    raise ValueError("No rows pass the significance filter p_corr < 0.001. "
                     "Relax the threshold or check the input.")


df_sig.sort_values(["peak", "p_corr"], ascending=[True, True], inplace=True)
topk_df = df_sig.groupby("peak", as_index=False).head(TOP_K)


peak_to_dis = topk_df.groupby("peak")["disease"].first().to_dict()
peaks_all   = [pk for pk in topk_df["peak"].unique() if peak_to_dis.get(pk) in VALID_DISEASES]
peaks_sorted = sorted(peaks_all, key=lambda pk: (DISEASE_ORDER.get(peak_to_dis.get(pk), 99), pk))

if len(peaks_sorted) < 3:
    raise ValueError(f"Need at least 3 peaks after filtering to draw a chord plot; got {len(peaks_sorted)}.")


peak_to_set = (
    topk_df.groupby("peak")["celltype"]
           .apply(lambda s: set(s.tolist()))
           .to_dict()
)


peak_to_weight = {}
if USE_WEIGHTED_JACCARD:
    for pk, grp in topk_df.groupby("peak"):
        ordered = grp["celltype"].tolist()   
        peak_to_weight[pk] = {ct: (TOP_K - i) for i, ct in enumerate(ordered)}


def jaccard_on_sets(peaks, pk2set):
    n = len(peaks)
    S = np.zeros((n, n), dtype=float)
    for i in range(n):
        Ai = pk2set.get(peaks[i], set())
        for j in range(i, n):
            Aj = pk2set.get(peaks[j], set())
            inter = len(Ai & Aj)
            union = len(Ai | Aj)
            s = (inter / union) if union > 0 else 0.0
            S[i, j] = S[j, i] = s
    return pd.DataFrame(S, index=peaks, columns=peaks)

def rank_weighted_jaccard(peaks, pk2w):
    """
    Weighted Jaccard on rank weights:
      s_ij = sum_c min(w_i(c), w_j(c)) / sum_c max(w_i(c), w_j(c))
    where weights w_i(c) are TOP_K..1 for ranks 1..TOP_K within each peak.
    """
    n = len(peaks)
    S = np.zeros((n, n), dtype=float)
    for i in range(n):
        wi = pk2w.get(peaks[i], {})
        keys_i = set(wi.keys())
        for j in range(i, n):
            wj = pk2w.get(peaks[j], {})
            keys_j = set(wj.keys())
            keys_union = keys_i | keys_j
            if not keys_union:
                s = 0.0
            else:
                num = sum(min(wi.get(c, 0), wj.get(c, 0)) for c in keys_union)
                den = sum(max(wi.get(c, 0), wj.get(c, 0)) for c in keys_union)
                s = 0.0 if den == 0 else (num / den)
            S[i, j] = S[j, i] = s
    return pd.DataFrame(S, index=peaks, columns=peaks)

sim_peak_df = (
    rank_weighted_jaccard(peaks_sorted, peak_to_weight)
    if USE_WEIGHTED_JACCARD else
    jaccard_on_sets(peaks_sorted, peak_to_set)
)


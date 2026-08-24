
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
import matplotlib.pyplot as plt


INPUT_XLSX = r"Data\disgenet_results_bioprocess.xlsx"
SIG_CUTOFF = 3                       
OUTPUT_PATH = "Data\similarity_bioprocess.pdf"
AXIS_RANGE = [-1.2, 1.2]
BASE_WIDTH = 10                    
SIM_THRESHOLD = 0.0                  

DISPLAY_NAMES = {"ASD":"Autism","BD":"Bipolar Disorder","SCZ":"Schizophrenia","MDD":"Major Depression"}

CUSTOM_COLORS =["#DDEEFF", "#99CCFF", "#3388CC", "#003366"]
cmap = LinearSegmentedColormap.from_list("custom", CUSTOM_COLORS)
norm = Normalize(vmin=0, vmax=1)

df = pd.read_excel(INPUT_XLSX)
expected = {"subtype1", "subtype2", "p_corr"}
if not expected.issubset(df.columns):
    raise ValueError(f"Expected columns {expected} in the Excel sheet.")

df["peak"]      = df["subtype1"].astype(str)  
df["disease"]   = df["peak"].str.split().str[0]
df["celltype"]  = df["subtype2"].astype(str)
df["neglog10_p"]= -np.log10(df["p_corr"].astype(float))


def level_from_neglog10(v, sig_cut=SIG_CUTOFF):
    if v < sig_cut: return 0
    if v < 10:      return 1
    if v < 20:      return 2
    if v < 30:      return 3
    return 4


peak_cell = (
    df.groupby(["peak", "celltype"])["neglog10_p"]
      .max()
      .reset_index()
)


neglog_mat_peak = peak_cell.pivot(index="peak", columns="celltype", values="neglog10_p").fillna(0.0)
level_mat_peak  = neglog_mat_peak.applymap(level_from_neglog10).astype(float)


peak_to_disease = df.groupby("peak")["disease"].first().to_dict()
valid_diseases = {"ASD","BD","SCZ","MDD"}
peaks_all = [pk for pk in level_mat_peak.index if peak_to_disease.get(pk) in valid_diseases]


disease_order = {"SCZ":0, "BD":1, "MDD":2, "ASD":3}
peaks_sorted = sorted(peaks_all, key=lambda pk: (disease_order.get(peak_to_disease.get(pk), 99), pk))
level_mat_peak = level_mat_peak.loc[peaks_sorted]

def weighted_jaccard_matrix(X: np.ndarray) -> np.ndarray:
    n, _ = X.shape
    S = np.zeros((n, n), dtype=float)
    for i in range(n):
        xi = X[i]
        for j in range(i, n):
            xj = X[j]
            num = np.minimum(xi, xj).sum()
            den = np.maximum(xi, xj).sum()
            s   = 0.0 if den == 0.0 else (num / den)
            S[i, j] = S[j, i] = s
    return S

sim_peak = weighted_jaccard_matrix(level_mat_peak.values)
sim_peak_df = pd.DataFrame(sim_peak, index=level_mat_peak.index, columns=level_mat_peak.index)



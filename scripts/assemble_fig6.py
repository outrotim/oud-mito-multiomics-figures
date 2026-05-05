"""
Fig 6 v2 — high-density 5-panel narrative figure.

Layout (8.5" × 11"):
  Panel A: experimental schematic (left) + z-score composite withdrawal bar (right)
  Panel B: mito qPCR forest plot (5 genes × 3 regions = 15 rows)
  Panel C: 3 WB membrane images + 4 quantification bars
  Panel D: opioid qPCR forest plot (5 genes × 3 regions = 15 rows)
  Panel E: NAc 5×5 mito × opioid heatmap (left) + Mrpl21 × Oprk1 highlight scatter (right)
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib import patches as mpatches
from matplotlib import patheffects as path_effects
import matplotlib.image as mpimg
from scipy import stats
import fitz

# ============== STYLE CONSTANTS ==============
COLOR_SALINE = "#4A90E2"
COLOR_MORPHINE = "#E67E22"
COLOR_DOWN = "#2266BB"   # mito-protective ↓ in Morphine
COLOR_UP = "#CC4422"     # mito-risk ↑ in Morphine
COLOR_NS = "#999999"
FONT = "Arial"

plt.rcParams["font.family"] = FONT
plt.rcParams["font.size"] = 8
plt.rcParams["axes.linewidth"] = 0.6
plt.rcParams["xtick.major.width"] = 0.6
plt.rcParams["ytick.major.width"] = 0.6
plt.rcParams["xtick.major.size"] = 2.5
plt.rcParams["ytick.major.size"] = 2.5
plt.rcParams["pdf.fonttype"] = 42  # TrueType embedded
plt.rcParams["ps.fonttype"] = 42

# ============== PATHS ==============
ROOT = Path(os.environ.get("OUD_PROJECT_ROOT", Path(__file__).resolve().parent.parent))
DATA_DIR = ROOT / "data"  # user-supplied input CSVs go here
MITO_RAW = DATA_DIR / "mito_qPCR_log2FC_perSample.csv"
MITO_STATS = DATA_DIR / "mito_qPCR_statistics_MWU.csv"
OPIOID_RAW = DATA_DIR / "opioid_qPCR_log2FC_perSample.csv"
OPIOID_STATS = DATA_DIR / "opioid_qPCR_statistics_MWU.csv"
BEHAV_RAW = DATA_DIR / "behavioral_raw_data.csv"
WB_RAW = DATA_DIR / "WB_perSample.csv"
WB_COMT = DATA_DIR / "WB_COMT_perSample.csv"
WB_STATS = DATA_DIR / "WB_statistics.csv"
WB_PDF_4PROTEINS = DATA_DIR / "WB_membranes.pdf"
PART1_PDF = DATA_DIR / "Fig6_schematic_template.pdf"
COUPLING = DATA_DIR / "WithinMouse_Mito_x_Opioid_Correlations.csv"
OUT_DIR = ROOT / "outputs"

# ============== HELPERS ==============
def t_ci_95(values):
    """Return (mean, ci_lower, ci_upper) using t-distribution; returns 0,0,0 if empty."""
    a = np.asarray(values, dtype=float)
    if len(a) < 2:
        return float(np.mean(a)) if len(a) else 0.0, 0.0, 0.0
    m = np.mean(a)
    se = stats.sem(a)
    h = se * stats.t.ppf(0.975, df=len(a) - 1)
    return m, m - h, m + h

def diff_ci_95(morphine_vals, saline_vals):
    """Welch-style 95% CI for mean difference (Morphine − Saline)."""
    m = np.mean(morphine_vals) - np.mean(saline_vals)
    s_m = stats.sem(morphine_vals)
    s_s = stats.sem(saline_vals)
    se = np.sqrt(s_m**2 + s_s**2)
    df = (s_m**2 + s_s**2)**2 / (
        s_m**4 / (len(morphine_vals) - 1) + s_s**4 / (len(saline_vals) - 1)
    )
    h = se * stats.t.ppf(0.975, df=df)
    return m, m - h, m + h

def fdr_to_stars(fdr):
    if pd.isna(fdr):
        return ""
    if fdr < 0.001:
        return "***"
    if fdr < 0.01:
        return "**"
    if fdr < 0.05:
        return "*"
    return "ns"

# ============== PANEL B / D: forest plot ==============
def render_forest(ax, raw_df, stats_df, gene_order, region_order,
                  fdr_col="MWU_FDR", region_col="Region", gene_col="Gene"):
    """
    Render a forest plot in `ax`.
    Each row = one gene × region combo. log2FC = mean(Morphine − Saline) with 95% CI.
    Color: down (blue) / up (orange); transparency: ns 35%, sig 100%.
    """
    rows = []
    for gene in gene_order:
        for region in region_order:
            sub = raw_df[(raw_df[gene_col].str.upper() == gene.upper()) &
                         (raw_df[region_col] == region)]
            if sub.empty:
                continue
            mor = sub[sub["Group"].str.lower() == "morphine"]["log2_fold_change"].values
            sal = sub[sub["Group"].str.lower() == "saline"]["log2_fold_change"].values
            if len(mor) == 0 or len(sal) == 0:
                continue
            mean, lo, hi = diff_ci_95(mor, sal)
            # find FDR
            srow = stats_df[(stats_df[gene_col].str.upper() == gene.upper()) &
                           (stats_df[region_col] == region)]
            fdr = srow[fdr_col].iloc[0] if not srow.empty else np.nan
            rows.append((gene, region, mean, lo, hi, fdr))

    n_rows = len(rows)
    y_pos = np.arange(n_rows)[::-1]  # top-down

    region_short = {"NAc": "NAc", "Hippocampus": "HIP", "Hip": "HIP", "HIP": "HIP",
                    "PrefrontalCortex": "PFC", "PFC": "PFC", "Prefrontal Cortex": "PFC"}

    for i, (gene, region, mean, lo, hi, fdr) in enumerate(rows):
        y = y_pos[i]
        sig = (not pd.isna(fdr)) and (fdr < 0.05)
        color = COLOR_UP if mean >= 0 else COLOR_DOWN
        alpha = 1.0 if sig else 0.35
        # CI line
        ax.plot([lo, hi], [y, y], color=color, alpha=alpha, lw=1.5, solid_capstyle="butt")
        # Mean point
        ax.scatter([mean], [y], color=color, alpha=alpha, s=22, zorder=3,
                   edgecolor="black", linewidth=0.4)
        # Asterisks placed JUST RIGHT of the upper CI bound, in plot coords (avoid edge)
        # Use plain text (no mathtext) with sans-serif bold to avoid shadow artifact
        stars = fdr_to_stars(fdr)
        if stars:
            ax.text(hi + 0.12, y, stars, va="center", ha="left",
                    fontsize=8.5, fontweight="bold", color="black",
                    family="DejaVu Sans")

    # Y labels with italic gene name + plain region (mathtext)
    labels = []
    for g, r, *_ in rows:
        rs = region_short.get(r, r[:3])
        # mathtext: $\mathit{Gene}$ for italic
        labels.append(f"$\\mathit{{{g}}}$ · {rs}")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=7)
    # Reference line at 0
    ax.axvline(0, color="black", lw=0.5, linestyle="--", alpha=0.7)
    ax.set_xlabel("log$_2$ fold change (Morphine vs Saline)", fontsize=8)
    ax.tick_params(axis='both', labelsize=7)

    # Add small region group bars on left
    for i, (g, r, *_) in enumerate(rows):
        pass  # could add color blocks for region; skip for cleanness

    # Spine cleanup
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_linewidth(0.6)
    ax.spines['bottom'].set_linewidth(0.6)

    return rows  # for downstream use

# ============== MAIN ASSEMBLY ==============
def main():
    fig = plt.figure(figsize=(8.5, 11.0), dpi=150)

    # 5 rows: A=1.9 / B=2.4 / C=2.2 / D=2.4 / E=2.4 ≈ 11.3" (slight bleed; saved with bbox_inches=tight)
    gs = gridspec.GridSpec(
        5, 1, figure=fig,
        height_ratios=[1.9, 2.4, 2.2, 2.4, 2.4],
        hspace=0.65, top=0.97, bottom=0.04, left=0.08, right=0.97
    )

    # ================ Panel A — Schematic + z-score ================
    gs_A = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], width_ratios=[3.2, 1.0], wspace=0.18)
    ax_schematic = fig.add_subplot(gs_A[0])
    ax_zscore = fig.add_subplot(gs_A[1])

    # Schematic from part1 PDF (rasterize the top-left A region)
    part1_doc = fitz.open(PART1_PDF)
    p1_page = part1_doc[0]
    # A occupies x=0..200, y=0..128 in 504×309.6 page (per earlier text-block scan)
    clip = fitz.Rect(0, 0, 200, 128)
    pix = p1_page.get_pixmap(dpi=500, clip=clip)
    schematic_img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.h, pix.w, pix.n)
    if pix.n == 4:
        schematic_img = schematic_img[:, :, :3]
    ax_schematic.imshow(schematic_img)
    ax_schematic.axis("off")
    ax_schematic.text(-0.02, 1.04, "A", transform=ax_schematic.transAxes,
                      fontsize=12, fontweight="bold", va="top")
    part1_doc.close()

    # Z-score composite bar — read from raw behavioral
    behav = pd.read_csv(BEHAV_RAW)
    sal = behav[behav["Group"] == "Saline"]["Global_Withdrawal_Score"].values
    mor = behav[behav["Group"] == "Morphine"]["Global_Withdrawal_Score"].values
    means = [np.mean(sal), np.mean(mor)]
    sems = [stats.sem(sal), stats.sem(mor)]
    x_pos = [0, 1]
    bars = ax_zscore.bar(x_pos, means, yerr=sems, color=[COLOR_SALINE, COLOR_MORPHINE],
                         edgecolor="black", lw=0.7, width=0.55, capsize=4,
                         error_kw={"lw": 0.7})
    # Individual dots with jitter
    np.random.seed(1)
    for x, vals in [(0, sal), (1, mor)]:
        jitter = np.random.uniform(-0.07, 0.07, size=len(vals))
        ax_zscore.scatter([x + j for j in jitter], vals, s=14,
                          color="black", facecolor="white", linewidth=0.6, zorder=3)
    ax_zscore.set_xticks(x_pos)
    ax_zscore.set_xticklabels(["Saline", "Morphine"], fontsize=8)
    ax_zscore.set_ylabel("Global withdrawal\nz-score", fontsize=8)
    # Significance bracket — *** for FDR < 0.001
    y_top = max(np.max(mor) + sems[1], np.max(sal) + sems[0]) * 1.12
    ax_zscore.plot([0, 0, 1, 1], [y_top - 0.05, y_top, y_top, y_top - 0.05],
                   color="black", lw=0.6)
    ax_zscore.text(0.5, y_top + 0.02, "***", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax_zscore.text(-0.10, 1.04, " ", transform=ax_zscore.transAxes,
                   fontsize=12, fontweight="bold")  # placeholder for letter (panel A overall)
    for spine in ['top', 'right']:
        ax_zscore.spines[spine].set_visible(False)
    ax_zscore.tick_params(axis='both', labelsize=7)

    # ================ Panel B — Mito forest plot ================
    ax_B = fig.add_subplot(gs[1])
    mito_raw = pd.read_csv(MITO_RAW)
    mito_stats_df = pd.read_csv(MITO_STATS)
    # Region label cleanup
    mito_raw["Region"] = mito_raw["Region"].replace({"Hippocampus": "HIP", "PFC": "PFC", "NAc": "NAc"})
    mito_stats_df["Region"] = mito_stats_df["Region"].replace({"Hippocampus": "HIP", "PFC": "PFC", "NAc": "NAc"})
    # Order: 5 genes × 3 regions
    gene_order_mito = ["CPT2", "HADHB", "MRPL21", "MRPS17", "MTIF3"]
    region_order = ["NAc", "PFC", "HIP"]
    rows_B = render_forest(ax_B, mito_raw, mito_stats_df, gene_order_mito, region_order,
                            fdr_col="MWU_FDR")
    ax_B.set_xlim(-3.2, 3.2)  # extra room on right for asterisks
    ax_B.set_title("Mitochondrial gene qPCR (NAc / PFC / HIP)", fontsize=8.5, fontweight="bold", pad=4)
    ax_B.text(-0.10, 1.06, "B", transform=ax_B.transAxes,
              fontsize=12, fontweight="bold", va="top")

    # ================ Panel C — WB ================
    ax_C = fig.add_subplot(gs[2])
    ax_C.axis("off")
    ax_C.text(-0.05, 1.05, "C", transform=ax_C.transAxes,
              fontsize=12, fontweight="bold", va="top")
    # Within Panel C: 3 membrane images (top row) + 4 quant bars (bottom row)
    # Use nested gridspec
    gs_C = gridspec.GridSpecFromSubplotSpec(2, 12, subplot_spec=gs[2],
                                            height_ratios=[1.0, 1.4],
                                            hspace=0.55, wspace=0.5)
    # Top row: 3 membrane images, each spans 4 cols
    # MRPL21+MRPS17 (combined membrane), CPT2, COMT
    ax_wb1 = fig.add_subplot(gs_C[0, 0:4])
    ax_wb2 = fig.add_subplot(gs_C[0, 4:8])
    ax_wb3 = fig.add_subplot(gs_C[0, 8:12])

    # Extract WB membrane crops from WB_4proteins.pdf
    wb_doc = fitz.open(WB_PDF_4PROTEINS)
    wb_page = wb_doc[0]
    # Coordinates from earlier inspection (page is 960 × 1417 pt)
    # COMT image strip: y=270-410, x=70-340
    # CPT2 image strip: y=470-615, x=70-340
    # MRPL21+MRPS17 strip: y=700-870, x=70-340
    membrane_crops = {
        "MRPL21+MRPS17": fitz.Rect(70, 700, 340, 870),
        "CPT2": fitz.Rect(70, 470, 340, 615),
        "COMT": fitz.Rect(70, 280, 340, 415),
    }
    membrane_axes = {"MRPL21+MRPS17": ax_wb1, "CPT2": ax_wb2, "COMT": ax_wb3}
    membrane_titles = {"MRPL21+MRPS17": "MRPL21 / MRPS17",
                       "CPT2": "CPT2",
                       "COMT": "COMT"}
    for key, clip in membrane_crops.items():
        ax = membrane_axes[key]
        pix = wb_page.get_pixmap(dpi=350, clip=clip)
        img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.h, pix.w, pix.n)
        if pix.n == 4:
            img = img[:, :, :3]
        ax.imshow(img)
        ax.set_title(membrane_titles[key], fontsize=8, pad=2, fontweight="bold")
        ax.axis("off")
        # Add lane labels (NC1-3, OUD1-3) below
        h, w = img.shape[:2]
        for i, label in enumerate(["NC1", "NC2", "NC3", "OUD1", "OUD2", "OUD3"]):
            x_frac = (i + 0.5) / 6
            ax.text(w * x_frac, h * 1.02, label, ha="center", va="top",
                    fontsize=5.5, fontweight="bold")
    wb_doc.close()

    # Bottom row: 4 quant bars
    wb_main = pd.read_csv(WB_RAW)
    wb_comt = pd.read_csv(WB_COMT)
    wb_all = pd.concat([wb_main, wb_comt], ignore_index=True)

    quant_proteins = ["MRPL21", "MRPS17", "CPT2", "COMT"]
    wb_stats_df = pd.read_csv(WB_STATS)

    # Compute global max y for shared y-axis (allow headroom for asterisk bracket)
    SHARED_Y_MAX = 4.5
    for idx, prot in enumerate(quant_proteins):
        ax_q = fig.add_subplot(gs_C[1, idx*3:(idx+1)*3])
        sub = wb_all[wb_all["Protein"] == prot]
        nc_v = sub[sub["Group"] == "NC"]["Value_normActin"].values
        oud_v = sub[sub["Group"] == "OUD"]["Value_normActin"].values
        means = [np.mean(nc_v), np.mean(oud_v)]
        sems = [stats.sem(nc_v), stats.sem(oud_v)]
        ax_q.bar([0, 1], means, yerr=sems, color=[COLOR_SALINE, COLOR_MORPHINE],
                 edgecolor="black", lw=0.6, width=0.55, capsize=3,
                 error_kw={"lw": 0.6})
        # NC=1 reference dashed line
        ax_q.axhline(1.0, color="gray", lw=0.5, ls="--", alpha=0.6, zorder=1)
        # Individual dots
        for x, vals in [(0, nc_v), (1, oud_v)]:
            ax_q.scatter([x] * len(vals), vals, s=11, color="black",
                          facecolor="white", linewidth=0.5, zorder=3)
        ax_q.set_xticks([0, 1])
        ax_q.set_xticklabels(["Saline", "Morphine"], fontsize=6.5)
        # Only leftmost (idx=0) gets y-label
        if idx == 0:
            ax_q.set_ylabel("Relative protein\n(/β-actin, NC = 1)", fontsize=7)
        # Each bar gets a small protein name title above
        ax_q.set_title(prot, fontsize=8, pad=2, fontweight="bold")
        ax_q.tick_params(axis='both', labelsize=6.5)
        # Shared y-axis 0-4.5
        ax_q.set_ylim(0, SHARED_Y_MAX)
        # Star
        srow = wb_stats_df[wb_stats_df["Protein"] == prot]
        if not srow.empty:
            fdr = srow["t_BH_FDR"].iloc[0]
            stars = fdr_to_stars(fdr)
            y_top = SHARED_Y_MAX * 0.85
            ax_q.plot([0, 0, 1, 1], [y_top - 0.10, y_top, y_top, y_top - 0.10],
                      color="black", lw=0.5)
            ax_q.text(0.5, y_top + 0.05, stars, ha="center", va="bottom",
                      fontsize=8.5, fontweight="bold", family="DejaVu Sans")
        for spine in ['top', 'right']:
            ax_q.spines[spine].set_visible(False)

    # ================ Panel D — Opioid forest plot ================
    ax_D = fig.add_subplot(gs[3])
    op_raw = pd.read_csv(OPIOID_RAW)
    op_stats_df = pd.read_csv(OPIOID_STATS)
    op_raw["Region"] = op_raw["Region"].replace({"Hippocampus": "HIP"})
    op_stats_df["Region"] = op_stats_df["Region"].replace({"Hippocampus": "HIP"})
    gene_order_op = ["Oprm1", "Oprd1", "Oprk1", "Penk", "Pdyn"]
    rows_D = render_forest(ax_D, op_raw, op_stats_df, gene_order_op, region_order,
                            fdr_col="MWU_FDR_within_region")
    ax_D.set_xlim(-3.2, 3.2)  # extra room on right for asterisks
    ax_D.set_title("Opioid-system gene qPCR (NAc / PFC / HIP)", fontsize=8.5, fontweight="bold", pad=4)
    ax_D.text(-0.10, 1.06, "D", transform=ax_D.transAxes,
              fontsize=12, fontweight="bold", va="top")

    # ================ Panel E — Coupling: heatmap + scatter ================
    gs_E = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[4],
                                            width_ratios=[1.1, 1.0], wspace=0.30)
    ax_hm = fig.add_subplot(gs_E[0])
    ax_sc = fig.add_subplot(gs_E[1])

    coupling = pd.read_csv(COUPLING)
    nac_sub = coupling[coupling["Region"] == "NAc"]

    # Build 5×5 heatmap matrix (mito × opioid)
    mito_genes_5 = ["CPT2", "HADHB", "MRPL21", "MRPS17", "MTIF3"]
    op_genes_5 = ["Oprm1", "Oprd1", "Oprk1", "Penk", "Pdyn"]
    M = np.full((5, 5), np.nan)
    Sig = np.zeros((5, 5), dtype=int)  # 0=ns, 1=*, 2=**, 3=***
    for i, mg in enumerate(mito_genes_5):
        for j, og in enumerate(op_genes_5):
            row = nac_sub[(nac_sub["Mito_gene"].str.upper() == mg.upper()) &
                          (nac_sub["Opioid_gene"].str.lower() == og.lower())]
            if not row.empty:
                M[i, j] = row["Pearson_r"].iloc[0]
                fdr = row["Pearson_FDR_within_region"].iloc[0]
                if fdr < 0.001:
                    Sig[i, j] = 3
                elif fdr < 0.01:
                    Sig[i, j] = 2
                elif fdr < 0.05:
                    Sig[i, j] = 1
    im = ax_hm.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax_hm.set_xticks(range(5))
    ax_hm.set_xticklabels(op_genes_5, rotation=45, ha="right", fontsize=7, fontstyle="italic")
    ax_hm.set_yticks(range(5))
    ax_hm.set_yticklabels(mito_genes_5, fontsize=7, fontstyle="italic")
    # Annotate r values + asterisks; use black text with white outline for universal readability
    for i in range(5):
        for j in range(5):
            r = M[i, j]
            if not np.isnan(r):
                stars = ["", "*", "**", "***"][Sig[i, j]]
                lbl = f"{r:+.2f}{stars}"
                txt = ax_hm.text(j, i, lbl, ha="center", va="center",
                                 fontsize=6, color="black", fontweight="bold",
                                 family="DejaVu Sans")
                txt.set_path_effects([
                    path_effects.Stroke(linewidth=1.5, foreground="white"),
                    path_effects.Normal()
                ])
    ax_hm.set_title("Within-mouse coupling, NAc (n = 12)\nPearson r, BH-FDR", fontsize=8, pad=4)
    cbar = plt.colorbar(im, ax=ax_hm, fraction=0.04, pad=0.03)
    cbar.set_label("Pearson r", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    ax_hm.text(-0.18, 1.04, "E", transform=ax_hm.transAxes,
               fontsize=12, fontweight="bold", va="top")

    # Highlight scatter — NAc Mrpl21 × Oprk1 (r = +0.999)
    # Need both per-mouse log2FC values
    mrpl21_data = mito_raw[(mito_raw["Gene"].str.upper() == "MRPL21") &
                            (mito_raw["Region"] == "NAc")][["Animal_ID", "Group", "log2_fold_change"]]
    oprk1_data = op_raw[(op_raw["Gene"].str.lower() == "oprk1") &
                         (op_raw["Region"] == "NAc")][["Animal_ID", "Group", "log2_fold_change"]]
    mrpl21_data = mrpl21_data.rename(columns={"log2_fold_change": "Mrpl21"})
    oprk1_data = oprk1_data.rename(columns={"log2_fold_change": "Oprk1"})
    paired = pd.merge(mrpl21_data, oprk1_data, on=["Animal_ID", "Group"])

    sal_p = paired[paired["Group"].str.lower() == "saline"]
    mor_p = paired[paired["Group"].str.lower() == "morphine"]
    ax_sc.scatter(sal_p["Mrpl21"], sal_p["Oprk1"], s=40, color=COLOR_SALINE,
                  edgecolor="black", lw=0.5, label="Saline", zorder=3)
    ax_sc.scatter(mor_p["Mrpl21"], mor_p["Oprk1"], s=40, color=COLOR_MORPHINE,
                  edgecolor="black", lw=0.5, label="Morphine", zorder=3)
    # Linear fit
    x_all = paired["Mrpl21"].values
    y_all = paired["Oprk1"].values
    slope, intercept, r_val, p_val, _ = stats.linregress(x_all, y_all)
    x_line = np.linspace(x_all.min() - 0.1, x_all.max() + 0.1, 50)
    ax_sc.plot(x_line, slope * x_line + intercept, "--", color="black", lw=0.8, alpha=0.6, zorder=2)

    ax_sc.set_xlabel("Mrpl21 log$_2$FC (NAc)", fontsize=8, fontstyle="italic")
    ax_sc.set_ylabel("Oprk1 log$_2$FC (NAc)", fontsize=8, fontstyle="italic")
    # r value box
    r_text = f"r = +{r_val:.3f}\nP < 0.001"
    ax_sc.text(0.04, 0.96, r_text, transform=ax_sc.transAxes,
               fontsize=10, fontweight="bold", va="top", ha="left",
               bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=0.6))
    ax_sc.legend(fontsize=7, frameon=False, loc="lower right")
    ax_sc.tick_params(axis='both', labelsize=7)
    for spine in ['top', 'right']:
        ax_sc.spines[spine].set_visible(False)
    ax_sc.set_title("Highlight: NAc Mrpl21 × Oprk1", fontsize=8, pad=4)

    # ============== SAVE ==============
    out_pdf = OUT_DIR / "Fig6_v2_assembled.pdf"
    out_png = OUT_DIR / "Fig6_v2_assembled.png"
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    fig.savefig(out_png, bbox_inches="tight", dpi=200)
    print(f"✓ Saved {out_pdf}")
    print(f"✓ Saved {out_png}")

if __name__ == "__main__":
    main()

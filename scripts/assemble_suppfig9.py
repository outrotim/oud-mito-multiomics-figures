"""
Supplementary Figure 9 — full per-region details for Fig 6.

Layout (8.5" × 11"):
  Panel A: 3 behavioral details (Jumps, Wet-Dog Shakes, %Body Weight Loss)
  Panel B: Mito qPCR per-region bars (3 regions × 5 genes = 15 bars)
  Panel C: Opioid qPCR per-region bars (3 regions × 5 genes = 15 bars)
  Panel D: PFC + HIP within-mouse coupling (2 heatmaps + 2 highlight scatter)

Style matched to Fig 6 v2: blue/orange Saline/Morphine, individual dots, MWU+BH-FDR asterisks.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import patheffects as path_effects
from scipy import stats

# ============== STYLE (match Fig 6 v2) ==============
COLOR_SALINE = "#4A90E2"
COLOR_MORPHINE = "#E67E22"
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 8
plt.rcParams["axes.linewidth"] = 0.6
plt.rcParams["xtick.major.width"] = 0.6
plt.rcParams["ytick.major.width"] = 0.6
plt.rcParams["xtick.major.size"] = 2.5
plt.rcParams["ytick.major.size"] = 2.5
plt.rcParams["pdf.fonttype"] = 42

ROOT = Path(os.environ.get("OUD_PROJECT_ROOT", Path(__file__).resolve().parent.parent))
DATA_DIR = ROOT / "data"  # user-supplied input CSVs go here
BEHAV_RAW = DATA_DIR / "behavioral_raw_data.csv"
BEHAV_STATS = DATA_DIR / "behavioral_8items_statistics.csv"
MITO_RAW = DATA_DIR / "mito_qPCR_log2FC_perSample.csv"
MITO_STATS = DATA_DIR / "mito_qPCR_statistics_MWU.csv"
OPIOID_RAW = DATA_DIR / "opioid_qPCR_log2FC_perSample.csv"
OPIOID_STATS = DATA_DIR / "opioid_qPCR_statistics_MWU.csv"
COUPLING = DATA_DIR / "WithinMouse_Mito_x_Opioid_Correlations.csv"
OUT_DIR = ROOT / "outputs"


def fdr_to_stars(fdr):
    if pd.isna(fdr): return ""
    if fdr < 0.001: return "***"
    if fdr < 0.01:  return "**"
    if fdr < 0.05:  return "*"
    return "ns"


def render_qpcr_bar(ax, raw, gene, region, fdr, *, gene_col="Gene", region_col="Region"):
    """One qPCR bar (Saline vs Morphine) with individual dots, error bar, asterisk."""
    sub = raw[(raw[gene_col].str.upper() == gene.upper()) & (raw[region_col] == region)]
    sal = sub[sub["Group"].str.lower() == "saline"]["log2_fold_change"].values
    mor = sub[sub["Group"].str.lower() == "morphine"]["log2_fold_change"].values
    means = [np.mean(sal), np.mean(mor)]
    sems = [stats.sem(sal), stats.sem(mor)]
    ax.bar([0, 1], means, yerr=sems, color=[COLOR_SALINE, COLOR_MORPHINE],
           edgecolor="black", lw=0.5, width=0.55, capsize=3, error_kw={"lw": 0.5})
    ax.axhline(0, color="black", lw=0.4, alpha=0.5)
    np.random.seed(7)
    for x, vals in [(0, sal), (1, mor)]:
        jitter = np.random.uniform(-0.06, 0.06, size=len(vals))
        ax.scatter([x + j for j in jitter], vals, s=8, color="black",
                   facecolor="white", linewidth=0.4, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["S", "M"], fontsize=6)
    ax.tick_params(axis='both', labelsize=6)
    # Asterisk above bars
    stars = fdr_to_stars(fdr)
    y_top = max(means[0] + sems[0], means[1] + sems[1], 0.1)
    y_bottom = min(means[0] - sems[0], means[1] - sems[1], 0)
    if stars:
        ax.text(0.5, y_top + 0.1, stars, ha="center", va="bottom",
                fontsize=7.5, fontweight="bold", family="DejaVu Sans")
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    return y_bottom, y_top


def panel_qpcr_grid(parent_gs, raw_df, stats_df, gene_order, regions, *,
                    fdr_col, gene_col="Gene", region_col="Region",
                    region_short_map=None):
    """
    Build a sub-grid of qPCR bars: 3 rows (regions) × 5 cols (genes) = 15 axes.
    Returns the list of axes.
    """
    inner = gridspec.GridSpecFromSubplotSpec(
        len(regions), len(gene_order),
        subplot_spec=parent_gs, wspace=0.35, hspace=0.55
    )
    axes = []
    fig = plt.gcf()
    # Compute global y-range across all bars for shared y-axis within panel
    raw_lo, raw_hi = 0, 0
    for region in regions:
        for gene in gene_order:
            sub = raw_df[(raw_df[gene_col].str.upper() == gene.upper()) & (raw_df[region_col] == region)]
            if not sub.empty:
                vals = sub["log2_fold_change"].values
                raw_lo = min(raw_lo, np.min(vals) - 0.2)
                raw_hi = max(raw_hi, np.max(vals) + 0.4)

    for ri, region in enumerate(regions):
        for gi, gene in enumerate(gene_order):
            ax = fig.add_subplot(inner[ri, gi])
            srow = stats_df[(stats_df[gene_col].str.upper() == gene.upper()) & (stats_df[region_col] == region)]
            fdr = srow[fdr_col].iloc[0] if not srow.empty else np.nan
            render_qpcr_bar(ax, raw_df, gene, region, fdr,
                           gene_col=gene_col, region_col=region_col)
            ax.set_ylim(raw_lo, raw_hi)
            # Top row: gene name (italic) as title
            if ri == 0:
                ax.set_title(f"$\\mathit{{{gene}}}$", fontsize=8, pad=2)
            # Leftmost col: region label as ylabel; y axis values
            if gi == 0:
                short = region_short_map.get(region, region) if region_short_map else region
                ax.set_ylabel(f"{short}\nlog$_2$FC", fontsize=7)
            else:
                ax.set_yticklabels([])
            axes.append(ax)
    return axes


def main():
    fig = plt.figure(figsize=(8.5, 11.0), dpi=150)

    # 4 rows: A=1.4 / B=3.4 / C=3.4 / D=2.6 ≈ 10.8"
    gs = gridspec.GridSpec(
        4, 1, figure=fig,
        height_ratios=[1.4, 3.4, 3.4, 2.6],
        hspace=0.85, top=0.96, bottom=0.04, left=0.08, right=0.97
    )

    # ============== Panel A — Behavioral details (3 bars) ==============
    behav = pd.read_csv(BEHAV_RAW)
    behav_stats = pd.read_csv(BEHAV_STATS)

    # Use the 3 detailed metrics (jumps, wet-dog, weight loss). Skip global z-score (in main Fig 6A).
    behav_items = [
        ("Jumps", "Jumps", "counts / 30 min"),
        ("Wet_Dog_Shakes", "Wet-dog shakes", "counts / 30 min"),
        ("Body_Weight_Loss_Percent", "% body-weight loss", "%"),
    ]
    gs_A = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0], wspace=0.45)
    for idx, (col, label, unit) in enumerate(behav_items):
        ax = fig.add_subplot(gs_A[idx])
        sal = behav[behav["Group"] == "Saline"][col].values
        mor = behav[behav["Group"] == "Morphine"][col].values
        means = [np.mean(sal), np.mean(mor)]
        sems = [stats.sem(sal), stats.sem(mor)]
        ax.bar([0, 1], means, yerr=sems, color=[COLOR_SALINE, COLOR_MORPHINE],
               edgecolor="black", lw=0.6, width=0.55, capsize=3, error_kw={"lw": 0.6})
        np.random.seed(idx + 1)
        for x, vals in [(0, sal), (1, mor)]:
            jitter = np.random.uniform(-0.07, 0.07, size=len(vals))
            ax.scatter([x + j for j in jitter], vals, s=12, color="black",
                       facecolor="white", linewidth=0.5, zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Saline", "Morphine"], fontsize=7)
        ax.set_ylabel(f"{label}\n({unit})", fontsize=7.5)
        ax.tick_params(axis='both', labelsize=6.5)
        # Stars
        srow = behav_stats[behav_stats["Variable_name"] == col]
        if not srow.empty:
            fdr = srow["MWU_BH_FDR"].iloc[0]
            stars = fdr_to_stars(fdr)
            y_top = max(means[0] + sems[0], means[1] + sems[1]) * 1.18
            ax.plot([0, 0, 1, 1], [y_top - y_top*0.03, y_top, y_top, y_top - y_top*0.03],
                    color="black", lw=0.5)
            ax.text(0.5, y_top * 1.02, stars, ha="center", va="bottom",
                    fontsize=9, fontweight="bold", family="DejaVu Sans")
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        if idx == 0:
            ax.text(-0.30, 1.15, "A", transform=ax.transAxes,
                    fontsize=12, fontweight="bold", va="top")

    # ============== Panel B — Mito per-region bars ==============
    mito_raw = pd.read_csv(MITO_RAW)
    mito_stats = pd.read_csv(MITO_STATS)
    mito_raw["Region"] = mito_raw["Region"].replace({"Hippocampus": "HIP"})
    mito_stats["Region"] = mito_stats["Region"].replace({"Hippocampus": "HIP"})
    region_order = ["NAc", "PFC", "HIP"]
    mito_genes = ["CPT2", "HADHB", "MRPL21", "MRPS17", "MTIF3"]
    region_short_map = {"NAc": "NAc", "PFC": "PFC", "HIP": "HIP"}

    # Wrap with header text via outer ghost axis
    ax_B_header = fig.add_subplot(gs[1], frameon=False)
    ax_B_header.set_xticks([]); ax_B_header.set_yticks([])
    # Panel B letter only; descriptive text → caption
    ax_B_header.text(-0.07, 1.20, "B", transform=ax_B_header.transAxes,
                     fontsize=12, fontweight="bold", va="top")
    panel_qpcr_grid(gs[1], mito_raw, mito_stats, mito_genes, region_order,
                   fdr_col="MWU_FDR", region_short_map=region_short_map)

    # ============== Panel C — Opioid per-region bars ==============
    op_raw = pd.read_csv(OPIOID_RAW)
    op_stats = pd.read_csv(OPIOID_STATS)
    op_raw["Region"] = op_raw["Region"].replace({"Hippocampus": "HIP"})
    op_stats["Region"] = op_stats["Region"].replace({"Hippocampus": "HIP"})
    op_genes = ["Oprm1", "Oprd1", "Oprk1", "Penk", "Pdyn"]

    ax_C_header = fig.add_subplot(gs[2], frameon=False)
    ax_C_header.set_xticks([]); ax_C_header.set_yticks([])
    # Panel C letter only; descriptive text → caption
    ax_C_header.text(-0.07, 1.20, "C", transform=ax_C_header.transAxes,
                     fontsize=12, fontweight="bold", va="top")
    panel_qpcr_grid(gs[2], op_raw, op_stats, op_genes, region_order,
                   fdr_col="MWU_FDR_within_region", region_short_map=region_short_map)

    # ============== Panel D — PFC and HIP within-mouse coupling ==============
    coupling = pd.read_csv(COUPLING)
    gs_D = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[3], wspace=0.35,
                                             width_ratios=[1, 1])
    mito_5 = ["CPT2", "HADHB", "MRPL21", "MRPS17", "MTIF3"]
    op_5 = ["Oprm1", "Oprd1", "Oprk1", "Penk", "Pdyn"]

    for di, region in enumerate(["PFC", "Hippocampus"]):
        ax_hm = fig.add_subplot(gs_D[di])
        sub = coupling[coupling["Region"] == region]
        M = np.full((5, 5), np.nan)
        Sig = np.zeros((5, 5), dtype=int)
        for i, mg in enumerate(mito_5):
            for j, og in enumerate(op_5):
                row = sub[(sub["Mito_gene"].str.upper() == mg.upper()) &
                          (sub["Opioid_gene"].str.lower() == og.lower())]
                if not row.empty:
                    M[i, j] = row["Pearson_r"].iloc[0]
                    fdr = row["Pearson_FDR_within_region"].iloc[0]
                    if fdr < 0.001: Sig[i, j] = 3
                    elif fdr < 0.01: Sig[i, j] = 2
                    elif fdr < 0.05: Sig[i, j] = 1
        im = ax_hm.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax_hm.set_xticks(range(5)); ax_hm.set_xticklabels(op_5, rotation=45, ha="right",
                                                          fontsize=7, fontstyle="italic")
        ax_hm.set_yticks(range(5)); ax_hm.set_yticklabels(mito_5, fontsize=7, fontstyle="italic")
        for i in range(5):
            for j in range(5):
                r = M[i, j]
                if not np.isnan(r):
                    stars = ["", "*", "**", "***"][Sig[i, j]]
                    txt = ax_hm.text(j, i, f"{r:+.2f}{stars}", ha="center", va="center",
                                     fontsize=6, color="black", fontweight="bold",
                                     family="DejaVu Sans")
                    txt.set_path_effects([
                        path_effects.Stroke(linewidth=1.5, foreground="white"),
                        path_effects.Normal()
                    ])
        region_label = "PFC" if region == "PFC" else "HIP"
        ax_hm.set_title(f"{region_label}: within-mouse coupling (n = 12)",
                        fontsize=8, pad=4)
        cbar = plt.colorbar(im, ax=ax_hm, fraction=0.04, pad=0.03)
        cbar.set_label("Pearson r", fontsize=7)
        cbar.ax.tick_params(labelsize=6)
        if di == 0:
            ax_hm.text(-0.20, 1.10, "D", transform=ax_hm.transAxes,
                       fontsize=12, fontweight="bold", va="top")

    # ============== Save ==============
    out_pdf = OUT_DIR / "SuppFig9_full_details.pdf"
    out_png = OUT_DIR / "SuppFig9_full_details.png"
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    fig.savefig(out_png, bbox_inches="tight", dpi=200)
    print(f"✓ Saved {out_pdf}")
    print(f"✓ Saved {out_png}")


if __name__ == "__main__":
    main()

"""
Fig 4 v2 — GTEx v10 brain co-expression (4 panels A-D, 8.5" × 9").

Panel A: Master heatmap, 13 mito × 11 opioid genes, mean Pearson r across 13 regions.
Panel B: Pathway schematic — two coupled axes (protective FAO ↔ μ/δ-OR; risk mito-translation ↔ KOR/PDYN/ABCB1) + COMT–dopamine bridge.
Panel C: Top-12 mito × opioid pairs by |mean r|, with cross-region SD bars (forest plot style).
Panel D: Per-region |r̄| comparison: FAO genes (CPT2/HADHB) vs mito-translation (MRPL21/MRPS17/MTIF3) vs other prioritized.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import patheffects as path_effects
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch

# ============== STYLE ==============
COLOR_PROTECTIVE = "#2266BB"   # blue — protective FAO axis
COLOR_RISK = "#CC4422"          # red — risk mito-translation axis
COLOR_BRIDGE = "#7B5BA6"        # purple — COMT/COMTD1 catecholamine bridge
COLOR_NEUTRAL = "#999999"

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 8
plt.rcParams["axes.linewidth"] = 0.6
plt.rcParams["xtick.major.width"] = 0.6
plt.rcParams["ytick.major.width"] = 0.6
plt.rcParams["xtick.major.size"] = 2.5
plt.rcParams["ytick.major.size"] = 2.5
plt.rcParams["pdf.fonttype"] = 42

# ============== DATA ==============
ROOT = Path(os.environ.get("OUD_PROJECT_ROOT", Path(__file__).resolve().parent.parent))
DATA_DIR = ROOT / "data"  # user-supplied input CSVs go here
GTEX = DATA_DIR / "Mito_x_Opioid_GTExV10_correlations.csv"
OUT_DIR = ROOT / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(GTEX)
print(f"Loaded {len(df)} pairs across {df['Region'].nunique()} regions")

# Mito and opioid gene order (consistent with manuscript narrative)
MITO_ORDER = ["CPT2", "HADHB",        # FAO axis
              "MRPL21", "MRPS17", "MTIF3",  # mito-translation axis
              "COMT", "COMTD1",       # catecholamine axis
              "QDPR", "DNA2", "MRPS7", "CHCHD2", "ERAL1", "MTHFD1L"]
OPIOID_ORDER = ["OPRM1", "OPRD1", "OPRK1", "PENK", "PDYN",  # core opioid
                "ABCB1",                                      # efflux
                "DRD2", "BDNF", "POMC", "CREB1", "SLC6A4"]   # related

# Group definitions for Panel D
FAO_GENES = ["CPT2", "HADHB"]
MITO_TX_GENES = ["MRPL21", "MRPS17", "MTIF3"]
CATECH_GENES = ["COMT", "COMTD1"]
OTHER_GENES = ["QDPR", "DNA2", "MRPS7", "CHCHD2", "ERAL1", "MTHFD1L"]

# Regions
REGIONS = sorted(df["Region"].unique())
REGION_SHORT = {
    "Amygdala": "AMY",
    "Anterior_cingulate_cortex_BA24": "ACC",
    "Caudate_basal_ganglia": "CAU",
    "Cerebellar_Hemisphere": "CBH",
    "Cerebellum": "CRBL",
    "Cortex": "CTX",
    "Frontal_Cortex_BA9": "FCx",
    "Hippocampus": "HIP",
    "Hypothalamus": "HYP",
    "Nucleus_accumbens_basal_ganglia": "NAc",
    "Putamen_basal_ganglia": "PUT",
    "Spinal_cord_cervical_c-1": "SC",
    "Substantia_nigra": "SN",
}


def fdr_to_stars(fdr):
    if pd.isna(fdr): return ""
    if fdr < 0.001: return "***"
    if fdr < 0.01:  return "**"
    if fdr < 0.05:  return "*"
    return ""


# ============== PANEL A: master heatmap ==============
def panel_A(ax):
    """Master heatmap: mito (rows) × opioid (cols); cell = mean r across 13 regions."""
    n_mito = len(MITO_ORDER)
    n_op = len(OPIOID_ORDER)
    M = np.full((n_mito, n_op), np.nan)
    Sig = np.zeros((n_mito, n_op), dtype=int)

    for i, mg in enumerate(MITO_ORDER):
        for j, og in enumerate(OPIOID_ORDER):
            sub = df[(df["Mito_gene"] == mg) & (df["Opioid_gene"] == og)]
            if not sub.empty:
                M[i, j] = sub["Pearson_r"].mean()
                # Significance: count regions BH-FDR < 0.05; use majority
                n_sig = (sub["Pearson_FDR_within_region"] < 0.05).sum()
                if n_sig >= len(sub) * 0.7:  # ≥70% regions sig
                    Sig[i, j] = 3
                elif n_sig >= len(sub) * 0.5:
                    Sig[i, j] = 2
                elif n_sig >= len(sub) * 0.3:
                    Sig[i, j] = 1

    im = ax.imshow(M, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(n_op))
    ax.set_xticklabels(OPIOID_ORDER, rotation=45, ha="right",
                       fontsize=7.5, fontstyle="italic")
    ax.set_yticks(range(n_mito))
    ax.set_yticklabels(MITO_ORDER, fontsize=7.5, fontstyle="italic")
    # Annotate r values with significance
    for i in range(n_mito):
        for j in range(n_op):
            r = M[i, j]
            if not np.isnan(r):
                stars = ["", "+", "*", "★"][Sig[i, j]]   # +: 30-50% sig regions, *: 50-70%, ★: ≥70%
                lbl = f"{r:+.2f}\n{stars}" if stars else f"{r:+.2f}"
                txt = ax.text(j, i, lbl, ha="center", va="center",
                              fontsize=5.5, color="black", fontweight="bold",
                              family="DejaVu Sans")
                txt.set_path_effects([
                    path_effects.Stroke(linewidth=1.2, foreground="white"),
                    path_effects.Normal()
                ])
    ax.set_title("Mean Pearson r across 13 GTEx v10 brain regions",
                 fontsize=9, pad=14, fontweight="bold")

    # Group annotation lines on left — moved to negative x with explicit clip_on=False
    for label, gene_list, ymin_idx, ymax_idx in [
        ("FAO", FAO_GENES, 0, 1),
        ("Mito-tx", MITO_TX_GENES, 2, 4),
        ("Catechol.", CATECH_GENES, 5, 6),
    ]:
        ann = ax.annotate('', xy=(-2.4, ymax_idx + 0.4), xytext=(-2.4, ymin_idx - 0.4),
                          arrowprops=dict(arrowstyle="-", color="black", lw=1.2),
                          xycoords='data')
        ann.set_annotation_clip(False)
        ax.text(-2.8, (ymin_idx + ymax_idx) / 2, label, rotation=90,
                ha="center", va="center", fontsize=6.5, fontweight="bold",
                clip_on=False)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Pearson r (mean over 13 regions)", fontsize=7)
    cbar.ax.tick_params(labelsize=6.5)


# ============== PANEL B: pathway schematic ==============
def panel_B(ax):
    """Pathway schematic — two coupled axes + COMT bridge."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

    # Box helper
    def box(x, y, w, h, label, color, txt_color="white", fontsize=7.5, fontweight="bold"):
        bb = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                            linewidth=0.8, edgecolor="black", facecolor=color)
        ax.add_patch(bb)
        ax.text(x + w/2, y + h/2, label, ha="center", va="center",
                fontsize=fontsize, fontweight=fontweight, color=txt_color)

    # Arrow helper
    def arrow(x1, y1, x2, y2, color="black", style="-|>", lw=1.3):
        a = FancyArrowPatch((x1, y1), (x2, y2),
                            arrowstyle=style, mutation_scale=10,
                            linewidth=lw, color=color, shrinkA=2, shrinkB=2)
        ax.add_patch(a)

    # ====== Left axis: PROTECTIVE (FAO) ======
    box(0.3, 7.5, 2.6, 1.0, "FAO ↑\nCPT2 · HADHB", COLOR_PROTECTIVE, fontsize=7.5)
    box(0.3, 4.2, 2.6, 1.0, "μ / δ-OR ↑\nOPRM1 · OPRD1", COLOR_PROTECTIVE, fontsize=7.5)
    arrow(1.6, 7.5, 1.6, 5.2, color=COLOR_PROTECTIVE, lw=1.5)
    ax.text(1.85, 6.35, "− r", color=COLOR_PROTECTIVE, fontsize=8, fontweight="bold",
            fontstyle="italic")
    box(0.3, 1.0, 2.6, 1.0, "↓ OUD risk", "#E8F4FF", txt_color="black",
        fontsize=8.5, fontweight="bold")
    arrow(1.6, 4.2, 1.6, 2.0, color=COLOR_PROTECTIVE, style="-|>", lw=1.5)

    # ====== Right axis: RISK (mito-translation) ======
    box(6.7, 7.3, 3.0, 1.3, "Mito-translation ↑\nMRPL21 · MRPS17\nMTIF3",
        COLOR_RISK, fontsize=7.0)
    box(7.0, 4.2, 2.6, 1.0, "KOR / PDYN ↑\nABCB1 ↑ efflux", COLOR_RISK, fontsize=7.5)
    arrow(8.2, 7.3, 8.2, 5.2, color=COLOR_RISK, lw=1.5)
    ax.text(8.45, 6.30, "+ r", color=COLOR_RISK, fontsize=8, fontweight="bold",
            fontstyle="italic")
    box(7.0, 1.0, 2.6, 1.0, "↑ OUD risk", "#FFEEEA", txt_color="black",
        fontsize=8.5, fontweight="bold")
    arrow(8.3, 4.2, 8.3, 2.0, color=COLOR_RISK, style="-|>", lw=1.5)

    # ====== Cross-axis inverse coupling ======
    arrow(2.9, 4.7, 7.0, 4.7, color="#888888", style="<->", lw=0.9)
    ax.text(4.95, 4.95, "cross-axis inverse coupling\n(e.g., CPT2 × PDYN)",
            ha="center", va="bottom", fontsize=6.5, fontstyle="italic", color="#666")

    # ====== COMT/COMTD1 bridge (middle) ======
    box(3.7, 8.4, 2.6, 0.9, "COMT · COMTD1\n(catecholamine bridge)",
        COLOR_BRIDGE, fontsize=6.5)
    arrow(5.0, 8.4, 5.0, 6.4, color=COLOR_BRIDGE, lw=1.0)
    box(3.7, 5.5, 2.6, 0.9, "DRD2 (substantia nigra)\nPCC ρ = +0.75",
        "#EFE7F7", txt_color="black", fontsize=6.5, fontweight="normal")

    ax.set_title("Two coupled regulatory axes in human brain (GTEx v10)",
                 fontsize=9, pad=4, fontweight="bold")


# ============== PANEL C: top-12 forest ==============
def panel_C(ax):
    """Top-12 strongest pairs by |mean r|."""
    # Compute mean and SD across regions per (mito, opioid)
    agg = (df.groupby(["Mito_gene", "Opioid_gene"])
              .agg(mean_r=("Pearson_r", "mean"),
                   sd_r=("Pearson_r", "std"),
                   n_sig=("Pearson_FDR_within_region",
                          lambda x: (x < 0.05).sum()),
                   n_regions=("Region", "count"))
              .reset_index())
    agg["abs_mean_r"] = agg["mean_r"].abs()
    top = agg.sort_values("abs_mean_r", ascending=False).head(12)
    top = top.sort_values("mean_r")  # bottom = most positive

    y = np.arange(len(top))
    SIG_X = 0.92  # x-position for sig.regions column
    for i, row in enumerate(top.itertuples()):
        color = COLOR_RISK if row.mean_r > 0 else COLOR_PROTECTIVE
        ax.errorbar(row.mean_r, i, xerr=row.sd_r, fmt='o', color=color,
                    ecolor=color, capsize=2, lw=1.0, markersize=5,
                    markeredgecolor="black", markeredgewidth=0.4)
        # frequency label (n significant regions)
        ax.text(SIG_X, i, f"{row.n_sig}/{row.n_regions}", va="center", ha="left",
                fontsize=6, color="#555", clip_on=False)

    labels = [f"$\\mathit{{{r.Mito_gene}}}$ × $\\mathit{{{r.Opioid_gene}}}$" for r in top.itertuples()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color="black", lw=0.4, ls="--", alpha=0.6)
    ax.set_xlabel("Mean Pearson r ± SD across 13 regions", fontsize=7.5)
    ax.set_xlim(-0.75, 0.85)
    ax.tick_params(axis='both', labelsize=7)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    # Column header for sig.regions, placed in axis-coords above plot area
    ax.text(SIG_X, len(top) - 0.4, "sig.\nregions", fontsize=6, color="#444",
            ha="left", va="bottom", fontstyle="italic", fontweight="bold",
            clip_on=False)
    ax.set_title("Top-12 strongest mito × opioid pairs", fontsize=9, pad=4, fontweight="bold")


# ============== PANEL D: per-region |r̄| bar ==============
def panel_D(ax):
    """Per-region average |r| broken down by gene group."""
    region_label_order = ["NAc", "PUT", "CAU", "FCx", "ACC", "CTX", "AMY", "HIP",
                         "HYP", "SN", "CRBL", "CBH", "SC"]
    region_full_for_label = {v: k for k, v in REGION_SHORT.items()}

    groups = {
        "FAO (CPT2, HADHB)": (FAO_GENES, COLOR_PROTECTIVE),
        "Mito-translation\n(MRPL21, MRPS17, MTIF3)": (MITO_TX_GENES, COLOR_RISK),
        "Other prioritized": (OTHER_GENES, COLOR_NEUTRAL),
    }

    # Compute |r̄| per group per region
    plot_data = {}  # group_name -> list of |r̄| in region order
    for gname, (gene_list, color) in groups.items():
        vals = []
        for r_short in region_label_order:
            r_full = region_full_for_label[r_short]
            sub = df[(df["Region"] == r_full) & (df["Mito_gene"].isin(gene_list))]
            mean_abs = sub["Pearson_r"].abs().mean() if not sub.empty else 0
            vals.append(mean_abs)
        plot_data[gname] = vals

    n_groups = len(groups)
    n_regions = len(region_label_order)
    x = np.arange(n_regions)
    width = 0.27

    legend_handles = []
    for gi, (gname, vals) in enumerate(plot_data.items()):
        color = list(groups.values())[gi][1]
        offset = (gi - 1) * width
        bars = ax.bar(x + offset, vals, width=width, color=color, label=gname,
                      edgecolor="black", lw=0.4)
        # Build proxy patches for clean legend (avoid bbox_inches truncation)
        legend_handles.append(Patch(facecolor=color, edgecolor="black", lw=0.4,
                                     label=gname))

    ax.set_xticks(x)
    ax.set_xticklabels(region_label_order, fontsize=6.5, rotation=45, ha="right")
    ax.set_ylabel("Mean |Pearson r| across opioid-system genes", fontsize=7.5)
    ax.tick_params(axis='both', labelsize=6.5)
    # Headroom: pad y-limit so the legend placed on top doesn't overlap bars
    ymax_bar = max(max(v) for v in plot_data.values())
    ax.set_ylim(0, ymax_bar * 1.50)
    # NAc-relevant region group: shade x range 0-2 (NAc, PUT, CAU)
    ax.axvspan(-0.5, 2.5, alpha=0.10, color="#FFE0B0", zorder=0)
    ax.text(1, ymax_bar * 1.05, "OUD-relevant\nstriatum", fontsize=6.5,
            ha="center", va="bottom", fontstyle="italic", color="#996633")

    # Legend INSIDE plot at upper-right, with white frame to avoid bar overlap
    ax.legend(handles=legend_handles, fontsize=6, frameon=True, framealpha=0.95,
              loc="upper right", bbox_to_anchor=(1.0, 1.0),
              ncol=1, handlelength=1.2, handletextpad=0.5,
              edgecolor="black")
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.set_title("Region-stratified co-expression strength by gene group",
                 fontsize=9, pad=4, fontweight="bold")


# ============== MAIN ==============
def main():
    fig = plt.figure(figsize=(8.5, 9.5), dpi=150)
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           height_ratios=[1.05, 1.0],
                           width_ratios=[1.0, 1.0],
                           hspace=0.55, wspace=0.45,
                           top=0.94, bottom=0.07, left=0.13, right=0.95)

    ax_A = fig.add_subplot(gs[0, 0])
    ax_B = fig.add_subplot(gs[0, 1])
    ax_C = fig.add_subplot(gs[1, 0])
    ax_D = fig.add_subplot(gs[1, 1])

    panel_A(ax_A)
    panel_B(ax_B)
    panel_C(ax_C)
    panel_D(ax_D)

    # Panel letters
    for ax, letter in [(ax_A, "A"), (ax_B, "B"), (ax_C, "C"), (ax_D, "D")]:
        ax.text(-0.12, 1.05, letter, transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="top")

    out_pdf = OUT_DIR / "Fig4_v2_assembled.pdf"
    out_png = OUT_DIR / "Fig4_v2_assembled.png"
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    fig.savefig(out_png, bbox_inches="tight", dpi=200)
    print(f"✓ Saved {out_pdf}")
    print(f"✓ Saved {out_png}")


if __name__ == "__main__":
    main()

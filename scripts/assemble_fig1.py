"""
Fig 1 v2 — Study design (8.5" × 8") workflow funnel + data-source strip.

Panel A: vertical workflow funnel showing 8 stages from 1,136 mito genes
         to in-vivo r = 0.999 coupling. Blue = in-silico (stages 1-5),
         orange = in-vivo (stages 6-8).

Panel B: data-source citation strip with all 5 datasets used (responds to R1.M4).
"""
import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

# ============== STYLE ==============
COLOR_INSILICO = "#3F7BB6"   # blue
COLOR_INVIVO = "#E87A2C"      # orange
COLOR_BRIDGE = "#7B5BA6"      # purple (cross-validation)
COLOR_HIGHLIGHT = "#C0392B"   # red highlight for r=0.999

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 8
plt.rcParams["pdf.fonttype"] = 42

ROOT = Path(os.environ.get("OUD_PROJECT_ROOT", Path(__file__).resolve().parent.parent))
OUT_DIR = ROOT / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def add_box(ax, x, y, w, h, lines, fill, *, txt_color="white",
            fontsize_main=9, fontsize_sub=7.5, edge="black", lw=0.8):
    """
    Place a rounded box with a main title line + optional sub-text lines.
    `lines` is either a single string or a list:
        ["Main title", "sub line 1", "sub line 2"]
    """
    bb = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06",
                        linewidth=lw, edgecolor=edge, facecolor=fill)
    ax.add_patch(bb)
    if isinstance(lines, str):
        lines = [lines]
    n = len(lines)
    if n == 1:
        ax.text(x + w/2, y + h/2, lines[0], ha="center", va="center",
                fontsize=fontsize_main, fontweight="bold", color=txt_color)
    else:
        ax.text(x + w/2, y + h - h*0.30, lines[0], ha="center", va="center",
                fontsize=fontsize_main, fontweight="bold", color=txt_color)
        for i, line in enumerate(lines[1:]):
            ax.text(x + w/2, y + h - h*(0.55 + i*0.25), line,
                    ha="center", va="center",
                    fontsize=fontsize_sub, color=txt_color)


def add_arrow(ax, x1, y1, x2, y2, *, color="black", lw=1.5, style="-|>"):
    a = FancyArrowPatch((x1, y1), (x2, y2),
                        arrowstyle=style, mutation_scale=14,
                        linewidth=lw, color=color, shrinkA=4, shrinkB=4)
    ax.add_patch(a)


def panel_A(ax):
    """Vertical 8-stage workflow funnel."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 16)
    ax.axis("off")

    # ===== Stage 1 — Discovery (widest) =====
    add_box(ax, 1.5, 14.5, 7.0, 1.1,
            ["1,136 mitochondrial protein-coding genes",
             "MitoCarta 3.0"],
            COLOR_INSILICO)

    # arrow to stage 2
    add_arrow(ax, 5.0, 14.5, 5.0, 13.7, color=COLOR_INSILICO)

    # ===== Stage 2 — Multi-omic instruments (3 sub-boxes) =====
    sub_w = 2.8
    sub_y = 12.4
    sub_h = 1.1
    add_box(ax, 0.5, sub_y, sub_w, sub_h,
            ["mQTL", "n = 1,980 indiv.", "→ 2,753 mQTLs (McRae 2018)"],
            COLOR_INSILICO, fontsize_main=8.5, fontsize_sub=6.0)
    add_box(ax, 3.6, sub_y, sub_w, sub_h,
            ["eQTL", "n = 31,684 indiv.", "→ 910 cis-eQTLs (Võsa 2021)"],
            COLOR_INSILICO, fontsize_main=8.5, fontsize_sub=6.0)
    add_box(ax, 6.7, sub_y, sub_w, sub_h,
            ["pQTL", "n = 54,219 indiv.", "→ 114 cis-pQTLs (Sun 2018)"],
            COLOR_INSILICO, fontsize_main=8.5, fontsize_sub=6.0)
    # 3 arrows merging into stage 3
    for x_start in [1.9, 5.0, 8.1]:
        add_arrow(ax, x_start, sub_y, 5.0, 11.6, color=COLOR_INSILICO, lw=1.0)

    # ===== Stage 3 — SMR + HEIDI + Coloc =====
    add_box(ax, 1.5, 10.4, 7.0, 1.2,
            ["SMR + HEIDI + Bayesian colocalization",
             "FinnGen Release 11 GWAS  ·  1,312 cases  /  445,828 controls"],
            COLOR_INSILICO, fontsize_main=9, fontsize_sub=7)
    add_arrow(ax, 5.0, 10.4, 5.0, 9.6, color=COLOR_INSILICO)

    # ===== Stage 4 — 13 prioritized genes =====
    add_box(ax, 2.0, 8.3, 6.0, 1.3,
            ["13 prioritized mitochondrial genes",
             "Tier 1 (n=5): CPT2 · MRPL21 · MRPS17 · MTIF3 · QDPR",
             "Tier 2/3 (n=8): COMT · COMTD1 · HADHB · …"],
            COLOR_INSILICO, fontsize_main=9, fontsize_sub=6.5)
    add_arrow(ax, 5.0, 8.3, 5.0, 7.5, color=COLOR_INSILICO)

    # ===== Stage 5 — Cross-validation (purple bridge) =====
    add_box(ax, 2.5, 6.3, 5.0, 1.2,
            ["In-silico cross-validation",
             "Deak et al. 2022 OUD GWAS",
             "MAGMA gene-set P = 0.011"],
            COLOR_BRIDGE, fontsize_main=9, fontsize_sub=7)
    add_arrow(ax, 5.0, 6.3, 5.0, 5.5, color=COLOR_BRIDGE)

    # ===== Stage 6 — GTEx co-expression =====
    add_box(ax, 1.5, 4.3, 7.0, 1.2,
            ["GTEx v10 brain co-expression",
             "13 brain regions  ·  Pearson + Spearman + MI  ·  BH-FDR"],
            COLOR_INVIVO, fontsize_main=9, fontsize_sub=7)
    add_arrow(ax, 5.0, 4.3, 5.0, 3.5, color=COLOR_INVIVO)

    # ===== Stage 7 — Mouse withdrawal =====
    add_box(ax, 2.0, 2.3, 6.0, 1.2,
            ["Mouse morphine / naloxone withdrawal",
             "C57BL/6J  ·  n = 6 / group  ·  qPCR + Western blot"],
            COLOR_INVIVO, fontsize_main=9, fontsize_sub=7)
    add_arrow(ax, 5.0, 2.3, 5.0, 1.6, color=COLOR_HIGHLIGHT, lw=2.0)

    # ===== Stage 8 — In-vivo coupling (HIGHLIGHT) =====
    add_box(ax, 2.5, 0.4, 5.0, 1.2,
            ["In-vivo mito × opioid coordinate regulation",
             "NAc Mrpl21 × Oprk1 :  r = +0.999 ★"],
            COLOR_HIGHLIGHT, fontsize_main=9, fontsize_sub=7.5,
            edge="black", lw=1.4)

    # ===== Side annotations: "in silico" vs "in vivo" =====
    ax.annotate('', xy=(-0.25, 14.0), xytext=(-0.25, 6.3),
                arrowprops=dict(arrowstyle="-", color=COLOR_INSILICO, lw=2.5),
                annotation_clip=False)
    ax.text(-0.55, 10.2, "in silico", rotation=90, ha="center", va="center",
            fontsize=9, color=COLOR_INSILICO, fontweight="bold",
            fontstyle="italic", clip_on=False)

    ax.annotate('', xy=(-0.25, 5.5), xytext=(-0.25, 0.4),
                arrowprops=dict(arrowstyle="-", color=COLOR_INVIVO, lw=2.5),
                annotation_clip=False)
    ax.text(-0.55, 2.95, "in vivo", rotation=90, ha="center", va="center",
            fontsize=9, color=COLOR_INVIVO, fontweight="bold",
            fontstyle="italic", clip_on=False)


def panel_B(ax):
    """Data sources strip (Reviewer 1.M4 response)."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 2)
    ax.axis("off")

    items = [
        ("FinnGen R11", "ICD-10 opioid-related disorder GWAS", "(1,312 cases / 445,828 ctrls)"),
        ("MitoCarta 3.0", "mitochondrial gene catalogue", "(1,136 genes)"),
        ("McRae 2018", "blood mQTL", "(n = 1,980; 2,753 mQTLs)"),
        ("Võsa 2021", "blood eQTL (eQTLGen)", "(n = 31,684; 910 cis-eQTLs)"),
        ("Sun 2018", "plasma pQTL (UK Biobank)", "(n = 54,219; 114 cis-pQTLs)"),
        ("GTEx v10", "brain co-expression", "(13 regions)"),
        ("Deak 2022", "European OUD GWAS", "(MAGMA replication)"),
    ]

    # Header
    ax.text(0.05, 1.85, "Data sources used:", fontsize=8.5, fontweight="bold",
            ha="left", va="top")

    # Layout: 7 items in 2 rows × 4 cols (with 1 empty)
    n_cols = 4
    col_w = 10 / n_cols
    for i, (cite, desc, count) in enumerate(items):
        col = i % n_cols
        row = i // n_cols
        x = col * col_w + 0.1
        y = 1.5 - row * 0.75
        ax.text(x, y, cite, fontsize=7, fontweight="bold", color=COLOR_INSILICO)
        ax.text(x, y - 0.18, desc, fontsize=6.5, color="#444")
        ax.text(x, y - 0.34, count, fontsize=6, color="#666", fontstyle="italic")

    # Bottom note: abbreviation definitions
    ax.text(0.05, 0.03,
            "mQTL: methylation quantitative trait loci · eQTL: expression QTL · "
            "pQTL: protein QTL · SMR: summary-data-based Mendelian randomization · "
            "HEIDI: heterogeneity in dependent instruments · BH-FDR: Benjamini–Hochberg false-discovery rate",
            fontsize=5.5, color="#444", style="italic", ha="left", va="bottom")


def main():
    fig = plt.figure(figsize=(8.5, 8.5), dpi=150)
    gs = gridspec.GridSpec(2, 1, figure=fig,
                           height_ratios=[7.0, 1.5],
                           hspace=0.18,
                           top=0.96, bottom=0.04, left=0.06, right=0.96)
    ax_A = fig.add_subplot(gs[0])
    ax_B = fig.add_subplot(gs[1])

    panel_A(ax_A)
    panel_B(ax_B)

    # Panel letters
    ax_A.text(-0.04, 1.00, "A", transform=ax_A.transAxes,
              fontsize=14, fontweight="bold", va="top")
    ax_B.text(-0.04, 1.05, "B", transform=ax_B.transAxes,
              fontsize=14, fontweight="bold", va="top")

    out_pdf = OUT_DIR / "Fig1_v2_assembled.pdf"
    out_png = OUT_DIR / "Fig1_v2_assembled.png"
    fig.savefig(out_pdf, bbox_inches="tight", dpi=300)
    fig.savefig(out_png, bbox_inches="tight", dpi=200)
    print(f"✓ Saved {out_pdf}")
    print(f"✓ Saved {out_png}")


if __name__ == "__main__":
    main()

# Mitochondrial Gene Networks in Opioid Use Disorder — Figure-Assembly Scripts

Public-release figure-assembly code accompanying the manuscript:

> **Mitochondrial Gene Networks in Opioid Use Disorder: Multi-omic Evidence for Pathway-Specific Risk and Resilience.** *(Submitted to Translational Psychiatry, 2025/2026.)*

These scripts reproduce **Figure 1**, **Figure 4**, **Figure 6**, and **Supplementary Figure 9** from the manuscript using user-supplied summary-statistics CSV files.

---

## Repository contents

```
.
├── README.md                 # this file
├── LICENSE                   # MIT
├── requirements.txt          # Python deps
└── scripts/
    ├── assemble_fig1.py      # Fig 1 — workflow funnel (no data input; pure schematic)
    ├── assemble_fig4.py      # Fig 4 — GTEx v10 mito × opioid co-expression (4 panels A–D)
    ├── assemble_fig6.py      # Fig 6 — mouse withdrawal model: forest plots + WB + within-mouse coupling (5 panels A–E)
    └── assemble_suppfig9.py  # Supp Fig 9 — full per-region bar charts + PFC/HIP coupling heatmaps
```

---

## Data availability — layered statement

This repository contains **only figure-assembly code**. Input data are obtained from the following sources, with different access conditions:

### 1. Publicly available reference datasets (open access)

The following data sources are publicly available from their respective consortia and are consumed (or summarised) by the scripts:

| Dataset | Source | Access |
|---|---|---|
| FinnGen R11 GWAS for ICD-10 endpoint F11 | https://www.finngen.fi | Public |
| MitoCarta 3.0 mitochondrial gene catalogue | https://www.broadinstitute.org/mitocarta | Public |
| Blood mQTL (McRae et al., 2018) | Genome Med 10:39 | Public via supplementary materials |
| Blood eQTL (Võsa et al., 2021 / eQTLGen) | https://www.eqtlgen.org/ | Public |
| Plasma pQTL (Sun et al., 2018; UK Biobank proteomics) | Open Targets Genetics / UK Biobank | Public via Open Targets |
| GTEx v10 brain expression (13 regions) | https://gtexportal.org/ | Public |
| Allen Human Brain Atlas microarray | https://human.brain-map.org/ | Public |
| Deak et al. (2022) European OUD GWAS | Mol Psychiatry 27:3970–3979 | Public via supplementary materials |

### 2. Mouse experimental data (released as journal supplementary materials)

Per-animal qPCR `log2` fold-change values, Western-blot quantifications and behavioural scores are not deposited in this repository. They are provided as supplementary CSV files alongside the journal publication. Download these and place them in a local `./data/` directory using the file names listed under "Expected input files" below.

### 3. Restricted / not publicly available

- **Western-blot raw membrane TIFF images** are retained by the authors and available upon reasonable request, in keeping with institutional policy.

---

## Expected input files (place in `./data/` or set `OUD_PROJECT_ROOT`)

`assemble_fig1.py` — **none** (pure schematic, no data input).

`assemble_fig4.py`:
- `Mito_x_Opioid_GTExV10_correlations.csv` — derived from public GTEx v10 normalised expression; 1,846 region × pair correlations (released as journal Supplementary Table 7 source CSV)

`assemble_fig6.py`:
- `mito_qPCR_log2FC_perSample.csv`
- `mito_qPCR_statistics_MWU.csv`
- `opioid_qPCR_log2FC_perSample.csv`
- `opioid_qPCR_statistics_MWU.csv`
- `behavioral_raw_data.csv`
- `WB_perSample.csv`
- `WB_COMT_perSample.csv`
- `WB_statistics.csv`
- `WB_membranes.pdf` (representative membrane images, journal supplementary)
- `Fig6_schematic_template.pdf` (mouse-dosing schematic, journal supplementary)
- `WithinMouse_Mito_x_Opioid_Correlations.csv`

`assemble_suppfig9.py`: same as `assemble_fig6.py` except `WB_membranes.pdf` and `Fig6_schematic_template.pdf` are not required.

All CSVs use a tidy long format with columns `Sample_ID`, `Region`, `Gene`, `Group`, `Animal_ID`, `log2_fold_change` (or analogous for behavioural / WB data). The expected schema for each file is documented in the journal supplementary materials.

---

## Usage

```bash
# 1. Install Python dependencies
pip install -r requirements.txt

# 2. Place input CSVs into the ./data directory
mkdir -p data
# (download the journal supplementary CSVs and place them here)

# 3. Run any figure script
python scripts/assemble_fig1.py    # Fig 1 — no data needed
python scripts/assemble_fig4.py    # Fig 4
python scripts/assemble_fig6.py    # Fig 6
python scripts/assemble_suppfig9.py  # Supp Fig 9

# Outputs go to ./outputs/
ls outputs/
# Fig1_v2_assembled.pdf, Fig1_v2_assembled.png, ...
```

You can override the project root by setting an environment variable:

```bash
export OUD_PROJECT_ROOT=/path/to/your/project
python scripts/assemble_fig4.py
```

---

## Caveats and reuse guidelines

- These scripts are designed to **reproduce the manuscript figures from summary-statistics CSV files**, not to re-run the underlying SMR + colocalization pipeline. The original SMR / HEIDI / coloc pipeline used standard tools (SMR v1.03, R `coloc` v5.1.0, MAGMA) on publicly available QTL data — please refer to the manuscript Methods for the upstream pipeline.
- Human GTEx v10 co-expression results are **steady-state** observations across donors and may not generalise to acute pharmacological perturbations.
- Mouse data were collected from male C57BL/6J only (n = 6 per group; n = 3 per group for WB). The corresponding analyses are not stratified by sex.
- The `assemble_fig6.py` workflow assumes the `Fig6_schematic_template.pdf` placeholder for the mouse-dosing diagram; if not provided, the schematic panel renders blank.

---

## License

Code in this repository is released under the **MIT License** (see `LICENSE`).

Derived numerical data (when included) are released under **CC BY 4.0** unless otherwise noted.

---

## Citation

> *(Citation to be updated upon publication; please refer to the published manuscript for the canonical reference.)*

If you use this code or the methods described in the manuscript, please cite both the article and this repository.

---

## Issues and contributions

Bug reports and reproducibility issues can be opened as GitHub issues. For questions about the underlying biology or data availability, please contact the corresponding authors as listed in the published manuscript.

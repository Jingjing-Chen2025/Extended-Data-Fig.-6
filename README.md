# Extended Data Fig. 6 — ESR1 at the ligand promoters, and tamoxifen in vivo

Code for Extended Data Fig. 6 of *A neurotrophic secretome from adipose-progenitor-like stromal cells promotes
prostate cancer lineage plasticity*. Seven panels supporting Fig. 5c–e.

| Panel | What it shows | Code here? |
|---|---|---|
| **a** | Transcription-factor occupancy at the *MDK* and *TNFSF12* promoters, over the full ReMap TF universe | yes |
| **b** | H3K27ac / H3K4me1 density across each promoter window, with the state factors' binding sites beneath | yes |
| **c** | sgESR1 knockdown efficiency | no — bench assay |
| **d** | Tamoxifen suppresses tumour growth | no — endpoint photograph |
| **e–g** | SomaScan plasma proteomics | partial, see below |

## Layout and order

```
R/common.R             paths (environment variables), helpers, palette
R/01_remap_accessible.R  ReMap peaks in the two promoter windows, intersected with adipose-lineage histone marks
R/02_panels_ab.R         draws panels a and b as one block
R/03_assemble_base.R     draws the title / panel letters / legend layer
python/compose.py        clips the artwork and merges it onto that layer, losslessly
reference/               the SomaScan analysis as held in the manuscript folder (see below)
data/state_TFs_fig1_2.csv  the 138 PI16+/POSTN+ state transcription factors used to label panels a and b
```

```bash
Rscript R/01_remap_accessible.R
Rscript R/02_panels_ab.R
Rscript R/03_assemble_base.R
ED6_BENCH="<path to the bench artwork PDF>" python python/compose.py   # -> figures/Extended Data Fig. 6.pdf
```

The output page is **612 × 960 pt (8.5 × 13.3 in)**. That is taller than the 180 × 240 mm a single Extended Data
figure is normally allowed. Panels a and b are intrinsically wide and tall, and shrinking them far enough to fit
seven panels on a shorter page made the type unreadable, so legibility was chosen over paper size. If the limit
has to be met, the natural cut is moving e–g into their own figure.

## What panels a and b are, and what they replaced

The previously published version of this analysis ranked 36 hand-picked transcription factors by a composite
"binding intensity". That quantity cannot be computed from the file it read: `remap2022_nr_macs2_hg38` is BED9
with **no signalValue**, so what was read as signal is column 7, `thickStart` — a genomic coordinate, which
correlates with peak start at r = 1.0000. No binding-intensity estimate is reported anywhere in this figure.

What the data do support is reported on two axes kept separate, because combining them is what made the earlier
metric opaque:

- **how many** distinct ReMap peaks a factor has in the −2,000/+500 window, and
- **how strongly marked** the chromatin around its best site is (Roadmap histone reads per kb).

A peak counts only if it falls in chromatin carrying an active mark in an adipose-lineage epigenome: **E063**
subcutaneous adipose, **E065** visceral adipose, **E025** adipose-derived MSC, marks H3K27ac and H3K4me1. E025
H3K27ac does not exist in the Roadmap release used — that epigenome ships H3K9ac instead — so ADMSC contributes
H3K4me1 only. This is stated rather than silently substituted.

Accessibility is sparse by construction: a 3.5 kb promoter receives few genome-wide reads (475 across all five
tracks), so it is reported as a rank, never as a quantitative signal.

The two *TNFSF12* promoter windows are plotted separately rather than as one span, because treating them as a
single interval would place the intergenic gap between them inside the plotted window.

## Panels e–g: the SomaScan code is incomplete

`reference/somascan_from_SomaScan.docx.R` is the analysis as it is stored in
`Nature Manuscript/Extended Data Fig. 5/SomaScan.docx`. **That Word document holds only part of the script** — it
stops mid-line inside the recovery-factor diagnostic plot, so the per-Source panel loop that writes the figure
PDF is missing. Nothing has been added or guessed. The file is kept as the record of the method, is not a build
step, and does not parse. Two lines were changed and only two: the absolute paths to the ADAT file and to the
aptamer-to-Source spreadsheet now read from `SOMASCAN_ADAT` and `SOMASCAN_MAP`.

The three panels are pages 4, 8 and 15 of `Source_comparison_RECOVERED_Normal_Vehicle_Tamoxifen.pdf`:
**e** Source ASPC, recovered absolute scale; **f** Source NEPC, recovered absolute scale; **g** Source PCa,
compositional (CLR). They were identified by matching axis ranges and *P* values against the per-Source pages,
because the panels in the artwork carry **no Source label of their own** — worth confirming against the original,
since those labels are now printed on the figure.

## Paths

No path is hard-coded. `R/common.R` reads `ED6_REPO`, `ED6_DATA` (working directory), `ED6_REMAP` (the ReMap 2022
non-redundant BED, bgzipped and tabix-indexed), `ED6_HIST` (a directory of Roadmap `<epigenome>_<mark>.bed`
files), and `ED6_OUT`. `python/compose.py` reads `ED6_REPO`, `ED6_OUT` and `ED6_BENCH`.

## Inputs that are NOT in this repository

- **ReMap 2022** non-redundant peaks (hg38) and the **Roadmap** tagAlign-derived BED files for E063, E065 and E025
  — both large public downloads.
- **The bench artwork** for panels c, d and e–g. `python/compose.py` clips these from an existing 612 × 792 pt
  page named by `ED6_BENCH`; the code that drew them is not in this repository, and for c and d there is none.
- **The upstream single-cell analysis.** `data/state_TFs_fig1_2.csv` is the finished product of the Figure 1–2
  state analysis and is read rather than recomputed; that chain lives in
  [Figure-5-new](https://github.com/Jingjing-Chen2025/Figure-5-new).

## Requirements

R >= 4.4 with `data.table`, `dplyr`, `tidyr`, `ggplot2`, `ggrepel`, `patchwork`, `GenomicRanges`, `Rsamtools`,
`TxDb.Hsapiens.UCSC.hg38.knownGene`, `org.Hs.eg.db`, `AnnotationDbi`.
Python >= 3.10 with `pypdf` (`pip install -r python/requirements.txt`).

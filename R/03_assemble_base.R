# EXTENDED DATA FIG. 6 - the title, panel letters and legend layer.
#
# The finished figure is this text layer with five pieces of artwork merged onto it by python/compose.py.
# Panels a and b come from R/02_panels_ab.R as one block, split by the compose step; c, d and e-g are clipped
# out of the existing bench artwork, which has no code in this repository (see README).
#
# Usage: Rscript R/03_assemble_base.R
suppressPackageStartupMessages({library(grid)})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
pt  <- function(x) unit(x, "bigpts")
txt <- function(s, x, y, size = 7, face = "plain", col = MUTED, hjust = "left", lineheight = 1.25)
  grid.text(s, x = pt(x), y = pt(y), just = c(hjust, "bottom"),
            gp = gpar(fontfamily = "Helvetica", fontsize = size, fontface = face, col = col, lineheight = lineheight))

f <- file.path(ED6_OUT, "ED6 base.pdf")
cairo_pdf(f, width = 612 / 72, height = 960 / 72, onefile = FALSE); grid.newpage()
txt("Extended Data Fig. 6 | ESR1 at the ligand promoters, and tamoxifen in vivo", 30, 936, 9.5, "bold", INK)
txt("Support for Fig. 5c-e: where the state transcription factors bind, that ESR1 was knocked down, and what tamoxifen does to tumours and to the circulating proteome.",
    30, 925, 6.6)
# a and b carry their own letters, drawn by R/02_panels_ab.R
txt("c", 30, 400, 8.5, "bold", INK);  txt("sgESR1 knockdown efficiency", 40, 400, 8, "bold", INK)
txt("d", 315, 400, 8.5, "bold", INK); txt("Tamoxifen suppresses tumour growth", 325, 400, 8, "bold", INK)
txt("e", 30, 213, 8.5, "bold", INK);  txt("Plasma, ASPC-source aptamers", 40, 213, 6.8, "bold", INK)
txt("f", 214, 213, 8.5, "bold", INK); txt("Plasma, NEPC-source aptamers", 224, 213, 6.8, "bold", INK)
txt("g", 398, 213, 8.5, "bold", INK); txt("Composition (CLR), PCa-source", 408, 213, 6.8, "bold", INK)
txt(paste(
  "a, Every transcription factor with at least one ReMap 2022 non-redundant peak in the -2,000/+500 promoter window of each ligand (584 at MDK, 396 at TNFSF12). Horizontal axis, number of distinct",
  "peaks; vertical axis, Roadmap histone read density over that factor's most strongly marked peak. PI16+ and POSTN+ mark the factors nominated for each stromal state. Axes are kept separate between",
  "promoters: the two windows are not contiguous. b, Upper, density of H3K27ac and H3K4me1 reads across each window in three adipose-lineage Roadmap epigenomes; curves are densities, so they show",
  "where reads concentrate, not coverage. Coverage is near-complete, so 99.7% of ReMap peaks sit in marked chromatin and restricting to marked chromatin retains every factor. Lower, the individual",
  "binding sites of the state-associated factors. The ReMap non-redundant file carries no signal values, so no binding-intensity estimate is reported anywhere in this figure. c, ESR1 mRNA after three",
  "independent sgRNAs against ESR1 versus scrambled control, with and without enzalutamide. d, Xenografts at endpoint from mice fed vehicle or tamoxifen diet (500 mg/kg tamoxifen diet, ad libitum).",
  "e-g, SomaScan plasma proteomics in mice without tumours (Normal) and in mice co-inoculated with adipose-MSC and 22Rv1 given vehicle or tamoxifen (n = 3 per group), with aptamers grouped by the",
  "cell source their target is assigned to. e, f, recovered absolute scale, after removing the per-sample dilution factor estimated from the floor-free high-abundance plateau; g, centred log-ratio",
  "composition, which is dilution-invariant. Two-sided t-tests, BH-adjusted.", sep = "\n"), 30, 8, 5.4)
dev.off()
cat("written", f, "\n")

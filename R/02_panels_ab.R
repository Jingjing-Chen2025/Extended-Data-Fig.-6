# EXTENDED DATA FIG. 6, PANELS a-b | TF occupancy at the MDK and TNFSF12 promoters, and where those binding
# sites sit relative to active chromatin in the adipose lineage.
#
# This replaces the previous panel a (36 hand-picked TFs ranked by a composite "binding intensity" that could
# not be computed from the file it read - remap2022_nr_macs2_hg38 is BED9 with no signalValue, so column 7,
# thickStart, was read as signal and correlates with peak start at r = 1.0000) and the former panel c (a rank
# plot, now dropped).
#
#   a  every TF with a ReMap peak in the window, on two axes kept separate: how MANY sites it has, and how
#      strongly marked the chromatin around its best site is. No composite, because combining the two is what
#      made the previous metric opaque.
#   b  Roadmap H3K27ac / H3K4me1 read density across each promoter window, with the actual binding sites of the
#      state-associated factors from panel a drawn underneath on the same coordinates. This is what makes the
#      panel informative: it shows WHERE each factor binds relative to the marks, rather than only that the
#      promoters are active.
#
# The two TNFSF12 promoter windows are plotted separately rather than as one span, because treating them as a
# single interval would place the intergenic gap between them inside the plotted window.
#
# Usage: Rscript R/02_panels_ab.R
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(ggrepel); library(patchwork); library(data.table)})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
V   <- wd("remap_accessible")
HIST <- Sys.getenv("ED6_HIST", wd("hist"))
OUT <- ED6_OUT
# palette comes from common.R
th <- theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(axis.line = element_line(linewidth = .3, colour = MUTED), axis.ticks = element_line(linewidth = .3, colour = MUTED),
        axis.text = element_text(colour = INK), axis.title = element_text(size = 6.2, colour = MUTED),
        plot.title = element_text(face = "bold", size = 8), plot.subtitle = element_text(size = 5.4, colour = MUTED),
        strip.background = element_rect(fill = "grey95", colour = NA), strip.text = element_text(size = 6, colour = INK),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 5.6), plot.title.position = "plot")
SETCOL <- c(`PI16+` = BLUE, `POSTN+` = RED, other = GREY)

R <- bind_rows(read.csv(file.path(V, "MDK_accessible_tf_ranking.csv")) %>% mutate(gene = "MDK"),
               read.csv(file.path(V, "TNFSF12_accessible_tf_ranking.csv")) %>% mutate(gene = "TNFSF12")) %>%
  mutate(set = ifelse(is.na(set) | set == "", "other", set),
         set = factor(set, c("PI16+", "POSTN+", "other")), gene = factor(gene, c("MDK", "TNFSF12")))

# ---------------- a
labme <- R %>% group_by(gene) %>% filter(set != "other" | rank_sites <= 5 | rank_access <= 3) %>% ungroup()
pa <- ggplot(R, aes(n_sites, reads_per_kb)) +
  geom_point(data = ~subset(.x, set == "other"), colour = GREY, size = .7, stroke = 0, alpha = .55) +
  geom_point(data = ~subset(.x, set != "other"), aes(colour = set), size = 1.7, stroke = 0) +
  geom_text_repel(data = labme, aes(label = TF, colour = set), size = 1.85, box.padding = .22,
                  min.segment.length = .25, segment.size = .2, segment.colour = "grey70", max.overlaps = 30, seed = 1) +
  facet_wrap(~ gene, scales = "free_x") +
  scale_colour_manual(values = SETCOL, name = NULL, breaks = c("PI16+", "POSTN+")) +
  labs(title = "a  Transcription-factor occupancy at the two ligand promoters, over the full ReMap TF universe",
       subtitle = paste("Every TF with at least one ReMap peak in the -2,000/+500 promoter window (584 at MDK, 396 at TNFSF12).",
                        "Axes are kept separate: how MANY sites a factor has, and how strongly the chromatin around its best site is marked.", sep = "\n"),
       x = "binding sites in the promoter window (distinct ReMap peaks)",
       y = expression("accessibility of best site (histone reads kb"^-1*")")) +
  th + theme(legend.position = "bottom")

# ---------------- b: chromatin marks, with the panel-a factors' binding sites drawn on the same coordinates
PK <- fread(file.path(V, "remap_peaks_with_marks.csv"))
WIN <- data.table(window = c("MDK  chr11", "TNFSF12  chr17 (distal)", "TNFSF12  chr17 (proximal)"),
                  chrom  = c("chr11", "chr17", "chr17"),
                  lo = c(46378756, 7546508, 7554036), hi = c(46382252, 7550310, 7556535))
assign_win <- function(ch, st, en) {
  vapply(seq_along(ch), function(i) {
    j <- which(WIN$chrom == ch[i] & en[i] > WIN$lo & st[i] < WIN$hi)
    if (length(j)) WIN$window[j[1]] else NA_character_ }, character(1))
}
reads <- rbindlist(lapply(list.files(HIST, "\\.bed$", full.names = TRUE), function(f) {
  d <- fread(f, header = FALSE, col.names = c("chrom", "start", "end")); d[, track := sub("\\.bed$", "", basename(f))][] }))
reads[, c("epigenome", "mark") := tstrsplit(track, "_", fixed = TRUE)]
reads[, window := assign_win(chrom, start, end)]
reads <- reads[!is.na(window)][WIN, on = "window"][, pos_kb := ((start + end) / 2 - lo) / 1000]

# R is a data.frame built with dplyr, so it needs converting before a data.table join
RT <- as.data.table(R)[, .(TF = TF, set = as.character(set))]
RT <- unique(RT[set != "other"])
STATE <- PK[, .(gene, chrom, start, end, TF)][RT, on = "TF", nomatch = 0]
STATE[, window := assign_win(chrom, start, end)]
STATE <- STATE[!is.na(window)][WIN, on = "window"][, pos_kb := ((start + end) / 2 - lo) / 1000]
STATE <- unique(STATE[, .(window, TF, set, pos_kb)])
ord <- STATE[, .(m = mean(pos_kb)), by = TF][order(m), TF]
STATE[, TF := factor(TF, ord)]
cat(sprintf("state-TF binding sites drawn on panel b: %d across %d factors", nrow(STATE), uniqueN(STATE$TF)))

pb1 <- ggplot(reads, aes(pos_kb, colour = mark, fill = mark)) +
  geom_density(adjust = .5, linewidth = .4, alpha = .18) +
  facet_grid(epigenome ~ window, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = c(H3K27ac = RED, H3K4me1 = AMBER), name = NULL, aesthetics = c("colour", "fill")) +
  labs(title = "b  Both promoters carry active marks across essentially their whole length in adipose lineage",
       subtitle = paste("Roadmap H3K27ac and H3K4me1 read density; curves are densities, so they show where reads concentrate, not coverage.",
                        "Coverage is near-complete - 33 of 35 100-bp bins at MDK and 36 of 39 and 25 of 25 at the two TNFSF12 windows contain a read -",
                        "so 99.7% of ReMap peaks sit in marked chromatin and restricting to marked chromatin retains every TF.",
                        "E025 (ADMSC) ships no H3K27ac in this Roadmap release and contributes H3K4me1 only.", sep = "\n"),
       x = NULL, y = "read density") +
  th + theme(legend.position = "right", axis.text.y = element_blank(), axis.ticks.y = element_blank(),
             axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(4, 6, 0, 3))

pb2 <- ggplot(STATE, aes(pos_kb, TF, colour = set)) +
  geom_point(shape = 124, size = 2.4, stroke = .8) +
  facet_grid(. ~ window, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = SETCOL, name = NULL, breaks = c("PI16+", "POSTN+")) +
  labs(x = "position in promoter window (kb from window start; window is -2,000 to +500 of the TSS)", y = NULL) +
  th + theme(legend.position = "right", strip.text = element_blank(), strip.background = element_blank(),
             axis.text.y = element_text(size = 5.4), panel.grid.major.y = element_line(linewidth = .2, colour = "grey92"),
             plot.margin = margin(0, 6, 3, 3))

fig <- pa / pb1 / pb2 + plot_layout(heights = c(1, .95, .8)) +
  plot_annotation(caption = paste(
    "a, Each point is one transcription factor with at least one ReMap 2022 non-redundant peak in the promoter window. Horizontal axis, number of distinct peaks; vertical axis, Roadmap histone read density over",
    "that factor's most strongly marked peak. PI16+ and POSTN+ sets are the factors nominated for each stromal state by the per-dataset CollecTRI hub analysis and by differential expression across all human stroma.",
    "b, Upper, density of H3K27ac and H3K4me1 reads across each promoter window in three adipose-lineage Roadmap epigenomes. Lower, the individual ReMap binding sites of the state-associated factors from a, on the",
    "same coordinates; each tick is one peak. The two TNFSF12 promoter windows are shown separately because they are not contiguous. The ReMap non-redundant file carries no signal values, so no binding-intensity",
    "estimate is reported anywhere in this figure; peak number, peak position and chromatin marking are the quantities the data support. Full per-TF tables in remap_accessible/.", sep = "\n"),
    theme = theme(plot.caption = element_text(size = 4.6, colour = MUTED, hjust = 0, lineheight = 1.35)))
f <- file.path(OUT, "ED6 panels ab - promoter occupancy.pdf")
ggsave(f, fig, width = 7.4, height = 7.6, device = cairo_pdf)
cat("written", f, "\n")

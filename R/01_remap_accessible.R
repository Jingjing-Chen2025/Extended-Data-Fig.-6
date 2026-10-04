# TF BINDING SITES AT THE MDK AND TNFSF12 PROMOTERS, RESTRICTED TO ADIPOSE-LINEAGE ACTIVE CHROMATIN
#
# The published Extended Data Fig. 5a ranks 36 hand-picked TFs by a composite "binding intensity" that cannot be
# computed from the file it reads: remap2022_nr_macs2_hg38 is BED9 with no signalValue, so the pipeline's
# V7 is thickStart, a genomic coordinate (r = 1.0000 with peak start). Its cell-type arm ran but changes nothing,
# because the histone marks are used only to decide whether a whole promoter is "accessible" and the globally
# computed TF scores are then copied into every tissue unchanged - all 60 rows of the by-celltype table carry
# identical values across adipose_sub, admsc and skeletal muscle.
#
# This script uses the epigenetic data the way the panel is captioned as using it: a ReMap peak counts only if it
# falls in chromatin carrying an active mark in an ADIPOSE-LINEAGE epigenome, and TFs are then ranked by how many
# such sites they have and how strongly marked those sites are.
#
#   epigenomes   E063 subcutaneous adipose, E065 visceral adipose, E025 adipose-derived MSC.
#   marks        H3K27ac and H3K4me1. NOTE: E025-H3K27ac does not exist in the Roadmap release held here (that
#                epigenome ships H3K9ac instead), so ADMSC contributes H3K4me1 only. Stated rather than silently
#                substituted.
#   ACCESSIBILITY  reads per kb of peak, from raw tagAlign. These are genome-wide read files and a 3.5 kb
#                promoter receives few of them - 475 reads across all five tracks here - so the measure is sparse
#                and is reported as a rank, never as a quantitative signal.
#   FREQUENCY    number of distinct ReMap peaks for that TF inside marked chromatin, and the number of distinct
#                biotypes those peaks come from.
#   RANKING      TFs are ranked on both axes and the two are reported side by side rather than combined into one
#                composite, because combining them is what made the published metric uninterpretable.
#
# Inputs  ReMap 2022 non-redundant peaks ($ED6_REMAP, bgzipped + tabix), Roadmap tagAlign-derived BED files in
#         $ED6_HIST, and the PI16+/POSTN+ state transcription factors in data/state_TFs_fig1_2.csv
# Output  $ED6_DATA/remap_accessible/{remap_peaks_with_marks.csv, MDK_accessible_tf_ranking.csv, TNFSF12_...}
# Usage   Rscript R/01_remap_accessible.R
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(GenomicRanges); library(Rsamtools)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db); library(AnnotationDbi)
})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
OUT <- wd("remap_accessible"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
BGZ <- ED6_REMAP
ED6_HIST <- Sys.getenv("ED6_HIST", wd("hist"))
set.seed(42); options(width = 200)

ez <- AnnotationDbi::select(org.Hs.eg.db, keys = c("MDK", "TNFSF12"), keytype = "SYMBOL", columns = "ENTREZID")
txb <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")
PROM <- do.call(c, lapply(seq_len(nrow(ez)), function(i) {
  p <- GenomicRanges::reduce(promoters(txb[[ez$ENTREZID[i]]], upstream = 2000, downstream = 500))
  mcols(p)$gene <- ez$SYMBOL[i]; p }))

raw <- Rsamtools::scanTabix(Rsamtools::TabixFile(BGZ), param = PROM)
pk <- rbindlist(lapply(seq_along(raw), function(i) {
  if (!length(raw[[i]])) return(NULL)
  d <- fread(text = paste(raw[[i]], collapse = "\n"), sep = "\t", header = FALSE)
  d[, gene := PROM$gene[i]][] }), fill = TRUE)
setnames(pk, c("V1","V2","V3","V4","V5","V6","V7","V8","V9"),
         c("chrom","start","end","name","score","strand","thickStart","thickEnd","itemRgb"))
pk[, TF := toupper(sub(":.*$", "", name))][, biotype := sub("^[^:]*:", "", name)][, width := end - start]
pk[, pid := .I]
GP <- GRanges(pk$chrom, IRanges(pk$start + 1, pk$end))
msg("ReMap peaks in the two promoters: %d over %d TFs", nrow(pk), uniqueN(pk$TF))

TRACKS <- c("E063_H3K27ac", "E063_H3K4me1", "E065_H3K27ac", "E065_H3K4me1", "E025_H3K4me1")
for (t in TRACKS) {
  f <- file.path(ED6_HIST, paste0(t, ".bed"))
  if (!file.exists(f) || file.size(f) == 0) { pk[[t]] <- 0L; warn_msg("track %s absent", t); next }
  r <- fread(f, header = FALSE, col.names = c("chrom", "start", "end"))
  gr <- GRanges(r$chrom, IRanges(r$start + 1, r$end))
  pk[[t]] <- GenomicRanges::countOverlaps(GP, gr)
}
pk[, reads_total := rowSums(.SD), .SDcols = TRACKS]
pk[, reads_per_kb := 1000 * reads_total / width]
pk[, marked := reads_total >= 1]
pk[, marked3 := reads_total >= 3]
msg("histone reads loaded: %s", paste(sprintf("%s=%d", TRACKS, colSums(pk[, ..TRACKS])), collapse = " "))
msg("peaks in marked chromatin: %d of %d (>=1 read), %d (>=3 reads)", sum(pk$marked), nrow(pk), sum(pk$marked3))
write.csv(pk[, .(gene, chrom, start, end, TF, biotype, width, score, reads_total, reads_per_kb, marked, marked3)],
          file.path(OUT, "remap_peaks_with_marks.csv"), row.names = FALSE)

# The PI16+ / POSTN+ state transcription factors, from the Figure 1-2 single-cell analysis. The table shipped in
# data/ is the finished product of that chain (see README); it is read here rather than recomputed.
ST <- read.csv(repo("data", "state_TFs_fig1_2.csv"), stringsAsFactors = FALSE)
PI16  <- unique(toupper(ST$gene[ST$state == "PI16_progenitor"]))
POSTN <- unique(toupper(ST$gene[ST$state == "POSTN_activated"]))
lab <- function(x) ifelse(x %in% PI16 & x %in% POSTN, "both", ifelse(x %in% PI16, "PI16+", ifelse(x %in% POSTN, "POSTN+", "")))

for (gn in c("MDK", "TNFSF12")) {
  d <- pk[gene == gn]
  acc <- d[marked == TRUE, .(n_sites = .N, n_biotypes = uniqueN(biotype),
                             reads_per_kb = round(max(reads_per_kb), 1),
                             best_track = TRACKS[which.max(colSums(.SD))]), by = TF, .SDcols = TRACKS]
  all_n <- d[, .(n_all = .N), by = TF]
  acc <- merge(acc, all_n, by = "TF")[, set := lab(TF)][order(-n_sites, -reads_per_kb)]
  acc[, rank_sites := .I][, rank_access := frank(-reads_per_kb, ties.method = "min")]
  write.csv(acc, file.path(OUT, sprintf("%s_accessible_tf_ranking.csv", gn)), row.names = FALSE)
  cat(sprintf("\n================ %s promoter: TFs with a binding site in adipose-lineage marked chromatin\n", gn))
  cat(sprintf("   %d of %d TFs retain >=1 site after the chromatin filter\n", nrow(acc), uniqueN(d$TF)))
  cat("\n-- top 20 by NUMBER of marked binding sites (frequency)\n")
  print(as.data.frame(acc[1:min(20, .N), .(TF, set, n_sites, n_all, n_biotypes, reads_per_kb, rank_access)]), row.names = FALSE)
  cat("\n-- top 12 by ACCESSIBILITY of the best site (reads per kb)\n")
  print(as.data.frame(acc[order(-reads_per_kb)][1:min(12, .N), .(TF, set, reads_per_kb, n_sites, n_biotypes, rank_sites)]), row.names = FALSE)
  cat("\n-- the state TFs, wherever they fall\n")
  s <- acc[set != ""][order(-n_sites)]
  if (nrow(s)) print(as.data.frame(s[, .(TF, set, n_sites, n_all, n_biotypes, reads_per_kb, rank_sites, rank_access)]), row.names = FALSE)
  else cat("   none of the state TFs retains a site in marked chromatin\n")
}
msg("written to %s", OUT)

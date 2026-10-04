# EXTENDED DATA FIG. 6, PANELS e-g | SomaScan plasma proteomics.
#
# INCOMPLETE - READ THIS FIRST. This is the SomaScan analysis exactly as it is held in
# `Nature Manuscript/Extended Data Fig. 5/SomaScan.docx`, extracted verbatim. That Word document stores only part
# of the script: it stops mid-line inside the recovery-factor diagnostic plot, so the per-Source panel loop that
# actually writes `Source_comparison_RECOVERED_Normal_Vehicle_Tamoxifen.pdf` is missing. Nothing has been added
# or guessed here. The three panels in the figure are pages 4, 8 and 15 of that PDF:
#     e  Source: ASPC       recovered absolute scale
#     f  Source: NEPC       recovered absolute scale
#     g  Source: PCa        compositional (CLR)
# They were identified by matching axis ranges and P values against the per-Source pages, because the panels in
# the artwork carry no Source label of their own. Confirm against the original before relying on them.
#
# The missing tail needs to be recovered from the author's working copy of the script.
#
# Two lines were changed, and only two: the absolute paths to the ADAT file and to the aptamer-to-Source mapping
# spreadsheet are read from environment variables instead of being hard-coded to a local home directory.
# Everything else is as stored in the document, including the point at which it stops.
#
# NOT RUNNABLE as it stands - it is kept here as the record of the analysis, not as a build step.

# ============================================================
# SomaScan ASPC 11K — RECOVER ~8-12x Vehicle dilution + report co-existing
#   xenograft-driven COMPOSITIONAL remodeling.
#   Per-SOURCE (cell-origin) comparison: Normal vs Vehicle vs Tamoxifen
#   Output: PDF + Excel
# ------------------------------------------------------------
# Groups (plasma):
#   Normal     : mouse WITHOUT tumors                          [Normal_1..3]
#   Vehicle    : adipose-MSC + 22Rv1 co-inoculated, vehicle    [Vehicle_1..3]
#   Tamoxifen  : adipose-MSC + 22Rv1 co-inoculated, tamoxifen  [Tam_1..3]
#
# Inputs:
#   1) ADAT : ...hybNorm.medNormInt.plateScale...adat (Sample + Buffer wells)
#   2) XLSX : SomaScan_ASPC_11K.xlsx (Sequence ID, Symbol, Source)
#
# MODEL (two effects co-exist; both correct):
#   (A) TECHNICAL DILUTION of Vehicle plasma -> per-sample multiplicative
#       factor (~8-12x) lowering ALL proteins equally. ESTABLISHED FIRST via
#       the ratio-vs-abundance (floor) diagnostic: the naive median ratio is
#       floor-biased toward 1; the TRUE factor is read off the HIGH-ABUNDANCE
#       PLATEAU (floor-free aptamers). RECOVER by multiplying Vehicle back up.
#   (B) XENOGRAFT-DRIVEN BIOLOGY -> genuine COMPOSITIONAL remodeling
#       (within-sample protein RATIOS change; tamoxifen partially reverses).
#       INDEPENDENT of (A): ratios are dilution-invariant.
#
#   => WORKFLOW: (Section 7) compute & justify the SCALING FACTOR from the
#      floor diagnostic; (Section 8) APPLY recovery using that factor.
#      Report BOTH recovered absolute scale (A) and CLR composition (B).
#
#   Caveat: human aptamers on MOUSE plasma -> conserved mouse orthologs
#   (host-proteome). Plateau factor is used as the dilution estimate.
#
# OUTPUTS:
#   PDF : floor-diagnostic (ratio vs abundance) + recovery-factor page;
#         recovered-scale overview + per-Source panels;
#         compositional (CLR) overview + per-Source panels.
#   XLSX: scaling-factor evidence, recovery factors, recovered per-sample +
#         group medians + fold-changes, CLR per-sample + effects + stats,
#         per-Source plateau factors, aptamer table, README.
# ============================================================
# Date: 2026-06-22
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(matrixStats)
  library(rstatix)
  library(ggpubr)
  library(SomaDataIO)
  library(openxlsx)
})
# -----------------------------
# OPTIONS
# -----------------------------
adat_path <- Sys.getenv("SOMASCAN_ADAT")   # ...hybNorm.medNormInt.plateScale....adat
map_path  <- Sys.getenv("SOMASCAN_MAP")    # SomaScan_ASPC_11K.xlsx (Sequence ID, Symbol, Source)
out_dir   <- dirname(map_path)
pdf_path  <- file.path(out_dir, "Source_comparison_RECOVERED_Normal_Vehicle_Tamoxifen.pdf")
xlsx_path <- file.path(out_dir, "Source_comparison_RECOVERED_Normal_Vehicle_Tamoxifen.xlsx")
GROUP_ORDER     <- c("Normal", "Vehicle", "Tamoxifen")
REF_GROUP       <- "Normal"
FOCUS_SOURCE    <- "ASPC"     # source used for the floor/scaling diagnostic
FLOOR_MULT      <- 2          # keep aptamers > FLOOR_MULT x buffer floor in ALL groups
PLATEAU_FRAC    <- 0.25       # top fraction (by Normal abundance) defines floor-free plateau
RECOVERY_BASIS  <- "global"   # "global" = one dilution factor per sample (recommended)
                              # "per_source" = recover each Source by its own plateau factor
MIN_APT_SRC     <- 2          # a Source needs >= this many kept aptamers to be analyzed
PAIRWISE_TEST   <- "t"        # "t" | "wilcox"
P_ADJUST        <- "BH"
eps             <- 1e-8
# -----------------------------
# 1) LOAD ADAT + 2) MAP + 3) HARMONIZE KEYS
# -----------------------------
if (!file.exists(adat_path)) stop("ADAT not found:\n", adat_path)
adat <- SomaDataIO::read_adat(adat_path)
meta_fields  <- SomaDataIO::getMeta(adat)
analyte_anno <- SomaDataIO::getAnalyteInfo(adat)
adat_df <- as.data.frame(adat, stringsAsFactors = FALSE); rownames(adat_df) <- NULL
message("ADAT: ", nrow(adat), " rows x ", ncol(adat), " cols")
if (!file.exists(map_path)) stop("Excel not found:\n", map_path)
map_raw <- readxl::read_excel(map_path)
nm <- gsub("\\s+", "", tolower(colnames(map_raw)))
seqid_col  <- colnames(map_raw)[nm %in% c("sequenceid","seqid","seq.id","sequence_id")]
symbol_col <- colnames(map_raw)[nm %in% c("symbol","genesymbol","gene")]
source_col <- colnames(map_raw)[nm %in% c("source","celltype","cell_type")]
if (length(seqid_col)  != 1) stop("Need a unique 'Sequence ID' column.")
if (length(source_col) != 1) stop("Need a unique 'Source' column.")
if (length(symbol_col) != 1) symbol_col <- NA_character_
map_df <- map_raw %>%
  transmute(seq_id = as.character(.data[[seqid_col]]),
            symbol = if (!is.na(symbol_col)) as.character(.data[[symbol_col]]) else NA_character_,
            source = as.character(.data[[source_col]])) %>%
  filter(!is.na(seq_id), !is.na(source), source != "")
message("Mapping rows: ", nrow(map_df),
        " | Sources: ", paste(sort(unique(map_df$source)), collapse = ", "))
analyte_cols <- grep("^seq\\.", colnames(adat), value = TRUE)
if (length(analyte_cols) == 0) analyte_cols <- setdiff(colnames(adat), meta_fields)
anno_keys <- analyte_anno %>% as.data.frame() %>% {
  df <- .
  apt_col <- intersect(c("AptName","Apt","aptName"), colnames(df))
  seq_col <- intersect(c("SeqId","SeqID","Seq.Id","seqid"), colnames(df))
  if (length(apt_col) == 0) apt_col <- colnames(df)[grepl("apt",   tolower(colnames(df)))][1]
  if (length(seq_col) == 0) seq_col <- colnames(df)[grepl("seqid", tolower(colnames(df)))][1]
  tibble(AptName = as.character(df[[apt_col]]), SeqId = as.character(df[[seq_col]]))
}
canon <- function(x){ x <- as.character(x); x <- sub("^seq[._]", "", x, ignore.case = TRUE); gsub("[._]", "-", x) }
anno_keys <- anno_keys %>% mutate(SeqId_canon = canon(SeqId))
map_df    <- map_df    %>% mutate(seq_id_canon = canon(seq_id))
col_source <- tibble(AptName = analyte_cols) %>%
  left_join(anno_keys %>% select(AptName, SeqId_canon), by = "AptName") %>%
  mutate(SeqId_canon = ifelse(is.na(SeqId_canon), canon(AptName), SeqId_canon)) %>%
  left_join(map_df %>% select(seq_id_canon, source, symbol),
            by = c("SeqId_canon" = "seq_id_canon"))
mapped_cols    <- col_source$AptName[!is.na(col_source$source)]
mapped_sources <- setNames(col_source$source[!is.na(col_source$source)], mapped_cols)
mapped_symbol  <- setNames(col_source$symbol[!is.na(col_source$source)],  mapped_cols)
message("Analytes mapped to a Source: ", length(mapped_cols),
        " across ", length(unique(mapped_sources)), " Sources")
# -----------------------------
# 4) GROUPS + BUFFER wells (for LOD)
# -----------------------------
sample_id_field <- "SampleId"
meta <- adat_df %>%
  filter(.data[["SampleType"]] == "Sample") %>%
  transmute(sample_id = as.character(.data[[sample_id_field]])) %>%
  mutate(group = case_when(
    grepl("^tam",     sample_id, ignore.case = TRUE) ~ "Tamoxifen",
    grepl("^vehicle", sample_id, ignore.case = TRUE) ~ "Vehicle",
    grepl("^normal",  sample_id, ignore.case = TRUE) ~ "Normal",
    TRUE ~ NA_character_)) %>%
  mutate(group = factor(group, levels = GROUP_ORDER))
if (any(is.na(meta$group)))
  stop("Unmatched SampleIds: ", paste(meta$sample_id[is.na(meta$group)], collapse = ", "))
message("Group sizes:\n", paste(capture.output(print(table(meta$group))), collapse = "\n"))
rfu_lin <- adat_df %>%
  mutate(sample_id = as.character(.data[[sample_id_field]])) %>%
  filter(sample_id %in% meta$sample_id) %>%
  select(sample_id, all_of(mapped_cols)) %>%
  column_to_rownames("sample_id") %>% as.matrix()
rfu_lin[rfu_lin <= 0] <- NA
samp_ids <- rownames(rfu_lin)
grp_of   <- setNames(as.character(meta$group)[match(samp_ids, meta$sample_id)], samp_ids)
# ============================================================
# 5) LOD / BUFFER FLOOR  (per-aptamer)
# ============================================================
have_buffer <- any(adat_df[["SampleType"]] == "Buffer")
if (have_buffer) {
  buf_mat <- adat_df %>% filter(.data[["SampleType"]] == "Buffer") %>%
    select(all_of(mapped_cols)) %>% as.matrix()
  buf_mat[buf_mat <= 0] <- NA
  floor_per_apt <- matrixStats::colMedians(buf_mat, na.rm = TRUE)
  floor_global  <- median(floor_per_apt, na.rm = TRUE)
  message("\nBuffer (LOD) floor: global median = ", round(floor_global, 1),
          " RFU (n_buffer = ", sum(adat_df[["SampleType"]] == "Buffer"), ")")
} else {
  floor_per_apt <- setNames(rep(quantile(rfu_lin[grp_of == REF_GROUP, ], 0.02, na.rm = TRUE),
                                ncol(rfu_lin)), colnames(rfu_lin))
  floor_global <- median(floor_per_apt, na.rm = TRUE)
  warning("No Buffer wells -> using 2nd percentile of Normal as floor proxy (", round(floor_global,1), ").")
}
# ============================================================
# 6) FLOOR FILTER (keep aptamers > FLOOR_MULT x floor in ALL groups)
# ============================================================
grp_med <- sapply(GROUP_ORDER, function(g)
  matrixStats::colMedians(rfu_lin[grp_of == g, , drop = FALSE], na.rm = TRUE))
colnames(grp_med) <- GROUP_ORDER
thresh <- FLOOR_MULT * floor_per_apt[rownames(grp_med)]
above_floor_all <- (grp_med[, "Normal"] > thresh) & (grp_med[, "Vehicle"] > thresh) &
                   (grp_med[, "Tamoxifen"] > thresh)
keep_cols <- rownames(grp_med)[above_floor_all]
message("\nFloor filter (", FLOOR_MULT, "x floor in ALL groups): kept ",
        length(keep_cols), " / ", length(mapped_cols), " mapped aptamers.")
ratio_tbl <- tibble(
  AptName   = rownames(grp_med),
  Source    = mapped_sources[rownames(grp_med)],
  Symbol    = mapped_symbol[rownames(grp_med)],
  Normal    = grp_med[, "Normal"],
  Vehicle   = grp_med[, "Vehicle"],
  Tamoxifen = grp_med[, "Tamoxifen"],
  Veh_over_Normal = grp_med[, "Vehicle"]   / grp_med[, "Normal"],
  Tam_over_Normal = grp_med[, "Tamoxifen"] / grp_med[, "Normal"],
  floor     = floor_per_apt[rownames(grp_med)],
  kept      = rownames(grp_med) %in% keep_cols) %>%
  arrange(Source, desc(Normal))
# ============================================================
# 7) *** SCALING FACTOR — established FIRST via floor diagnostic ***
#    "Tie each ASPC protein's Veh/Normal ratio to its ABUNDANCE (Normal)."
#    Floor artifact -> low-abundance ratios pulled toward 1; high-abundance
#    plateau reveals the TRUE scaling factor (~8-12x). We:
#      (a) test the abundance~ratio dependence (Spearman),
#      (b) flag near-floor aptamers (per-aptamer LOD),
#      (c) READ the scaling factor off the floor-free high-abundance plateau,
#      (d) derive per-sample dilution factors used by Section 8 for recovery.
# ============================================================
focus_floor <- ratio_tbl %>% filter(Source == FOCUS_SOURCE) %>%
  mutate(log10_abundance = log10(Normal),
         near_floor_Veh  = Vehicle < FLOOR_MULT * floor,
         floor_free      = (Vehicle > FLOOR_MULT * floor) & (Normal > FLOOR_MULT * floor),
         status          = ifelse(floor_free, "floor-free", "near floor")) %>%
  arrange(Normal)
# (a) abundance ~ ratio dependence (floor signature if strongly negative)
spear <- with(focus_floor,
              suppressWarnings(cor(log10_abundance, Veh_over_Normal,
                                   use = "complete.obs", method = "spearman")))
# (b) near-floor count (per-aptamer LOD)
n_near  <- sum(focus_floor$near_floor_Veh, na.rm = TRUE)
n_focus <- nrow(focus_floor)
# (c) TRUE scaling factor from the floor-free high-abundance plateau
plateau_focus <- focus_floor %>% filter(floor_free) %>%
  arrange(desc(Normal)) %>%
  slice_head(n = max(5, ceiling(PLATEAU_FRAC * sum(focus_floor$floor_free))))
plateau_ratio_med <- median(plateau_focus$Veh_over_Normal, na.rm = TRUE)
plateau_ratio_lo  <- as.numeric(quantile(plateau_focus$Veh_over_Normal, 0.25, na.rm = TRUE))
plateau_ratio_hi  <- as.numeric(quantile(plateau_focus$Veh_over_Normal, 0.75, na.rm = TRUE))
naive_ratio_med   <- median(focus_floor$Veh_over_Normal, na.rm = TRUE)
message("\n=== *** SCALING FACTOR (floor diagnostic; ", FOCUS_SOURCE, ") *** ===")
message("Spearman(log10 abundance, Veh/Normal) = ", round(spear, 3),
        "  (strong negative => detection-floor artifact)")
message("Near-floor in Vehicle (< ", FLOOR_MULT, "x per-aptamer LOD): ",
        n_near, " / ", n_focus, "  (", round(100*n_near/n_focus), "% background noise)")
message("Naive median ratio (ALL ", n_focus, ", floor-biased): ", round(naive_ratio_med, 3),
        "  (~", round(1/naive_ratio_med, 1), "x lower)")
message("PLATEAU ratio (top ", nrow(plateau_focus), " floor-free): ", round(plateau_ratio_med, 3),
        "  [IQR ", round(plateau_ratio_lo, 3), "-", round(plateau_ratio_hi, 3), "]")
message("==> TRUE shared scaling factor ~ ",
        round(1/plateau_ratio_hi, 1), "-", round(1/plateau_ratio_lo, 1),
        "x lower in Vehicle  (point est. ~", round(1/plateau_ratio_med, 1), "x)")
# Evidence table exported to Excel
scaling_evidence <- tibble(
  metric = c("Spearman(abundance,ratio)", "near_floor_Vehicle_n", "n_focus_total",
             "naive_median_ratio", "naive_fold_lower",
             "plateau_ratio_median", "plateau_ratio_Q1", "plateau_ratio_Q3",
             "TRUE_fold_lower_point", "TRUE_fold_lower_low", "TRUE_fold_lower_high",
             "buffer_floor_global_RFU"),
  value  = c(round(spear,3), n_near, n_focus,
             round(naive_ratio_med,3), round(1/naive_ratio_med,1),
             round(plateau_ratio_med,3), round(plateau_ratio_lo,3), round(plateau_ratio_hi,3),
             round(1/plateau_ratio_med,1), round(1/plateau_ratio_hi,1), round(1/plateau_ratio_lo,1),
             round(floor_global,1)))
# (d) PER-SAMPLE dilution factors (used for recovery in Section 8)
# Built from the floor-free high-abundance plateau across the KEPT panel,
# referenced to Normal (=1). This generalizes the ASPC plateau logic above
# to all kept aptamers so recovery is robust.
plateau_cols_global <- ratio_tbl %>% filter(kept) %>%
  arrange(desc(Normal)) %>%
  slice_head(n = max(5, ceiling(PLATEAU_FRAC * sum(ratio_tbl$kept)))) %>% pull(AptName)
norm_ref_plateau <- matrixStats::colMedians(
  rfu_lin[grp_of == REF_GROUP, plateau_cols_global, drop = FALSE], na.rm = TRUE)
dilution_factor <- apply(rfu_lin[, plateau_cols_global, drop = FALSE], 1, function(x)
  median(norm_ref_plateau / x, na.rm = TRUE))                 # >1 means "scale up"
dilution_factor[!is.finite(dilution_factor) | dilution_factor <= 0] <- 1
rec_tbl <- tibble(sample_id = samp_ids,
                  group = factor(grp_of[samp_ids], levels = GROUP_ORDER),
                  dilution_factor = as.numeric(dilution_factor[samp_ids]))
rec_group <- rec_tbl %>% group_by(group) %>%
  summarise(mean_factor = mean(dilution_factor, na.rm = TRUE),
            sd_factor   = sd(dilution_factor,   na.rm = TRUE),
            min_factor  = min(dilution_factor,  na.rm = TRUE),
            max_factor  = max(dilution_factor,  na.rm = TRUE), n = n(), .groups = "drop")
veh_fac <- rec_group %>% filter(group == "Vehicle") %>% pull(mean_factor)
message("\n=== PER-SAMPLE dilution factors (plateau-based; Normal=1) ===")
print(rec_tbl, n = Inf); print(rec_group, n = Inf)
message("==> Vehicle scaling factor ~", round(veh_fac, 1),
        "x (consistent with the ASPC plateau ~8-12x). Used for recovery in Section 8.")
# Per-Source plateau factor (reporting + per_source recovery)
plateau_by_src <- ratio_tbl %>% filter(kept) %>% group_by(Source) %>%
  arrange(desc(Normal), .by_group = TRUE) %>%
  mutate(rank_ab = row_number(), n_src = n()) %>%
  filter(rank_ab <= pmax(3, ceiling(PLATEAU_FRAC * n_src))) %>%
  summarise(plateau_Veh_over_Normal = median(Veh_over_Normal, na.rm = TRUE),
            plateau_Tam_over_Normal = median(Tam_over_Normal, na.rm = TRUE),
            n_plateau = n(), .groups = "drop") %>%
  mutate(Vehicle_recovery_factor   = 1 / plateau_Veh_over_Normal,
         Tamoxifen_recovery_factor = 1 / plateau_Tam_over_Normal)
# Floor-diagnostic plot (saved later into the PDF)
p_floor <- ggplot(focus_floor, aes(log10_abundance, Veh_over_Normal)) +
  annotate("rect", xmin = log10(min(plateau_focus$Normal)), xmax = Inf,
           ymin = -Inf, ymax = Inf, alpha = 0.07, fill = "#3B7DD8") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = plateau_ratio_med, linetype = "dotted",
             color = "#3B7DD8", linewidth = 0.8) +
  geom_point(aes(color = status), size = 2.6, alpha = 0.85) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.7) +
  scale_color_manual(values = c(`floor-free` = "#222222", `near floor` = "#C0392B"), name = NULL) +
  labs(title = paste0(FOCUS_SOURCE, " aptamers: Veh/Normal ratio vs abundance (scaling factor)"),
       subtitle = paste0("Spearman = ", round(spear, 2),
                        "  |  plateau (blue) ~", round(1/plateau_ratio_hi,1), "-",
                        round(1/plateau_ratio_lo,1), "x lower  |  red = near-floor (artifactual)"),
       x = "log10(Normal RFU)  [abundance]", y = "Veh / Normal ratio") +
  theme_bw()
# ============================================================
# 8) *** APPLY RECOVERY ***  using the scaling factor from Section 7
# ============================================================
if (RECOVERY_BASIS == "global") {
  rfu_rec <- sweep(rfu_lin, 1, dilution_factor[rownames(rfu_lin)], "*")
} else {
  rfu_rec <- rfu_lin
  for (s in unique(mapped_sources)) {
    cs <- intersect(names(mapped_sources)[mapped_sources == s], colnames(rfu_rec))
    if (!length(cs)) next
    vf <- plateau_by_src %>% filter(Source == s) %>% pull(Vehicle_recovery_factor)
    tf <- plateau_by_src %>% filter(Source == s) %>% pull(Tamoxifen_recovery_factor)
    if (!length(vf)) vf <- 1; if (!length(tf)) tf <- 1
    fac <- setNames(rep(1, length(samp_ids)), samp_ids)
    fac[grp_of == "Vehicle"]   <- ifelse(is.finite(vf) && vf > 0, vf, 1)
    fac[grp_of == "Tamoxifen"] <- ifelse(is.finite(tf) && tf > 0, tf, 1)
    rfu_rec[, cs] <- sweep(rfu_rec[, cs, drop = FALSE], 1, fac[rownames(rfu_rec)], "*")
  }
}
# Recovered group medians per kept aptamer (sanity: Veh/Normal -> ~1; remaining
# deviation from 1 = the co-existing BIOLOGY, not dilution)
rec_grp_med <- sapply(GROUP_ORDER, function(g)
  matrixStats::colMedians(rfu_rec[grp_of == g, keep_cols, drop = FALSE], na.rm = TRUE))
colnames(rec_grp_med) <- GROUP_ORDER
recovered_apt <- tibble(
  AptName = keep_cols, Source = mapped_sources[keep_cols], Symbol = mapped_symbol[keep_cols],
  Normal_rec = rec_grp_med[, "Normal"], Vehicle_rec = rec_grp_med[, "Vehicle"],
  Tamoxifen_rec = rec_grp_med[, "Tamoxifen"]) %>%
  mutate(Veh_over_Normal_rec = Vehicle_rec / Normal_rec,
         Tam_over_Normal_rec = Tamoxifen_rec / Normal_rec)
message("\nPost-recovery Vehicle/Normal (kept aptamers): median = ",
        round(median(recovered_apt$Veh_over_Normal_rec, na.rm = TRUE), 2),
        " (≈1 = dilution removed; spread around 1 = xenograft biology).")
# Per-Source recovered comparison (log2 recovered RFU)
log2_rec <- log2(rfu_rec[, keep_cols, drop = FALSE] + eps)
kept_src_counts <- table(mapped_sources[keep_cols])
analyze_sources <- sort(names(kept_src_counts)[kept_src_counts >= MIN_APT_SRC])
message("\nSources analyzed (>= ", MIN_APT_SRC, " kept aptamers): ",
        paste0(analyze_sources, " (", kept_src_counts[analyze_sources], ")", collapse = ", "))
score_rec <- sapply(analyze_sources, function(s) {
  cs <- intersect(names(mapped_sources)[mapped_sources == s], colnames(log2_rec))
  rowMeans(log2_rec[, cs, drop = FALSE], na.rm = TRUE)
})
score_rec <- as.matrix(score_rec); rownames(score_rec) <- rownames(log2_rec)
rec_long <- as.data.frame(score_rec) %>% rownames_to_column("sample_id") %>%
  left_join(meta, by = "sample_id") %>%
  pivot_longer(-c(sample_id, group), names_to = "Source", values_to = "log2_recovered_RFU") %>%
  mutate(group = factor(group, levels = GROUP_ORDER))
# ============================================================
# 9) CO-EXISTING BIOLOGY: CLR (compositional; dilution-invariant) per Source
# ============================================================
log_keep <- log(rfu_lin[, keep_cols, drop = FALSE] + eps)
clr_mat  <- log_keep - rowMeans(log_keep, na.rm = TRUE)
score_clr <- sapply(analyze_sources, function(s) {
  cs <- intersect(names(mapped_sources)[mapped_sources == s], colnames(clr_mat))
  rowMeans(clr_mat[, cs, drop = FALSE], na.rm = TRUE)
})
score_clr <- as.matrix(score_clr); rownames(score_clr) <- rownames(clr_mat)
clr_long <- as.data.frame(score_clr) %>% rownames_to_column("sample_id") %>%
  left_join(meta, by = "sample_id") %>%
  pivot_longer(-c(sample_id, group), names_to = "Source", values_to = "CLR") %>%
  mutate(group = factor(group, levels = GROUP_ORDER))
# ============================================================
# 10) STATS per Source (recovered AND CLR)
# ============================================================
stat_block <- function(long_df, value_col) {
  kw <- long_df %>% group_by(Source) %>%
    rstatix::kruskal_test(reformulate("group", value_col)) %>%
    rstatix::adjust_pvalue(method = P_ADJUST) %>%
    rstatix::add_significance("p.adj") %>% ungroup()
  pw <- long_df %>% group_by(Source) %>%
    { if (PAIRWISE_TEST == "t")
        rstatix::t_test(., reformulate("group", value_col), p.adjust.method = P_ADJUST)
      else
        rstatix::wilcox_test(., reformulate("group", value_col), p.adjust.method = P_ADJUST) } %>%
    rstatix::add_significance("p.adj") %>% ungroup()
  list(kw = kw, pw = pw)
}
S_rec <- stat_block(rec_long, "log2_recovered_RFU")
rec_means <- rec_long %>% group_by(Source, group) %>%
  summarise(mean_log2 = mean(log2_recovered_RFU, na.rm = TRUE),
            sd_log2 = sd(log2_recovered_RFU, na.rm = TRUE), n = n(), .groups = "drop")
rec_fc <- rec_means %>% select(Source, group, mean_log2) %>%
  pivot_wider(names_from = group, values_from = mean_log2) %>%
  mutate(log2FC_Vehicle_vs_Normal    = Vehicle   - Normal,
         log2FC_Tamoxifen_vs_Normal  = Tamoxifen - Normal,
         log2FC_Tamoxifen_vs_Vehicle = Tamoxifen - Vehicle,
         FC_Vehicle_vs_Normal        = 2^log2FC_Vehicle_vs_Normal,
         FC_Tamoxifen_vs_Normal      = 2^log2FC_Tamoxifen_vs_Normal,
         FC_Tamoxifen_vs_Vehicle     = 2^log2FC_Tamoxifen_vs_Vehicle)
S_clr <- stat_block(clr_long, "CLR")
clr_means <- clr_long %>% group_by(Source, group) %>%
  summarise(mean_CLR = mean(CLR, na.rm = TRUE), sd_CLR = sd(CLR, na.rm = TRUE),
            n = n(), .groups = "drop")
clr_fc <- clr_means %>% select(Source, group, mean_CLR) %>%
  pivot_wider(names_from = group, values_from = mean_CLR) %>%
  mutate(Vehicle_vs_Normal    = Vehicle   - Normal,
         Tamoxifen_vs_Normal  = Tamoxifen - Normal,
         Tamoxifen_vs_Vehicle = Tamoxifen - Vehicle,
         Tam_reverses_Vehicle = sign(Tamoxifen_vs_Normal) != sign(Vehicle_vs_Normal) |
                                (abs(Tamoxifen_vs_Normal) < abs(Vehicle_vs_Normal)))
message("\n=== RECOVERED-scale per-Source fold-changes (dilution removed) ===")
print(rec_fc %>% select(Source, log2FC_Vehicle_vs_Normal, log2FC_Tamoxifen_vs_Normal,
                        FC_Vehicle_vs_Normal, FC_Tamoxifen_vs_Normal), n = Inf)
message("\n=== CO-EXISTING COMPOSITION (CLR) per-Source effect sizes ===")
print(clr_fc %>% select(Source, Vehicle_vs_Normal, Tamoxifen_vs_Normal,
                        Tamoxifen_vs_Vehicle, Tam_reverses_Vehicle), n = Inf)
# ============================================================
# 11) PLOTS -> PDF
# ============================================================
theme_set(theme_bw(base_size = 12))
grp_cols <- c(Normal = "#E64B35", Vehicle = "#3C9E40", Tamoxifen = "#3B7DD8")
panel_one <- function(long_df, value_col, pw_tab, kw_tab, src, ylab, title_tag, subtag = "") {
  d <- long_df %>% filter(Source == src) %>% rename(val = !!value_col)
  if (nrow(d) == 0 || all(is.na(d$val))) return(NULL)
  ymax <- max(d$val, na.rm = TRUE); yrng <- diff(range(d$val, na.rm = TRUE))
  yrng <- ifelse(is.finite(yrng) && yrng > 0, yrng, abs(ymax) + 1)
  pw <- pw_tab %>% filter(Source == src) %>%
    mutate(group1 = factor(group1, levels = GROUP_ORDER),
           group2 = factor(group2, levels = GROUP_ORDER),
           step = row_number(), y.position = ymax + step * 0.13 * yrng)
  kwp <- kw_tab %>% filter(Source == src) %>% pull(p.adj)
  ggplot(d, aes(group, val, color = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.25) +
    geom_jitter(width = 0.10, size = 2.2, alpha = 0.9) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "gray25") +
    scale_color_manual(values = grp_cols) +
    ggpubr::stat_pvalue_manual(pw, label = "p = {round(p,3)}",
                               tip.length = 0.01, bracket.size = 0.3, size = 3) +
    labs(title = paste0(title_tag, " — Source: ", src),
         subtitle = paste0("KW p.adj = ", ifelse(length(kwp)>0, signif(kwp,3), NA), subtag),
         y = ylab, x = NULL) +
    theme(legend.position = "none")
}
# Recovery-factor diagnostic
p_factor <- rec_tbl %>%
  ggplot(aes(group, dilution_factor, color = group)) +
  geom_jitter(width = 0.08, size = 2.8, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = grp_cols) +
  labs(title = "Per-sample DILUTION/SCALING factor (pl*


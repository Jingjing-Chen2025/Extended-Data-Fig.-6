# Shared configuration for the Extended Data Fig. 6 scripts.
#
# Every path is read from an environment variable so the scripts run unchanged on another machine. The defaults
# are the paths on the machine the figure was made on.
#
#   ED6_REPO     this repository                        (default: the parent of this file's directory)
#   ED6_DATA     analysis working directory             where intermediate results are written and re-read
#   ED6_REMAP    ReMap 2022 non-redundant peaks, bgzipped and tabix-indexed (hg38)
#   ED6_OUT      where the finished figure is written
#
# Large public reference files (ReMap, Roadmap tagAlign) are NOT in this repository. See README.md.
f5_path <- function(var, default) { v <- Sys.getenv(var, ""); if (nzchar(v)) v else default }
ED6_REPO <- f5_path("ED6_REPO", normalizePath(file.path(
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), mustWork = FALSE))
ED6_DATA     <- f5_path("ED6_DATA",     "/Volumes/Jingjing_Chen/references/controls")
ED6_REMAP    <- f5_path("ED6_REMAP",    "/Volumes/Jingjing_Chen/remap2022_nr_macs2_hg38_v1_0.bed.gz.bgz")
ED6_OUT      <- f5_path("ED6_OUT",      file.path(ED6_REPO, "figures"))
dir.create(ED6_OUT, showWarnings = FALSE, recursive = TRUE)
repo <- function(...) file.path(ED6_REPO, ...)
wd   <- function(...) file.path(ED6_DATA, ...)

`%||%`   <- function(a, b) if (is.null(a) || length(a) == 0) b else a
msg      <- function(...) cat("[INFO]", sprintf(...), "\n")
warn_msg <- function(...) cat("[WARN]", sprintf(...), "\n")
read_set <- function(f) toupper(trimws(readLines(f)))

RED <- "#C6553B"; BLUE <- "#3B6FA8"; INK <- "#2B2B2B"; MUTED <- "#6B6B6B"
GREY <- "#C9C9C9"; AMBER <- "#D9A441"

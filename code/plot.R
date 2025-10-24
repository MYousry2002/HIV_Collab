#!/usr/bin/env Rscript
# plot.R — clustered heatmap for count matrix (rows = sequences, cols = TF families)

# Rscript /projectnb/vcres/myousry/HIV_Colab/code/plot.R   --input ../results/fimo_hiv2_clean/counts.tsv   --output ../results/fimo_hiv2_clean/heatmap.png   --cmap viridis   --figwidth 35   --figheight 50   --dpi 300   --pdf --exclude_col ZNF --cellheight 8
# Rscript /projectnb/vcres/myousry/HIV_Colab/code/plot.R   --input ../results/fimo_hiv1_clean/counts.tsv   --output ../results/fimo_hiv1_clean/heatmap.png   --cmap viridis   --figwidth 30   --figheight 100   --dpi 300   --pdf --exclude_col ZNF --cellheight 8

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(gtable)
})

# ---------- helpers ----------
zscore_rows <- function(m) {
  rsd <- matrixStats::rowSds(m, na.rm=TRUE)
  rsd[rsd == 0] <- NA
  sweep(sweep(m, 1, rowMeans(m, na.rm=TRUE), "-"), 1, rsd, "/")
}
zscore_cols <- function(m) {
  csd <- matrixStats::colSds(m, na.rm=TRUE)
  csd[csd == 0] <- NA
  sweep(sweep(m, 2, colMeans(m, na.rm=TRUE), "-"), 2, csd, "/")
}

clip_percentile <- function(m, p) {
  # symmetric clipping for z-scored data; otherwise clip both tails by percentile
  lo <- quantile(m, probs = (100 - p)/100, na.rm=TRUE)
  hi <- quantile(m, probs = p/100, na.rm=TRUE)
  m[m < lo] <- lo
  m[m > hi] <- hi
  m
}

make_palette <- function(name) {
  # basic named palettes; fall back to viridis-like if unknown
  nm <- tolower(name)
  if (nm %in% c("vlag","viko","coolwarm","rdblu","rdylbu","bwr","blue-red")) {
    colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
  } else if (nm %in% c("magma","inferno","plasma","viridis")) {
    # approximate with Brewer if viridisLite not installed
    if (!requireNamespace("viridisLite", quietly = TRUE)) {
      colorRampPalette(brewer.pal(9, "YlGnBu"))(255)
    } else {
      switch(nm,
        magma   = viridisLite::magma(255),
        inferno = viridisLite::inferno(255),
        plasma  = viridisLite::plasma(255),
        viridis = viridisLite::viridis(255)
      )
    }
  } else {
    # default
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      viridisLite::viridis(255)
    } else {
      colorRampPalette(brewer.pal(9, "YlGnBu"))(255)
    }
  }
}

# ---------- args ----------
opt_list <- list(
  make_option(c("--input"), type="character", help="Input TSV; first column 'sequence_name'"),
  make_option(c("--output"), type="character", help="Output PNG path"),
  make_option(c("--pdf"), action="store_true", default=FALSE, help="Also save PDF"),
  make_option(c("--log1p"), action="store_true", default=FALSE, help="Apply log1p"),
  make_option(c("--zscore"), type="character", default="none",
              help="Z-score: none|rows|cols"),
  make_option(c("--metric"), type="character", default="correlation",
              help="Distance metric: correlation|euclidean|cosine"),
  make_option(c("--method"), type="character", default="average",
              help="Linkage: average|complete|single|ward.D|ward.D2"),
  make_option(c("--cmap"), type="character", default="viridis",
              help="Colormap name (viridis, magma, RdBu/vlag, etc.)"),
  make_option(c("--figwidth"), type="double", default=20, help="Figure width (inches)"),
  make_option(c("--figheight"), type="double", default=22, help="Figure height (inches)"),
  make_option(c("--dpi"), type="integer", default=300, help="DPI for PNG"),
  make_option(c("--xtick_rot"), type="integer", default=75, help="Column label angle"),
  make_option(c("--ytick_size"), type="double", default=7.5, help="Row label font size"),
  make_option(c("--cellheight"), type="double", default=NA_real_,
              help="Row cell height in pixels for pheatmap (NA = auto)"),
  make_option(c("--min_cellheight"), type="double", default=12,
              help="If --cellheight is NA, auto-compute but not below this many pixels per row"),
  make_option(c("--treeheight_row"), type="integer", default=120, help="Row dendrogram height"),
  make_option(c("--treeheight_col"), type="integer", default=60, help="Column dendrogram height"),
  make_option(c("--clip_percentile"), type="double", default=NA_real_,
              help="Clip values to [100-p, p] (e.g., 99)."),
  make_option(c("--exclude_col"), type="character", default=NULL,
              help="Column name to exclude from heatmap (e.g., ZNF)")
  ,make_option(c("--norm"), type="character", default="none",
              help="Normalization: none|row_prop|col_prop|row_cpm|col_cpm|tfidf|row_l2"),
  make_option(c("--max_row_labels"), type="integer", default=0,
              help="Maximum number of row labels to draw (0 = no limit; labels are evenly thinned if >0)"),
  make_option(c("--label_every"), type="integer", default=1,
              help="Draw every Nth row label (1 = show all; overrides --max_row_labels if set)")
  ,make_option(c("--pad_bottom"), type="double", default=18,
              help="Extra bottom padding in points to prevent label clipping (default: 18)")
  ,make_option(c("--pad_right"), type="double", default=12,
              help="Extra right padding in points to prevent clipping (default: 12)")
)
opt <- parse_args(OptionParser(option_list = opt_list))

if (is.null(opt$input) || is.null(opt$output)) {
  stop("Provide --input and --output.")
}

# ---------- load ----------
dt <- data.table::fread(opt$input, sep="\t", header=TRUE, data.table=FALSE, check.names=FALSE)
if (ncol(dt) < 2) stop("Input must have at least two columns.")
# Normalize first column name and ensure it's 'sequence_name'
first_col <- trimws(names(dt)[1])
if (tolower(first_col) != "sequence_name") {
  message("[warn] first column renamed to 'sequence_name' (was: ", first_col, ")")
}
names(dt)[1] <- "sequence_name"
# Use the real sequence IDs as row names
rownames(dt) <- as.character(dt$sequence_name)
dt$sequence_name <- NULL

# ensure numeric (preserve row names)
mat <- as.matrix(data.frame(lapply(dt, function(x) as.numeric(as.character(x))), check.names = FALSE))
rownames(mat) <- rownames(dt)

# exclude specified column if provided
if (!is.null(opt$exclude_col) && opt$exclude_col %in% colnames(mat)) {
  message("[info] excluding column: ", opt$exclude_col)
  mat <- mat[, setdiff(colnames(mat), opt$exclude_col), drop=FALSE]
}

# ---- normalization (before log1p/zscore) ----
normalize_matrix <- function(m, method) {
  method <- tolower(method)
  if (method == "none") return(m)
  if (method == "row_prop") {
    rs <- rowSums(m, na.rm=TRUE); rs[rs == 0] <- NA
    return(sweep(m, 1, rs, "/"))
  }
  if (method == "col_prop") {
    cs <- colSums(m, na.rm=TRUE); cs[cs == 0] <- NA
    return(sweep(m, 2, cs, "/"))
  }
  if (method == "row_cpm") {
    rs <- rowSums(m, na.rm=TRUE); rs[rs == 0] <- NA
    return(sweep(m, 1, rs, "/") * 1e6)
  }
  if (method == "col_cpm") {
    cs <- colSums(m, na.rm=TRUE); cs[cs == 0] <- NA
    return(sweep(m, 2, cs, "/") * 1e6)
  }
  if (method == "row_l2") {
    l2 <- sqrt(rowSums(m^2, na.rm=TRUE)); l2[l2 == 0] <- NA
    return(sweep(m, 1, l2, "/"))
  }
  if (method == "tfidf") {
    # TF: per-row proportions
    tf_den <- rowSums(m, na.rm=TRUE); tf_den[tf_den == 0] <- NA
    tf <- sweep(m, 1, tf_den, "/")
    # IDF: downweight ubiquitous columns
    df <- colSums(m > 0, na.rm=TRUE); N <- nrow(m)
    idf <- log(1 + N / (1 + df))
    tfidf <- sweep(tf, 2, idf, "*")
    # optional L2 per row for stability
    l2 <- sqrt(rowSums(tfidf^2, na.rm=TRUE)); l2[l2 == 0] <- NA
    tfidf <- sweep(tfidf, 1, l2, "/")
    return(tfidf)
  }
  stop(paste0("Unknown --norm: ", method))
}
if (!is.null(opt$norm) && tolower(opt$norm) != "none") {
  message("[info] normalization: ", opt$norm)
}
mat <- normalize_matrix(mat, opt$norm)

# transform
if (isTRUE(opt$log1p)) {
  mat <- log1p(mat)
}

# zscore
if (opt$zscore == "rows") {
  if (!requireNamespace("matrixStats", quietly=TRUE)) stop("Install matrixStats for row/col z-scoring.")
  mat <- zscore_rows(mat)
} else if (opt$zscore == "cols") {
  if (!requireNamespace("matrixStats", quietly=TRUE)) stop("Install matrixStats for row/col z-scoring.")
  mat <- zscore_cols(mat)
}

# drop all-NA rows/cols and zero-variance
row_var <- apply(mat, 1, function(v) var(v, na.rm=TRUE))
col_var <- apply(mat, 2, function(v) var(v, na.rm=TRUE))
keep_r <- !is.na(row_var) & row_var > 0
keep_c <- !is.na(col_var) & col_var > 0
mat <- mat[keep_r, keep_c, drop=FALSE]
mat <- mat[ rowSums(is.na(mat)) < ncol(mat), , drop=FALSE]
mat <- mat[ , colSums(is.na(mat)) < nrow(mat), drop=FALSE]
if (nrow(mat) == 0 || ncol(mat) == 0) stop("Matrix empty after filtering.")

# clipping
if (!is.na(opt$clip_percentile)) {
  mat <- clip_percentile(mat, opt$clip_percentile)
}

# --- order columns by average presence (descending) ---
# Presence is proxied by column mean after normalization/log/zscore/clipping above
if (ncol(mat) > 1) {
  col_means <- colMeans(mat, na.rm = TRUE)
  ord_cols  <- order(col_means, decreasing = TRUE, na.last = NA)
  mat <- mat[, ord_cols, drop = FALSE]
}

# --- build row label vector (thin if too many) ---
if (is.null(rownames(mat))) {
  message("[warn] input had no row names; synthesizing from row index")
  rownames(mat) <- paste0("row_", seq_len(nrow(mat)))
}
all_row_names <- rownames(mat)
n_rows <- nrow(mat)

# Determine thinning step: prefer explicit --label_every, else derive from --max_row_labels
lab_every_input <- suppressWarnings(as.integer(opt$label_every))
if (!is.na(lab_every_input) && lab_every_input > 0L) {
  lab_every <- lab_every_input
} else if (!is.na(opt$max_row_labels) && opt$max_row_labels > 0L) {
  lab_every <- ceiling(n_rows / as.integer(opt$max_row_labels))
  if (is.na(lab_every) || lab_every < 1L) lab_every <- 1L
} else {
  lab_every <- 1L
}
if (n_rows <= 1L) lab_every <- 1L

show_idx <- rep(FALSE, n_rows)
if (n_rows >= 1L) show_idx[seq.int(1L, n_rows, by = lab_every)] <- TRUE
row_labels <- ifelse(show_idx, all_row_names, "")
if (lab_every > 1L) {
  message("[info] thinning row labels: showing 1 every ", lab_every, " rows (", sum(show_idx), "/", n_rows, ").")
}

# distances
dist_fun <- function(x, metric) {
  if (metric == "correlation") {
    as.dist(1 - stats::cor(t(x), use="pairwise.complete.obs", method="pearson"))
  } else if (metric == "cosine") {
    # cosine distance = 1 - cosine similarity (rows are observations)
    xx <- x
    norms <- sqrt(rowSums(xx^2, na.rm=TRUE))
    norms[norms == 0] <- NA
    xx <- sweep(xx, 1, norms, "/")
    # similarity between rows
    sim <- xx %*% t(xx)
    sim[is.na(sim)] <- 0
    as.dist(1 - sim)
  } else {
    dist(x, method = metric)
  }
}

row_dist <- if (nrow(mat) > 1) dist_fun(mat, opt$metric) else NULL
col_dist <- NULL
col_hclust <- NULL

row_hclust <- if (!is.null(row_dist)) hclust(row_dist, method = opt$method) else NA

# colors & label params
pal <- make_palette(opt$cmap)

# pheatmap only accepts specific angles for column labels; coerce to nearest allowed
angle_allowed <- c(0, 45, 90, 270, 315)
xtick_rot_num <- suppressWarnings(as.numeric(opt$xtick_rot))
if (is.na(xtick_rot_num)) xtick_rot_num <- 45
# choose nearest angle (simple absolute difference is fine here)
angle_col_val <- as.character(angle_allowed[which.min(abs(xtick_rot_num - angle_allowed))])

# Compute cellheight if not provided: allocate pixels based on figure height and dendrogram
cellheight_to_use <- opt$cellheight
if (is.na(cellheight_to_use)) {
  total_px <- opt$figheight * opt$dpi
  # subtract row dendrogram + margins (~80 px safety)
  avail_px <- max(100, total_px - (opt$treeheight_row + 80))
  cellheight_to_use <- max(opt$min_cellheight, floor(avail_px / max(1, nrow(mat))))
}

ph <- pheatmap(mat,
               color = pal,
               cluster_rows = if (!is.null(row_dist)) row_hclust else FALSE,
               cluster_cols = FALSE,
               clustering_distance_rows = opt$metric,
               clustering_distance_cols = opt$metric,
               clustering_method = opt$method,
               treeheight_row = opt$treeheight_row,
               treeheight_col = opt$treeheight_col,
               angle_col = angle_col_val,
               fontsize_row = opt$ytick_size,
               fontsize_col = 12,
               cellheight = cellheight_to_use,
               show_rownames = TRUE,
               labels_row = row_labels,
               show_colnames = TRUE,
               border_color = NA,
               silent = TRUE)
gt <- ph$gtable
if (!is.null(opt$pad_bottom) && opt$pad_bottom > 0) {
  gt <- gtable::gtable_add_rows(gt, grid::unit(opt$pad_bottom, "pt"), pos = nrow(gt))
}
if (!is.null(opt$pad_right) && opt$pad_right > 0) {
  gt <- gtable::gtable_add_cols(gt, grid::unit(opt$pad_right, "pt"), pos = ncol(gt))
}
png(filename = opt$output, width = opt$figwidth, height = opt$figheight, units = "in", res = opt$dpi)
grid::grid.newpage(); grid::grid.draw(gt); dev.off()

if (isTRUE(opt$pdf)) {
  pdf_file <- sub("\\.png$", ".pdf", opt$output)
  pdf(file = pdf_file, width = opt$figwidth, height = opt$figheight, useDingbats = FALSE, paper = "special")
  grid::grid.newpage(); grid::grid.draw(gt); dev.off()
}

message("[ok] saved: ", opt$output)
if (isTRUE(opt$pdf)) message("[ok] saved: ", sub("\\.png$", ".pdf", opt$output))
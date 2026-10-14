#!/usr/bin/env Rscript

# 01_build_merged_table.R  (Step 1 of 3: Build Complete Master Merged Table)
# -----------------------------------------------------------------------------
# Merges genome_overview_qc_*.tsv + amr_reconciled_*.tsv + optional .treefile
# into a complete master CSV (all genes & raw ST/Serotype preserved).
# Also writes a <out>.genes_catalog.tsv helper file so you can easily copy-paste
# gene names into a .txt file for 02_filter_merged_table.R.
#
# Usage:
#   Rscript 01_build_merged_table.R \
#       --overview .\genome_overview_qc_TN.tsv \
#       --amr      .\amr_reconciled_TN.tsv \
#       --tree     .\core_genome.treefile \
#       --out      .\merged_Ecoli_TN_full.csv
# -----------------------------------------------------------------------------

options(warn = 1)

log_msg <- function(...) {
  cat(paste0(..., "\n"))
  flush.console()
}

parse_args <- function(args) {
  opt <- list(
    overview     = NULL,
    amr          = NULL,
    tree         = NULL,
    gene_label   = "auto",
    sample_col   = NULL,
    st_col       = NULL,
    serotype_col = NULL,
    out          = "merged_tables_full.csv"
  )
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    take <- function() {
      if (i + 1 > length(args)) stop("Argument ", a, " requires a value.", call. = FALSE)
      i <<- i + 1
      args[i]
    }
    if      (a %in% c("--overview", "--metadata", "--accession", "-m", "-a")) opt$overview <- take()
    else if (a %in% c("--amr", "--rgi", "-r")) opt$amr          <- take()
    else if (a %in% c("--tree", "-t"))         opt$tree         <- take()
    else if (a == "--gene-label")              opt$gene_label   <- take()
    else if (a == "--sample-col")              opt$sample_col   <- take()
    else if (a == "--st-col")                  opt$st_col       <- take()
    else if (a == "--serotype-col")            opt$serotype_col <- take()
    else if (a %in% c("--out", "-o"))          opt$out          <- take()
    else if (a %in% c("--help", "-h")) {
      log_msg("Usage: Rscript 01_build_merged_table.R --overview <TSV> --amr <TSV> [--tree <NWK>] --out <CSV>")
      quit(save = "no", status = 0)
    }
    else stop("Unknown argument: '", a, "'.", call. = FALSE)
    i <- i + 1
  }
  opt
}

read_collate_table <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- sub("\r$", "", lines)
  lines <- lines[!grepl("^\\s*#", lines) & nzchar(trimws(lines))]
  if (length(lines) < 2) stop("File has no header/data rows: ", path, call. = FALSE)

  if (grepl("\t", lines[1])) {
    split_lines <- strsplit(lines, "\t", fixed = TRUE)
    hdr <- split_lines[[1]]
    max_cols <- length(hdr)
    body <- lapply(split_lines[-1], function(x) {
      if (length(x) < max_cols) c(x, rep("", max_cols - length(x)))
      else x[seq_len(max_cols)]
    })
    df <- as.data.frame(do.call(rbind, body), stringsAsFactors = FALSE)
    colnames(df) <- hdr
    return(df)
  } else {
    con <- textConnection(lines)
    on.exit(close(con))
    read.csv(con, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
  }
}

norm_name <- function(x) {
  x <- as.character(x)
  x <- sub("\\.(gff3?|fasta|fa|fna|tsv|csv)$", "", x, ignore.case = TRUE)
  trimws(x)
}

parse_reconciled_cell <- function(v) {
  v_chr <- trimws(as.character(v))
  v_chr[is.na(v_chr) | v_chr %in% c("", "NA", "N/A", "-")] <- "0,0"
  if (all(v_chr %in% c("0", "1", "2"))) return(v_chr)

  is_xy <- grepl("^[0-9NA]+,[0-9NA]+$", v_chr)
  if (any(is_xy)) {
    parts <- strsplit(v_chr, ",", fixed = TRUE)
    vapply(parts, function(p) {
      x <- suppressWarnings(as.numeric(p[1]))
      y <- suppressWarnings(as.numeric(p[2]))
      if (is.na(x)) x <- 0
      if (is.na(y)) y <- 0
      if (x > 0 && y > 0) "2"
      else if (x > 0 || y > 0) "1"
      else "0"
    }, character(1))
  } else {
    ifelse(v_chr %in% c("0", "0.0", "0,0", "Absent", "absent", "FALSE", "no"), "0", "2")
  }
}

tryCatch({
  log_msg("=== [01] Starting Master Table Build ===")
  opt <- parse_args(commandArgs(trailingOnly = TRUE))

  if (is.null(opt$overview) || is.null(opt$amr)) {
    stop("Both --overview and --amr are required!", call. = FALSE)
  }
  for (nm in c("overview", "amr", "tree")) {
    p <- opt[[nm]]
    if (!is.null(p) && !file.exists(p)) stop(sprintf("File for --%s not found: %s", nm, p), call. = FALSE)
  }

  # 1. Read Overview
  ov_raw <- read_collate_table(opt$overview)
  sample_col <- if (!is.null(opt$sample_col)) opt$sample_col else "Sample"
  st_col     <- if (!is.null(opt$st_col))     opt$st_col     else if ("MLST_ST" %in% names(ov_raw)) "MLST_ST" else "ST"
  sero_col   <- if (!is.null(opt$serotype_col)) opt$serotype_col else if ("Serotype" %in% names(ov_raw)) "Serotype" else "K_locus"

  clean_cat <- function(v) {
    v <- trimws(as.character(v))
    ifelse(is.na(v) | v %in% c("", "NA", "N/A", "-", "None", "null", "unknown", "Unknown"), "Unknown", v)
  }

  st_vec <- clean_cat(ov_raw[[st_col]])
  is_num_st <- st_vec != "Unknown" & !grepl("^ST", st_vec, ignore.case = TRUE) & grepl("^[0-9]+$", st_vec)
  st_vec[is_num_st] <- paste0("ST", st_vec[is_num_st])

  meta_df <- data.frame(
    sample     = norm_name(ov_raw[[sample_col]]),
    ST         = st_vec,
    Serotype   = if (sero_col %in% names(ov_raw)) clean_cat(ov_raw[[sero_col]]) else "Unknown",
    Serotype_O = if ("Serotype_O" %in% names(ov_raw)) clean_cat(ov_raw[["Serotype_O"]]) else "Unknown",
    Serotype_H = if ("Serotype_H" %in% names(ov_raw)) clean_cat(ov_raw[["Serotype_H"]]) else "Unknown",
    Overall_QC = if ("Overall_QC" %in% names(ov_raw)) ov_raw[["Overall_QC"]] else "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  meta_df <- meta_df[!duplicated(meta_df$sample), , drop = FALSE]

  if (!is.null(opt$tree)) {
    if (requireNamespace("ape", quietly = TRUE)) {
      tree_tips <- norm_name(ape::read.tree(opt$tree)$tip.label)
    } else {
      nwk_txt <- paste(readLines(opt$tree, warn = FALSE), collapse = "")
      m <- gregexpr("(?<=\\(|,)[^():,;]+(?=:)", nwk_txt, perl = TRUE)
      tree_tips <- norm_name(regmatches(nwk_txt, m)[[1]])
    }
    meta_df$in_tree <- meta_df$sample %in% tree_tips
  }

  # 2. Read AMR Table
  amr_raw <- read_collate_table(opt$amr)
  if (all(c("gene_key", "abricate_gene", "amrfinder_gene") %in% names(amr_raw))) {
    display_gene <- if (opt$gene_label == "gene_key") {
      amr_raw$gene_key
    } else if (opt$gene_label == "abricate_gene") {
      ifelse(nzchar(amr_raw$abricate_gene), amr_raw$abricate_gene, amr_raw$amrfinder_gene)
    } else {
      ifelse(nzchar(amr_raw$amrfinder_gene), amr_raw$amrfinder_gene,
             ifelse(nzchar(amr_raw$abricate_gene), amr_raw$abricate_gene, amr_raw$gene_key))
    }
    amr_raw$.display_gene <- make.unique(trimws(display_gene))
  } else {
    amr_raw$.display_gene <- make.unique(trimws(amr_raw[[1]]))
  }

  col_names_norm <- norm_name(colnames(amr_raw))
  sample_cols_keep <- colnames(amr_raw)[col_names_norm %in% meta_df$sample]

  amr_mat <- as.matrix(amr_raw[, sample_cols_keep, drop = FALSE])
  rownames(amr_mat) <- amr_raw$.display_gene
  amr_t <- as.data.frame(t(amr_mat), stringsAsFactors = FALSE, check.names = FALSE)
  gene_names_final <- colnames(amr_t)

  for (g in gene_names_final) {
    amr_t[[g]] <- parse_reconciled_cell(amr_t[[g]])
  }
  amr_t$sample <- norm_name(rownames(amr_t))
  rownames(amr_t) <- NULL

  idx_match <- match(meta_df$sample, amr_t$sample)
  amr_aligned <- amr_t[idx_match, gene_names_final, drop = FALSE]
  for (g in gene_names_final) {
    amr_aligned[[g]][is.na(amr_aligned[[g]]) | amr_aligned[[g]] == ""] <- "0"
  }

  merged_out <- cbind(meta_df, amr_aligned, stringsAsFactors = FALSE)

  # 3. Write Master CSV & Gene Catalog Helper TSV
  out_full <- normalizePath(opt$out, winslash = "/", mustWork = FALSE)
  write.csv(merged_out, file = out_full, row.names = FALSE, quote = TRUE)

  # Save a helper catalog of all gene columns so you can easily build your .txt filter list
  catalog_path <- sub("\\.csv$", ".genes_catalog.tsv", out_full, ignore.case = TRUE)
  if (all(c("gene_key", "abricate_gene", "amrfinder_gene") %in% names(amr_raw))) {
    cat_df <- data.frame(
      column_name     = amr_raw$.display_gene,
      gene_key        = amr_raw$gene_key,
      abricate_gene   = amr_raw$abricate_gene,
      amrfinder_gene  = amr_raw$amrfinder_gene,
      type            = if ("amrfinder_type" %in% names(amr_raw)) amr_raw$amrfinder_type else "",
      drug_class      = if ("amrfinder_class" %in% names(amr_raw)) amr_raw$amrfinder_class else "",
      n_concordant    = if ("n_samples_concordant" %in% names(amr_raw)) amr_raw$n_samples_concordant else "",
      stringsAsFactors = FALSE
    )
    write.table(cat_df, file = catalog_path, sep = "\t", row.names = FALSE, quote = FALSE)
    log_msg("  [Catalog] Saved gene reference list to: ", catalog_path)
  }

  log_msg(sprintf("=== [01] SUCCESS! Master table saved: %s (%d samples, %d genes) ===",
                  out_full, nrow(merged_out), length(gene_names_final)))
}, error = function(e) {
  log_msg("\n[FATAL ERROR] ", conditionMessage(e))
  quit(save = "no", status = 1)
})
#!/usr/bin/env Rscript

# 02_filter_merged_table.R  (Step 2 of 3: Filter Gene Columns & Top [X] ST / Serotype)
# -----------------------------------------------------------------------------
# Manipulates the master merged table without overwriting it:
#   1. Filters AMR columns using a .txt inclusion list (--genes-list) and/or
#      minimum isolate prevalence (--min-gene-count).
#   2. Collapses ST into the top [X] most frequent types + "Others" + "Unknown" (--top-st).
#   3. Collapses Serotype into the top [Y] most frequent types + "Others" + "Unknown" (--top-serotype).
#
# Usage:
#   Rscript 02_filter_merged_table.R \
#       --in           .\merged_Ecoli_TN_full.csv \
#       --genes-list   .\selected_genes_TN.txt \
#       --top-st       7 \
#       --top-serotype 8 \
#       --out          .\merged_Ecoli_TN_filtered.csv
#
# Options:
#   --in, -i PATH             Input master CSV from 01_build_merged_table.R [required]
#   --out, -o PATH            Output filtered CSV (will not overwrite --in) [required]
#   --genes-list, -g PATH     .txt file of gene columns to INCLUDE (preserves order)
#   --exclude-genes, -x PATH  .txt file of gene columns to EXCLUDE
#   --min-gene-count INT      Drop genes present in fewer than N isolates [default: 0]
#   --concordant-only-count   When using --min-gene-count, count only concordant ('2') hits
#   --top-st INT              Color only top X common STs; lump rest into "Others" ("Unknown" kept separate)
#   --top-serotype INT        Color only top Y common Serotypes; lump rest into "Others"
#   --serotype-source COL     Which column to use as Serotype ("Serotype" [O:H] or "Serotype_O" [O-antigen])
# -----------------------------------------------------------------------------

options(warn = 1)

log_msg <- function(...) {
  cat(paste0(..., "\n"))
  flush.console()
}

parse_args <- function(args) {
  opt <- list(
    input           = NULL,
    out             = NULL,
    genes_list      = NULL,
    exclude_genes   = NULL,
    min_gene_count  = 0L,
    concordant_only = FALSE,
    top_st          = NA_integer_,
    top_serotype    = NA_integer_,
    serotype_source = "Serotype"
  )
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    take <- function() {
      if (i + 1 > length(args)) stop("Argument ", a, " requires a value.", call. = FALSE)
      i <<- i + 1
      args[i]
    }
    if      (a %in% c("--in", "--input", "--merged", "-i", "-m")) opt$input     <- take()
    else if (a %in% c("--out", "-o"))                opt$out             <- take()
    else if (a %in% c("--genes-list", "-g"))         opt$genes_list      <- take()
    else if (a %in% c("--exclude-genes", "-x"))      opt$exclude_genes   <- take()
    else if (a == "--min-gene-count")                opt$min_gene_count  <- as.integer(take())
    else if (a == "--concordant-only-count")         opt$concordant_only <- TRUE
    else if (a == "--top-st")                        opt$top_st          <- as.integer(take())
    else if (a == "--top-serotype")                  opt$top_serotype    <- as.integer(take())
    else if (a == "--serotype-source")               opt$serotype_source <- take()
    else if (a %in% c("--help", "-h")) {
      log_msg("Usage: Rscript 02_filter_merged_table.R --in <MASTER.csv> [--genes-list <GENES.txt>] [--top-st X] [--top-serotype Y] --out <FILTERED.csv>")
      quit(save = "no", status = 0)
    }
    else stop("Unknown argument: '", a, "'.", call. = FALSE)
    i <- i + 1
  }
  opt
}

read_gene_list_file <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- sub("#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(character(0))
  tokens <- unlist(strsplit(lines, "[,\\t]+"))
  unique(trimws(tokens[nzchar(trimws(tokens))]))
}

# Normalize gene string for fuzzy fallback matching (strips 'bla', hyphens, parentheses, quotes)
 strip_gene_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- sub("^bla", "", x)
  gsub("[^a-z0-9]+", "", x)
}

# Collapse a categorical vector to Top X most common values + "Others" + "Unknown"
collapse_to_top_x <- function(vec, top_x, label_name = "Category") {
  v <- trimws(as.character(vec))
  v[is.na(v) | v %in% c("", "NA", "N/A", "-", "none", "None", "null", "unknown", "Unknown", "Untypeable")] <- "Unknown"

  if (is.na(top_x) || top_x <= 0) return(v)

  known_vals <- v[!v %in% c("Unknown", "Others")]
  if (length(known_vals) == 0) return(v)

  freq_table <- sort(table(known_vals), decreasing = TRUE)
  total_unique_known <- length(freq_table)

  if (total_unique_known <= top_x) {
    log_msg(sprintf("  [%s] Only %d known types present (<= top %d); keeping all.", label_name, total_unique_known, top_x))
    return(v)
  }

  top_names <- names(freq_table)[seq_len(top_x)]
  collapsed <- ifelse(v == "Unknown", "Unknown",
                      ifelse(v %in% top_names, v, "Others"))

  top_summary <- paste(sprintf("%s (n=%d)", top_names, as.integer(freq_table[top_names])), collapse = ", ")
  n_others    <- sum(collapsed == "Others")
  n_unknown   <- sum(collapsed == "Unknown")

  log_msg(sprintf("  [%s] Top %d kept: %s", label_name, top_x, top_summary))
  log_msg(sprintf("  [%s] Grouped %d remaining types into 'Others' (n=%d isolates) | 'Unknown' (n=%d isolates)",
                  label_name, total_unique_known - top_x, n_others, n_unknown))
  collapsed
}

tryCatch({
  log_msg("=== [02] Starting Merged Table Filtering ===")
  opt <- parse_args(commandArgs(trailingOnly = TRUE))

  if (is.null(opt$input) || is.null(opt$out)) {
    stop("Both --in and --out are required!", call. = FALSE)
  }
  if (!file.exists(opt$input)) stop("Input file not found: ", opt$input, call. = FALSE)

  in_full  <- normalizePath(opt$input, winslash = "/", mustWork = TRUE)
  out_full <- normalizePath(opt$out, winslash = "/", mustWork = FALSE)

  if (tolower(in_full) == tolower(out_full)) {
    stop("--out path must be different from --in path so the original master table is not overwritten!", call. = FALSE)
  }

  df <- read.csv(in_full, stringsAsFactors = FALSE, check.names = FALSE)
  known_meta_cols <- c(
    "sample", "ST", "MLST_ST", "Serotype", "Serotype_O", "Serotype_H",
    "Overall_QC", "CheckM_Status", "Assembly_Status", "Size_Status", "in_tree",
    "accession", "strain", "geo_location", "O-type", "O_type", "country", "host"
  )
  meta_cols <- intersect(colnames(df), known_meta_cols)
  gene_cols <- setdiff(colnames(df), known_meta_cols)

  log_msg(sprintf("  [Loaded] %d isolates, %d metadata columns, %d AMR gene columns.",
                  nrow(df), length(meta_cols), length(gene_cols)))

  # 1. Optional: Switch Serotype source (e.g., use Serotype_O instead of full O:H Serotype)
  if (opt$serotype_source != "Serotype" && opt$serotype_source %in% colnames(df)) {
    log_msg(sprintf("  [Serotype] Using '%s' as active Serotype column.", opt$serotype_source))
    df$Serotype <- df[[opt$serotype_source]]
  }

  # 2. Collapse ST and Serotype to Top [X] + "Others" + "Unknown"
  if ("ST" %in% colnames(df)) {
    df$ST <- collapse_to_top_x(df$ST, opt$top_st, "ST")
  }
  if ("Serotype" %in% colnames(df)) {
    df$Serotype <- collapse_to_top_x(df$Serotype, opt$top_serotype, "Serotype")
  }

  # 3. Filter AMR Gene Columns by .txt Inclusion List (--genes-list)
  kept_genes <- gene_cols

  if (!is.null(opt$genes_list)) {
    if (!file.exists(opt$genes_list)) stop("--genes-list file not found: ", opt$genes_list, call. = FALSE)
    wanted <- read_gene_list_file(opt$genes_list)
    if (length(wanted) > 0) {
      cols_lower <- tolower(kept_genes)
      cols_fuzzy <- strip_gene_key(kept_genes)

      matched <- character(0)
      missing <- character(0)

      for (w in wanted) {
        if (w %in% kept_genes) {
          matched <- c(matched, w)
        } else {
          idx <- which(cols_lower == tolower(w))
          if (length(idx) == 0) idx <- which(cols_fuzzy == strip_gene_key(w))
          if (length(idx) > 0) {
            matched <- c(matched, kept_genes[idx[1]])
          } else {
            missing <- c(missing, w)
          }
        }
      }
      matched <- unique(matched)
      if (length(matched) == 0) {
        stop("None of the genes in --genes-list matched the columns in ", opt$input,
             "\n  Available gene columns (first 20): ", paste(head(kept_genes, 20), collapse = ", "), call. = FALSE)
      }
      if (length(missing) > 0) {
        log_msg("  [--genes-list] Note - requested gene(s) not found in table: ", paste(missing, collapse = ", "))
      }
      log_msg(sprintf("  [--genes-list] Filtered gene columns: kept %d / %d genes (ordered as in %s).",
                      length(matched), length(kept_genes), basename(opt$genes_list)))
      kept_genes <- matched
    }
  }

  # 4. Exclude specific genes (--exclude-genes)
  if (!is.null(opt$exclude_genes)) {
    if (!file.exists(opt$exclude_genes)) stop("--exclude-genes file not found: ", opt$exclude_genes, call. = FALSE)
    excl <- read_gene_list_file(opt$exclude_genes)
    if (length(excl) > 0) {
      drop_mask <- (tolower(kept_genes) %in% tolower(excl)) | (strip_gene_key(kept_genes) %in% strip_gene_key(excl))
      log_msg(sprintf("  [--exclude-genes] Removed %d gene(s).", sum(drop_mask)))
      kept_genes <- kept_genes[!drop_mask]
    }
  }

  # 5. Optional minimum isolate count filter (--min-gene-count)
  if (!is.na(opt$min_gene_count) && opt$min_gene_count > 0) {
    pos_values <- if (opt$concordant_only) c("2") else c("1", "2")
    counts <- vapply(kept_genes, function(g) sum(as.character(df[[g]]) %in% pos_values, na.rm = TRUE), integer(1))
    keep_mask <- counts >= opt$min_gene_count
    log_msg(sprintf("  [--min-gene-count %d] Kept %d / %d genes.", opt$min_gene_count, sum(keep_mask), length(kept_genes)))
    kept_genes <- kept_genes[keep_mask]
  }

  df_out <- df[, c(meta_cols, kept_genes), drop = FALSE]
  write.csv(df_out, file = out_full, row.names = FALSE, quote = TRUE)

  log_msg(sprintf("=== [02] SUCCESS! Saved filtered table: %s (%d isolates, %d AMR genes) ===",
                  out_full, nrow(df_out), length(kept_genes)))
}, error = function(e) {
  log_msg("\n[FATAL ERROR] ", conditionMessage(e))
  quit(save = "no", status = 1)
})
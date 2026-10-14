#!/usr/bin/env Rscript

# 03_plot_amr_tree_heatmap.R  (Step 3 of 3: Render Tree + ST + Serotype + AMR Heatmap)
# -----------------------------------------------------------------------------
# Usage:
#   Rscript 03_plot_amr_tree_heatmap.R \
#       --merged     .\merged_Ecoli_TN_filtered.csv \
#       --tree       .\core_genome.treefile \
#       --out-prefix .\Ecoli_TN_Clean_AMR
#
# Key Spacing & Layout Options:
#   --label-offset-frac NUM   Gap from furthest tip to Serotype strip (fraction of tree width) [default: 0.48]
#   --panel-gap-frac NUM      Gap between Serotype strip and AMR heatmap (fraction of tree width) [default: 0.045]
#   --strip-width-mult NUM    Width of ST/Serotype strip relative to 1 AMR column [default: 1.3]
#   --heatmap-width NUM       Total width of AMR heatmap (in tree-width multiples) [default: auto-scaled]
#   --st-as-strip             Draw ST as a 2nd vertical strip next to Serotype
#   --AMR-hmap-contrast       Binary mode: 0 = Absent, 1 or 2 = Present
#
# AMR Legend Tool Name Options:
#   --tool1 TEXT              First AMR tool name [default: "ABRicate"]
#   --tool2 TEXT              Second AMR tool name [default: "AMRFinderPlus"]
#   --amr-lab-0 TEXT          Override label for '0' (Absent)
#   --amr-lab-1 TEXT          Override label for '1' (Single tool)
#   --amr-lab-2 TEXT          Override label for '2' (Reconciled in both)
# -----------------------------------------------------------------------------

options(warn = 1)

log_msg <- function(...) {
  cat(paste0(..., "\n"))
  flush.console()
}

parse_args <- function(args) {
  opt <- list(
    merged            = "merged_tables.csv",
    tree              = NULL,
    out_prefix        = "AMR_tree_heatmap",
    genes_list        = NULL,
    top_st            = NA_integer_,
    top_serotype      = NA_integer_,
    midpoint_root     = TRUE,
    st_as_strip       = FALSE,
    edge_scale        = 3,
    label_offset_frac = 0.48,   # Generous spacing so right-most clade tip labels never touch Serotype strip
    label_offset      = NA_real_,
    panel_gap_frac    = 0.045,  # Clean gap between Serotype strip and AMR heatmap
    strip_width_mult  = 1.3,
    heatmap_width     = NA_real_,
    tip_nudge_frac    = 0.012,
    tip_size          = 3.1,
    tip_point_size    = 2.7,
    gene_font_size    = 3.2,
    cell_border       = "white",
    amr_contrast      = FALSE,
    # Tool names & custom AMR legend labels
    tool1             = "ABRicate",
    tool2             = "AMRFinderPlus",
    amr_lab_0         = NULL,
    amr_lab_1         = NULL,
    amr_lab_2         = NULL,
    # General legend & canvas settings
    legend_pos        = "right",
    legend_box        = "vertical",
    legend_title_size = 11,
    legend_text_size  = 9.5,
    legend_key_size   = 0.42,
    st_ncol           = NA_integer_,
    serotype_ncol     = NA_integer_,
    title_st          = "Sequence Type (ST)",
    title_serotype    = "Serotype",
    title_amr         = "AMR gene",
    fig_width         = NA_real_,
    fig_height        = 18,
    bottom_margin     = 5
  )
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    take <- function() {
      if (i + 1 > length(args)) stop("Argument ", a, " requires a value.", call. = FALSE)
      i <<- i + 1
      args[i]
    }
    if      (a %in% c("--merged", "--in", "-m", "-i")) opt$merged       <- take()
    else if (a %in% c("--tree", "-t"))        opt$tree              <- take()
    else if (a %in% c("--out-prefix", "-o"))  opt$out_prefix        <- take()
    else if (a %in% c("--genes-list", "-g"))  opt$genes_list        <- take()
    else if (a %in% c("--top-st", "--top-n-st")) opt$top_st         <- as.integer(take())
    else if (a %in% c("--top-serotype", "--top-n-serotypes")) opt$top_serotype <- as.integer(take())
    else if (a == "--midpoint-root")          opt$midpoint_root     <- TRUE
    else if (a == "--no-midpoint-root")       opt$midpoint_root     <- FALSE
    else if (a == "--st-as-strip")            opt$st_as_strip       <- TRUE
    else if (a == "--edge-scale")             opt$edge_scale        <- as.numeric(take())
    else if (a == "--label-offset-frac")      opt$label_offset_frac <- as.numeric(take())
    else if (a == "--label-offset")           opt$label_offset      <- as.numeric(take())
    else if (a == "--panel-gap-frac")         opt$panel_gap_frac    <- as.numeric(take())
    else if (a %in% c("--strip-width-mult", "--otype-width-mult")) opt$strip_width_mult <- as.numeric(take())
    else if (a == "--heatmap-width")          opt$heatmap_width     <- as.numeric(take())
    else if (a == "--tip-nudge-frac")         opt$tip_nudge_frac    <- as.numeric(take())
    else if (a == "--tip-size")               opt$tip_size          <- as.numeric(take())
    else if (a == "--tip-point-size")         opt$tip_point_size    <- as.numeric(take())
    else if (a == "--gene-font-size")         opt$gene_font_size    <- as.numeric(take())
    else if (a == "--cell-border")            opt$cell_border       <- take()
    else if (a == "--AMR-hmap-contrast")      opt$amr_contrast      <- TRUE
    else if (a == "--tool1")                  opt$tool1             <- take()
    else if (a == "--tool2")                  opt$tool2             <- take()
    else if (a == "--amr-lab-0")              opt$amr_lab_0         <- take()
    else if (a == "--amr-lab-1")              opt$amr_lab_1         <- take()
    else if (a == "--amr-lab-2")              opt$amr_lab_2         <- take()
    else if (a == "--legend-pos")             opt$legend_pos        <- take()
    else if (a == "--legend-box")             opt$legend_box        <- take()
    else if (a == "--legend-title-size")      opt$legend_title_size <- as.numeric(take())
    else if (a == "--legend-text-size")       opt$legend_text_size  <- as.numeric(take())
    else if (a == "--legend-key-size")        opt$legend_key_size   <- as.numeric(take())
    else if (a %in% c("--st-ncol", "--country-ncol"))     opt$st_ncol       <- as.integer(take())
    else if (a %in% c("--serotype-ncol", "--otype-ncol")) opt$serotype_ncol <- as.integer(take())
    else if (a %in% c("--title-st", "--title-country"))   opt$title_st      <- take()
    else if (a %in% c("--title-serotype", "--title-otype")) opt$title_serotype <- take()
    else if (a == "--title-amr")              opt$title_amr         <- take()
    else if (a == "--fig-width")              opt$fig_width         <- as.numeric(take())
    else if (a == "--fig-height")             opt$fig_height        <- as.numeric(take())
    else if (a == "--bottom-margin")          opt$bottom_margin     <- as.numeric(take())
    else if (a %in% c("--help", "-h")) {
      log_msg("Usage: Rscript 03_plot_amr_tree_heatmap.R --merged <CSV> --tree <NWK> --out-prefix <PREFIX>")
      quit(save = "no", status = 0)
    }
    else stop("Unknown argument: ", a, call. = FALSE)
    i <- i + 1
  }
  opt
}

norm_name <- function(x) {
  x <- as.character(x)
  x <- sub("\\.(gff3?|fasta|fa|fna|tsv|csv)$", "", x, ignore.case = TRUE)
  trimws(x)
}

strip_gene_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- sub("^bla", "", x)
  gsub("[^a-z0-9]+", "", x)
}

collapse_to_top_x <- function(vec, top_x) {
  v <- trimws(as.character(vec))
  v[is.na(v) | v %in% c("", "NA", "N/A", "-", "none", "None", "null", "unknown", "Unknown", "Untypeable")] <- "Unknown"
  if (is.na(top_x) || top_x <= 0) return(v)
  known <- v[!v %in% c("Unknown", "Others")]
  if (length(unique(known)) <= top_x) return(v)
  top_names <- names(sort(table(known), decreasing = TRUE))[seq_len(top_x)]
  ifelse(v == "Unknown", "Unknown", ifelse(v %in% top_names, v, "Others"))
}

order_levels_clean <- function(vec) {
  u <- unique(as.character(vec))
  special <- intersect(c("Others", "Unknown", "unknown", "Untypeable", "-"), u)
  reg <- setdiff(u, special)
  if (length(reg) > 1) {
    freq <- table(vec[vec %in% reg])
    reg <- names(sort(freq, decreasing = TRUE))
  }
  c(reg, special)
}

build_discrete_palette <- function(levels_vec, base_pal) {
  special <- c(
    "Others"     = "#D3D3D3", # Light Grey
    "Unknown"    = "#8C8C8C", # Medium-Dark Grey
    "unknown"    = "#8C8C8C",
    "Untypeable" = "#BDBDBD",
    "-"          = "#CCCCCC"
  )
  regular <- setdiff(levels_vec, names(special))
  reg_cols <- if (length(regular) == 0) {
    character(0)
  } else if (length(regular) <= length(base_pal)) {
    base_pal[seq_along(regular)]
  } else {
    colorRampPalette(base_pal)(length(regular))
  }
  cols <- setNames(reg_cols, regular)
  for (sp in intersect(levels_vec, names(special))) {
    cols[sp] <- special[sp]
  }
  cols
}

tryCatch({
  log_msg("=== [03] Starting Tree + ST + Serotype + AMR Plot ===")
  opt <- parse_args(commandArgs(trailingOnly = TRUE))

  suppressPackageStartupMessages({
    library(tidyverse)
    library(ggtree)
    library(treeio)
    library(RColorBrewer)
    library(ape)
    library(ggnewscale)
  })

  if (is.null(opt$tree) || !file.exists(opt$tree)) stop("Tree file not found: ", opt$tree, call. = FALSE)
  if (!file.exists(opt$merged)) stop("Merged CSV file not found: ", opt$merged, call. = FALSE)

  out_pdf <- normalizePath(paste0(opt$out_prefix, ".pdf"), winslash = "/", mustWork = FALSE)
  out_png <- normalizePath(paste0(opt$out_prefix, ".png"), winslash = "/", mustWork = FALSE)

  # 1. Read Merged CSV
  merged_raw <- read.csv(opt$merged, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"sample" %in% names(merged_raw)) merged_raw$sample <- merged_raw[[1]]
  if (!"ST" %in% names(merged_raw))     merged_raw$ST     <- "Unknown"
  if (!"Serotype" %in% names(merged_raw)) merged_raw$Serotype <- "Unknown"

  merged_raw$sample   <- norm_name(merged_raw$sample)
  merged_raw$ST       <- collapse_to_top_x(merged_raw$ST, opt$top_st)
  merged_raw$Serotype <- collapse_to_top_x(merged_raw$Serotype, opt$top_serotype)
  merged_raw          <- merged_raw[!duplicated(merged_raw$sample), , drop = FALSE]

  known_meta_cols <- c(
    "sample", "ST", "MLST_ST", "Serotype", "Serotype_O", "Serotype_H", "Serotype_Plot",
    "Overall_QC", "CheckM_Status", "Assembly_Status", "Size_Status", "in_tree",
    "accession", "strain", "geo_location", "O-type", "O_type", "country", "host"
  )
  amr_cols <- setdiff(names(merged_raw), known_meta_cols)

  if (!is.null(opt$genes_list)) {
    if (!file.exists(opt$genes_list)) stop("--genes-list file not found: ", opt$genes_list, call. = FALSE)
    lines <- readLines(opt$genes_list, warn = FALSE, encoding = "UTF-8")
    lines <- trimws(sub("#.*$", "", lines))
    wanted <- unique(trimws(unlist(strsplit(lines[nzchar(lines)], "[,\\t]+"))))
    if (length(wanted) > 0) {
      cols_lower <- tolower(amr_cols)
      cols_fuzzy <- strip_gene_key(amr_cols)
      matched <- character(0)
      for (w in wanted) {
        if (w %in% amr_cols) matched <- c(matched, w)
        else {
          idx <- which(cols_lower == tolower(w))
          if (length(idx) == 0) idx <- which(cols_fuzzy == strip_gene_key(w))
          if (length(idx) > 0) matched <- c(matched, amr_cols[idx[1]])
        }
      }
      amr_cols <- unique(matched)
    }
  }

  if (length(amr_cols) == 0) stop("No AMR gene columns left to plot!", call. = FALSE)
  amr_wide <- merged_raw[, amr_cols, drop = FALSE]
  rownames(amr_wide) <- merged_raw$sample

  # 2. Load Tree & Align Tips
  tree <- ape::read.tree(opt$tree)
  if (opt$midpoint_root && requireNamespace("phytools", quietly = TRUE)) {
    tree <- tryCatch(phytools::midpoint.root(tree), error = function(e) tree)
  }
  tree$tip.label <- norm_name(tree$tip.label)

  matched_tips <- intersect(tree$tip.label, merged_raw$sample)
  if (length(matched_tips) == 0) stop("No tree tips matched sample names in merged table!", call. = FALSE)
  if (length(matched_tips) < length(tree$tip.label)) {
    tree <- ape::keep.tip(tree, matched_tips)
  }
  amr_wide <- amr_wide[tree$tip.label, , drop = FALSE]

  log_msg(sprintf("  [Plot] %d isolates | %d ST groups | %d Serotype groups | %d AMR gene columns",
                  nrow(amr_wide), length(unique(merged_raw$ST)), length(unique(merged_raw$Serotype)), ncol(amr_wide)))

  # 3. Scale Bar Calculation
  tree_height_raw <- max(ape::node.depth.edgelength(tree), na.rm = TRUE)
  if (tree_height_raw <= 0 || is.infinite(tree_height_raw)) tree_height_raw <- 0.01
  suggested_width <- tree_height_raw * 0.25
  exponent <- floor(log10(suggested_width))
  scalebar_unscaled <- round(suggested_width / (10^exponent)) * (10^exponent)
  if (scalebar_unscaled == 0) scalebar_unscaled <- 10^exponent

  if (!is.null(tree$edge.length) && any(!is.na(tree$edge.length))) {
    tree$edge.length <- tree$edge.length * opt$edge_scale
  }
  scalebar_width <- scalebar_unscaled * opt$edge_scale

  # 4. Prepare Metadata & High-Contrast Palettes
  meta_plot <- merged_raw[match(tree$tip.label, merged_raw$sample), , drop = FALSE]
  st_levels   <- order_levels_clean(meta_plot$ST)
  sero_levels <- order_levels_clean(meta_plot$Serotype)

  meta_plot$ST       <- factor(meta_plot$ST, levels = st_levels)
  meta_plot$Serotype <- factor(meta_plot$Serotype, levels = sero_levels)

  st_colors   <- build_discrete_palette(st_levels, c(brewer.pal(8, "Dark2"), brewer.pal(9, "Set1")))
  sero_colors <- build_discrete_palette(sero_levels, c(brewer.pal(12, "Paired"), brewer.pal(8, "Set2")))

  # 5. Format AMR Matrix
  if (opt$amr_contrast) {
    amr_wide[] <- lapply(amr_wide, function(x) {
      factor(ifelse(as.character(x) %in% c("1", "2"), "1", "0"), levels = c("0", "1"))
    })
  } else {
    amr_wide[] <- lapply(amr_wide, function(x) {
      v <- as.character(x)
      v[!v %in% c("0", "1", "2")] <- "0"
      factor(v, levels = c("0", "1", "2"))
    })
  }

  # 6. Build Base Tree
  n_tips <- length(tree$tip.label)
  tree_plot <- ggtree(tree, size = 0.45) %<+% meta_plot +
    theme_tree() +
    annotate("segment",
             x = 0, xend = scalebar_width,
             y = n_tips + 0.8, yend = n_tips + 0.8,
             color = "black", linewidth = 0.5) +
    annotate("text",
             x = scalebar_width / 2, y = n_tips + 1.3,
             label = format(scalebar_unscaled, scientific = FALSE, trim = TRUE),
             size = 3.2, color = "black", vjust = 0)

  # 7. Proportional Offsets & Auto-Scaled Heatmap Width
  tree_x_range <- diff(range(tree_plot$data$x, na.rm = TRUE))
  if (tree_x_range == 0 || is.infinite(tree_x_range)) tree_x_range <- 1

  heatmap_width <- if (!is.na(opt$heatmap_width)) {
    opt$heatmap_width
  } else {
    max(1.2, min(12, ncol(amr_wide) * 0.12))
  }
  single_col_w <- heatmap_width / max(1, ncol(amr_wide))
  strip_col_w  <- max(0.06, min(0.18, single_col_w * opt$strip_width_mult))

  label_offset <- if (!is.na(opt$label_offset)) opt$label_offset else tree_x_range * opt$label_offset_frac
  panel_gap    <- tree_x_range * opt$panel_gap_frac
  tip_nudge    <- tree_x_range * opt$tip_nudge_frac

  st_ncol   <- if (!is.na(opt$st_ncol))       opt$st_ncol       else if (length(st_levels) > 14) 2L else 1L
  sero_ncol <- if (!is.na(opt$serotype_ncol)) opt$serotype_ncol else if (length(sero_levels) > 14) 2L else 1L
  cell_border <- if (is.na(opt$cell_border) || opt$cell_border %in% c("none", "")) NULL else opt$cell_border

  current_offset <- label_offset
  gheat <- tree_plot

  if (opt$st_as_strip) {
    st_df <- data.frame(ST = meta_plot$ST, row.names = meta_plot$sample, check.names = FALSE)
    gheat <- gheatmap(gheat,
                      st_df,
                      offset            = current_offset,
                      width             = strip_col_w,
                      colnames          = TRUE,
                      colnames_angle    = 90,
                      colnames_offset_y = -0.5,
                      font.size         = opt$gene_font_size,
                      hjust             = 1,
                      color             = cell_border) +
      scale_fill_manual(
        name   = opt$title_st,
        values = st_colors,
        guide  = guide_legend(order = 1, ncol = st_ncol, override.aes = list(shape = NA, color = NA))
      ) +
      ggnewscale::new_scale_fill()

    current_offset <- current_offset + (strip_col_w * tree_x_range) + (panel_gap * 0.5)
  }

  # Layer 1: Serotype Strip
  sero_df <- data.frame(Serotype = meta_plot$Serotype, row.names = meta_plot$sample, check.names = FALSE)
  gheat <- gheatmap(gheat,
                    sero_df,
                    offset            = current_offset,
                    width             = strip_col_w,
                    colnames          = TRUE,
                    colnames_angle    = 90,
                    colnames_offset_y = -0.5,
                    font.size         = opt$gene_font_size,
                    hjust             = 1,
                    color             = cell_border) +
    scale_fill_manual(
      name   = opt$title_serotype,
      values = sero_colors,
      guide  = guide_legend(order = 2, ncol = sero_ncol, override.aes = list(shape = NA, color = NA))
    ) +
    ggnewscale::new_scale_fill()

  amr_offset <- current_offset + (strip_col_w * tree_x_range) + panel_gap

  # Layer 2: AMR Heatmap with Detailed Tool Labels
  lab_0 <- if (!is.null(opt$amr_lab_0)) opt$amr_lab_0 else sprintf("Absent in both (%s & %s)", opt$tool1, opt$tool2)
  lab_1 <- if (!is.null(opt$amr_lab_1)) opt$amr_lab_1 else sprintf("Single tool only (%s or %s)", opt$tool1, opt$tool2)
  lab_2 <- if (!is.null(opt$amr_lab_2)) opt$amr_lab_2 else sprintf("Reconciled in both (%s & %s)", opt$tool1, opt$tool2)

  amr_scale <- if (opt$amr_contrast) {
    scale_fill_manual(
      name   = opt$title_amr,
      values = c("0" = "#E8E8E8", "1" = "#111111"),
      labels = c("0" = lab_0, "1" = sprintf("Present (%s / %s)", opt$tool1, opt$tool2)),
      drop   = FALSE,
      guide  = guide_legend(order = 3, override.aes = list(shape = NA, color = "grey70"))
    )
  } else {
    scale_fill_manual(
      name   = opt$title_amr,
      values = c("0" = "#E8E8E8", "1" = "#888888", "2" = "#111111"),
      labels = c("0" = lab_0, "1" = lab_1, "2" = lab_2),
      drop   = FALSE,
      guide  = guide_legend(order = 3, override.aes = list(shape = NA, color = "grey70"))
    )
  }

  gheat <- gheatmap(gheat,
                    amr_wide,
                    offset            = amr_offset,
                    width             = heatmap_width,
                    colnames_angle    = 90,
                    colnames_offset_y = -0.5,
                    font.size         = opt$gene_font_size,
                    hjust             = 1,
                    color             = cell_border) +
    amr_scale

  # Layer 3: Tip Dots (ST) & Tip Labels
  tip_coords <- gheat$data %>%
    filter(isTip) %>%
    select(node, x, y, label) %>%
    left_join(meta_plot %>% select(sample, ST, Serotype), by = c("label" = "sample"))

  final_plot <- gheat +
    geom_point(data        = tip_coords,
               aes(x = x, y = y, color = ST),
               size        = opt$tip_point_size,
               inherit.aes = FALSE,
               show.legend = !opt$st_as_strip) +
    geom_text(data        = tip_coords,
              aes(x = x, y = y, label = label),
              hjust       = 0,
              nudge_x     = tip_nudge,
              size        = opt$tip_size,
              inherit.aes = FALSE) +
    scale_color_manual(
      name   = opt$title_st,
      values = st_colors,
      guide  = if (opt$st_as_strip) "none" else guide_legend(
        order        = 1,
        ncol         = st_ncol,
        override.aes = list(size = opt$tip_point_size + 1.2, stroke = 0)
      )
    ) +
    coord_cartesian(clip = "off") +
    theme(
      legend.position   = opt$legend_pos,
      legend.box        = opt$legend_box,
      legend.key.height = unit(opt$legend_key_size, "cm"),
      legend.key.width  = unit(opt$legend_key_size, "cm"),
      legend.title      = element_text(size = opt$legend_title_size, face = "bold"),
      legend.text       = element_text(size = opt$legend_text_size),
      legend.spacing.y  = unit(0.25, "cm"),
      plot.margin       = margin(t = 0.5, r = 1.2, b = opt$bottom_margin, l = 0.5, unit = "cm")
    )

  fig_w <- if (!is.na(opt$fig_width)) opt$fig_width else max(15, min(28, 12 + ncol(amr_wide) * 0.17))

  ggsave(out_pdf, final_plot, width = fig_w, height = opt$fig_height, units = "in", limitsize = FALSE)
  ggsave(out_png, final_plot, width = fig_w, height = opt$fig_height, units = "in", dpi = 300, limitsize = FALSE)

  log_msg("==================================================================")
  log_msg("SUCCESS! Saved plots:")
  log_msg("  - PDF : ", out_pdf)
  log_msg("  - PNG : ", out_png)
  log_msg("==================================================================")
}, error = function(e) {
  log_msg("\n[FATAL ERROR] ", conditionMessage(e))
  quit(save = "no", status = 1)
})
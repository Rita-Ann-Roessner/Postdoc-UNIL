# Computing motifs for AF3 benchmark (pLDDT-annotated V/J profiles)
# Rita Ann Roessner (ritaann.roessner@unil.ch)
# 2026-04-21

library(MixTCRviz)
library(ggpubr)

# Save internal MixTCRviz plotVJ before source() may shadow it
internal_plotVJ <- MixTCRviz:::plotVJ

# =============================================================================
# Config — change these for each new dataset
# =============================================================================
motif_file              <- "AF3_class_I_motif.txt"   # pre-filtered input TCRs
baseline_file           <- "AF3_class_I_baseline.txt"                               # NULL → use baseline_HomoSapiens
factor_file             <- "AF3_class_I_factor.csv"
plddt_file              <- "AF3_class_I_cdr_plddts.csv"

plotting_functions_path <- "plotting_functions.R"

use_internal_baseline   <- TRUE    # TRUE: MixTCRviz::baseline_HomoSapiens
                                   # FALSE: use baseline_file as baseline

label_input             <- "Input"
label_baseline          <- "Baseline"

plddt_limits            <- c(0.7, 1.0)
plddt_midpoint          <- 0.85

out_vj_plddt            <- "motif_auc_class_I_plddt.png"
out_vj_plain            <- "motif_auc_class_I_plain.png"
# =============================================================================

source(plotting_functions_path)

# -- Utility ------------------------------------------------------------------

make_info <- function(gene, input_name = "Input", baseline_name = "Baseline",
                      model = "Model_default") {
  x <- c(gene, input_name, baseline_name, model)
  names(x) <- c("gene", "input1.name", "baseline.name", "model")
  x
}

parse_num_list <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return(numeric(0))
  x <- gsub("\\[|\\]", "", x)
  vals <- unlist(strsplit(x, "[,[:space:]]+"))
  vals <- vals[nzchar(vals)]
  as.numeric(vals)
}

# Merges a new weighted count vector into the running total (union of gene names).
accumulate_weighted <- function(running, current) {
  if (is.null(running)) return(current)
  all_names <- union(names(running), names(current))
  v1 <- setNames(rep(0, length(all_names)), all_names)
  v2 <- v1
  v1[names(running)] <- running
  v2[names(current)] <- current
  v1 + v2
}

# -- pLDDT per V/J gene -------------------------------------------------------

compute_plddt_vj <- function(df_plddt) {
  make_plddt <- function(gene_col, cols) {
    cols <- intersect(cols, colnames(df_plddt))
    if (length(cols) == 0 || !gene_col %in% colnames(df_plddt)) return(NULL)
    sapply(split(df_plddt, df_plddt[[gene_col]]), function(d) {
      mean(unlist(lapply(cols, function(col) unlist(lapply(d[[col]], parse_num_list)))),
           na.rm = TRUE)
    })
  }
  list(
    V = list(
      TRA = make_plddt("TRAV", c("CDR1A", "CDR2A", "CDR3AV")),
      TRB = make_plddt("TRBV", c("CDR1B", "CDR2B", "CDR3BV"))
    ),
    J = list(
      TRA = make_plddt("TRAJ", c("CDR3AJ")),
      TRB = make_plddt("TRBJ", c("CDR3BJ"))
    )
  )
}

# -- Weighted V/J aggregation -------------------------------------------------

# Runs MixTCRviz on df_sub, then sums factor-weighted V/J counts across models.
# Iterates over names(motif$stat) so keys always match what MixTCRviz created.
# If df_sub has no 'model' column, falls back to stat[[1]] with no weighting.
build_weighted_stats <- function(df_sub, df_factors) {
  motif      <- MixTCRviz(input1 = df_sub)
  stat_names <- names(motif$stat)
  
  if (!("model" %in% colnames(df_sub)) || length(stat_names) == 1) {
    es <- motif$stat[[1]]
    return(list(total_countV = es$countV, total_countJ = es$countJ))
  }

  total_countV <- list(TRA = NULL, TRB = NULL)
  total_countJ <- list(TRA = NULL, TRB = NULL)

  for (model in stat_names) {
    factor_row <- df_factors$factor[df_factors$model == model]
    factor     <- if (length(factor_row) > 0 && !is.na(factor_row[1])) factor_row[1] else 1
    es         <- motif$stat[[model]]

    for (ch in c("TRA", "TRB")) {
      total_countV[[ch]] <- accumulate_weighted(total_countV[[ch]], es$countV[[ch]] * factor)
      total_countJ[[ch]] <- accumulate_weighted(total_countJ[[ch]], es$countJ[[ch]] * factor)
    }
  }

  list(total_countV = total_countV, total_countJ = total_countJ)
}

# -- V/J plot assembly --------------------------------------------------------

build_and_save_vj <- function(stats_input, rep_countV, rep_countJ,
                               rep_sdV = NULL, rep_sdJ = NULL,
                               plddt_vj, label_input, label_baseline,
                               out_file, use_plddt = TRUE, comp.baseline = TRUE,
                               n.input = NULL, plddt_limits = c(0.5, 1), plddt_midpoint = 0.75) {
  make_vj_info <- function(gene) make_info(gene, input_name  = label_input,
                                                  baseline_name = label_baseline)
  if (use_plddt) {
    pVA <- plotVJ(stats_input$total_countV[["TRA"]], rep_countV[["TRA"]],
                  plddt.es = plddt_vj$V[["TRA"]], sd.rep = rep_sdV[["TRA"]],
                  info = make_vj_info("TRAV"), show.plddt.legend = TRUE,
                  n.input = n.input, plddt_limits = plddt_limits, plddt_midpoint = plddt_midpoint)
    pJA <- plotVJ(stats_input$total_countJ[["TRA"]], rep_countJ[["TRA"]],
                  plddt.es = plddt_vj$J[["TRA"]], sd.rep = rep_sdJ[["TRA"]],
                  info = make_vj_info("TRAJ"),
                  n.input = n.input, plddt_limits = plddt_limits, plddt_midpoint = plddt_midpoint)
    pVB <- plotVJ(stats_input$total_countV[["TRB"]], rep_countV[["TRB"]],
                  plddt.es = plddt_vj$V[["TRB"]], sd.rep = rep_sdV[["TRB"]],
                  info = make_vj_info("TRBV"),
                  n.input = n.input, plddt_limits = plddt_limits, plddt_midpoint = plddt_midpoint)
    pJB <- plotVJ(stats_input$total_countJ[["TRB"]], rep_countJ[["TRB"]],
                  plddt.es = plddt_vj$J[["TRB"]], sd.rep = rep_sdJ[["TRB"]],
                  info = make_vj_info("TRBJ"),
                  n.input = n.input, plddt_limits = plddt_limits, plddt_midpoint = plddt_midpoint)
  } else {
    plain_label <- if (!is.null(n.input)) paste0(label_input, " (", n.input, ")") else label_input
    make_vj_info_plain <- function(gene) make_info(gene, input_name = plain_label,
                                                         baseline_name = label_baseline)
    pVA <- internal_plotVJ(count.es = stats_input$total_countV[["TRA"]],
                           count.rep = rep_countV[["TRA"]], sd.rep = rep_sdV[["TRA"]],
                           info = make_vj_info_plain("TRAV"), comp.baseline = comp.baseline,
                           print.size = FALSE)
    pJA <- internal_plotVJ(count.es = stats_input$total_countJ[["TRA"]],
                           count.rep = rep_countJ[["TRA"]], sd.rep = rep_sdJ[["TRA"]],
                           info = make_vj_info_plain("TRAJ"), comp.baseline = comp.baseline,
                           print.size = FALSE)
    pVB <- internal_plotVJ(count.es = stats_input$total_countV[["TRB"]],
                           count.rep = rep_countV[["TRB"]], sd.rep = rep_sdV[["TRB"]],
                           info = make_vj_info_plain("TRBV"), comp.baseline = comp.baseline,
                           print.size = FALSE)
    pJB <- internal_plotVJ(count.es = stats_input$total_countJ[["TRB"]],
                           count.rep = rep_countJ[["TRB"]], sd.rep = rep_sdJ[["TRB"]],
                           info = make_vj_info_plain("TRBJ"), comp.baseline = comp.baseline,
                           print.size = FALSE)
  }

  p_all <- ggarrange(pVA, pJA, pVB, pJB, ncol = 4, nrow = 1)
  ggsave(out_file, plot = p_all, width = 16, height = 4, dpi = 300)
  invisible(p_all)
}

# =============================================================================
# Main
# =============================================================================

df         <- read.table(motif_file, header = TRUE, sep = "\t")
df_factors <- read.csv(factor_file)
df_plddt   <- read.csv(plddt_file)

plddt_vj    <- compute_plddt_vj(df_plddt)
stats_input <- build_weighted_stats(df, df_factors)

if (use_internal_baseline) {
  bl         <- MixTCRviz::baseline_HomoSapiens
  rep_countV <- bl$countV
  rep_countJ <- bl$countJ
  rep_sdV    <- bl$sdV
  rep_sdJ    <- bl$sdJ
} else {
  df_baseline  <- read.table(baseline_file, header = TRUE, sep = "\t")
  stats_bl     <- build_weighted_stats(df_baseline, df_factors)
  rep_countV   <- stats_bl$total_countV
  rep_countJ   <- stats_bl$total_countJ
  rep_sdV      <- NULL
  rep_sdJ      <- NULL
}

build_and_save_vj(stats_input, rep_countV, rep_countJ, rep_sdV, rep_sdJ,
                  plddt_vj, label_input, label_baseline,
                  out_vj_plddt, use_plddt = TRUE,  comp.baseline = use_internal_baseline,
                  n.input = nrow(df), plddt_limits = plddt_limits, plddt_midpoint = plddt_midpoint)
build_and_save_vj(stats_input, rep_countV, rep_countJ, rep_sdV, rep_sdJ,
                  plddt_vj, label_input, label_baseline,
                  out_vj_plain, use_plddt = FALSE, comp.baseline = use_internal_baseline,
                  n.input = nrow(df))

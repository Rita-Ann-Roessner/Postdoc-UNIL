# Computing motifs for AF3 benchmark (pLDDT-annotated)
# Rita Ann Roessner (ritaann.roessner@unil.ch)
# 2026-04-21

library(MixTCRviz)
library(ggpubr)

# Save internal MixTCRviz plotVJ before source() may shadow it
internal_plotVJ <- MixTCRviz:::plotVJ

# =============================================================================
# Config — change these for each new dataset
# =============================================================================
motif_file              <- "../AF3_class_I_motif.txt"
baseline_file           <- "../AF3_class_I_baseline.txt"
plddt_file              <- "../AF3_class_I_cdr_plddts.csv"

plotting_functions_path <- "plotting_functions.R"

use_internal_baseline   <- TRUE    # TRUE: MixTCRviz::baseline_HomoSapiens
                                   # FALSE: use baseline_file as baseline
renormVJ                <- TRUE    # weight baseline CDR3/length by input P(VJ)

label_input             <- "Input"
label_baseline          <- "Baseline"

plddt_limits            <- c(0.7, 1.0)
plddt_midpoint          <- 0.85

out_full_plddt          <- "../motif_auc_class_I_plddt.png"
out_full_plain          <- "../motif_auc_class_I_plain.png"
out_full_plddt_novj     <- "../motif_auc_class_I_plddt_novj.png"
out_full_plain_novj     <- "../motif_auc_class_I_plain_novj.png"
# =============================================================================

source(plotting_functions_path)

aa.list <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
N.aa    <- length(aa.list)

# -- Utility ------------------------------------------------------------------

make_info <- function(gene, input_name = "Input", baseline_name = "Baseline",
                      model = "Model_default") {
  x <- c(gene, input_name, baseline_name, model)
  names(x) <- c("gene", "input1.name", "baseline.name", "model")
  x
}

make_chain_info <- function(chain, input_name = "Input", baseline_name = "Baseline",
                             model = "Model_default") {
  x <- c(chain, input_name, baseline_name, model)
  names(x) <- c("chain", "input1.name", "baseline.name", "model")
  x
}

parse_num_list <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return(numeric(0))
  x <- gsub("\\[|\\]", "", x)
  vals <- unlist(strsplit(x, "[,[:space:]]+"))
  vals <- vals[nzchar(vals)]
  as.numeric(vals)
}

# Accumulate named vectors (union of names)
accumulate <- function(running, current) {
  if (is.null(running)) return(current)
  all_names <- union(names(running), names(current))
  v1 <- setNames(rep(0, length(all_names)), all_names)
  v2 <- v1
  v1[names(running)] <- running
  v2[names(current)] <- current
  v1 + v2
}

# Accumulate 2D tables (union of row/col names)
accumulate_table <- function(running, current) {
  if (is.null(running)) return(current)
  all_r <- union(rownames(running), rownames(current))
  all_c <- union(colnames(running), colnames(current))
  m <- matrix(0, nrow = length(all_r), ncol = length(all_c), dimnames = list(all_r, all_c))
  m[rownames(running), colnames(running)] <- m[rownames(running), colnames(running)] + running
  m[rownames(current), colnames(current)] <- m[rownames(current), colnames(current)] + current
  m
}

# Accumulate list of 2D tables (countVJ.L: list per length of V×J tables)
accumulate_table_list <- function(running, current) {
  if (is.null(running)) return(current)
  for (lc in union(names(running), names(current))) {
    if (!lc %in% names(running)) {
      running[[lc]] <- current[[lc]]
    } else if (lc %in% names(current)) {
      running[[lc]] <- accumulate_table(running[[lc]], current[[lc]])
    }
  }
  running
}

# Accumulate list of matrices (countCDR3.L: list per length of 20×L matrices)
accumulate_mat_list <- function(running, current) {
  if (is.null(running)) return(current)
  for (lc in union(names(running), names(current))) {
    if (!lc %in% names(running)) {
      running[[lc]] <- current[[lc]]
    } else if (lc %in% names(current)) {
      running[[lc]] <- running[[lc]] + current[[lc]]
    }
  }
  running
}

# -- pLDDT helpers ------------------------------------------------------------

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
    V = list(TRA = make_plddt("TRAV", c("CDR1A", "CDR2A", "CDR3AV")),
             TRB = make_plddt("TRBV", c("CDR1B", "CDR2B", "CDR3BV"))),
    J = list(TRA = make_plddt("TRAJ", c("CDR3AJ")),
             TRB = make_plddt("TRBJ", c("CDR3BJ")))
  )
}

# Mean pLDDT per CDR3 length (for plotLD)
make_plddt_by_length <- function(df_plddt, col, prefix = "L_") {
  vals_list <- lapply(df_plddt[[col]], parse_num_list)
  cdr3_len  <- sapply(vals_list, length)
  cdr3_mean <- sapply(vals_list, function(x) if (length(x) == 0) NA_real_ else mean(x, na.rm = TRUE))
  keep      <- cdr3_len > 0 & !is.na(cdr3_mean)
  out_raw   <- tapply(cdr3_mean[keep], cdr3_len[keep], mean, na.rm = TRUE)
  out       <- as.numeric(out_raw)
  names(out) <- paste0(prefix, names(out_raw))
  out
}

# pLDDT-weighted AA matrices per CDR3 length (for plotCDR3)
make_plddt_cdr3_matrices <- function(df, seq_col, plddt_col, aa_order = aa.list) {
  keep   <- !is.na(df[[seq_col]]) & df[[seq_col]] != "" &
            !is.na(df[[plddt_col]]) & df[[plddt_col]] != ""
  df     <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) return(list())
  seqs   <- as.character(df[[seq_col]])
  plddts <- lapply(df[[plddt_col]], parse_num_list)
  lens   <- nchar(seqs)
  ok     <- mapply(function(s, p) nchar(s) == length(p), seqs, plddts)
  seqs   <- seqs[ok]; plddts <- plddts[ok]; lens <- lens[ok]
  if (length(seqs) == 0) return(list())
  out <- list()
  for (L in sort(unique(lens))) {
    idx     <- which(lens == L)
    sum_mat <- matrix(0, nrow = length(aa_order), ncol = L,
                      dimnames = list(aa_order, paste0("pos", seq_len(L))))
    n_mat   <- sum_mat
    for (i in seq_along(idx)) {
      aa_vec <- strsplit(seqs[idx[i]], "", fixed = TRUE)[[1]]
      p_vec  <- plddts[[idx[i]]]
      for (pos in seq_len(L)) {
        aa <- aa_vec[pos]
        if (!is.na(aa) && aa %in% aa_order && !is.na(p_vec[pos])) {
          sum_mat[aa, pos] <- sum_mat[aa, pos] + p_vec[pos]
          n_mat[aa, pos]   <- n_mat[aa, pos]   + 1
        }
      }
    }
    avg_mat            <- sum_mat
    avg_mat[n_mat > 0] <- sum_mat[n_mat > 0] / n_mat[n_mat > 0]
    out[[paste0("L_", L)]] <- avg_mat
  }
  out
}

# -- Baseline helpers ---------------------------------------------------------

get_baseline_length_for_chain <- function(es, baseline_model, chain_key, info, renormVJ = TRUE) {
  info_out <- info
  if (!renormVJ) {
    return(list(countL.rep = baseline_model$countL[[chain_key]],
                sd.rep     = baseline_model$sdL[[chain_key]],
                info       = info_out))
  }
  has_input_vj     <- length(es$countVJ[[chain_key]]) > 0
  has_baseline_lvj <- !is.null(baseline_model$countL.VJ[[chain_key]])
  if (has_input_vj && has_baseline_lvj) {
    info_out["baseline.name"] <- paste(info_out["baseline.name"], "P(VJ)", sep = " | ")
    bs <- weighted_countL(baseline_model$countL.VJ[[chain_key]], es$countVJ[[chain_key]])
  } else {
    bs <- baseline_model$countL[[chain_key]]
  }
  list(countL.rep = bs, sd.rep = baseline_model$sdL[[chain_key]], info = info_out)
}

get_baseline_cdr3_for_chain <- function(es, baseline_model, chain_key, info, renormVJ = TRUE) {
  info_out <- info
  if (!renormVJ) {
    return(list(countCDR3.rep = baseline_model$countCDR3.L[[chain_key]], info = info_out))
  }
  has_input_vj     <- max(sapply(es$countVJ.L[[chain_key]], length)) > 0
  has_baseline_vjl <- !is.null(baseline_model$countCDR3.VJL[[chain_key]])
  if (has_input_vj && has_baseline_vjl) {
    info_out["baseline.name"] <- paste(info_out["baseline.name"], "P(VJ)", sep = " | ")
    bs <- weighted_countCDR3(baseline_model$countCDR3.VJL[[chain_key]], es$countVJ.L[[chain_key]])
  } else {
    bs <- baseline_model$countCDR3.L[[chain_key]]
  }
  list(countCDR3.rep = bs, info = info_out)
}

# -- Motif aggregation -----------------------------------------------------

build_stats <- function(df_sub) {
  motif      <- MixTCRviz(input1 = df_sub)
  stat_names <- names(motif$stat)
  print('pups')
  print(colnames(df_sub))
  if (!("model" %in% colnames(df_sub)) || length(stat_names) == 1) {
    es <- motif$stat[[1]]
    return(list(total_countV = es$countV, total_countJ = es$countJ,
                countL       = es$countL, countVJ      = es$countVJ,
                countVJ.L    = es$countVJ.L, countCDR3.L = es$countCDR3.L))
  }

  total_countV      <- list(TRA = NULL, TRB = NULL)
  total_countJ      <- list(TRA = NULL, TRB = NULL)
  total_countL      <- list(TRA = NULL, TRB = NULL)
  total_countVJ     <- list(TRA = NULL, TRB = NULL)
  total_countVJ.L   <- list(TRA = NULL, TRB = NULL)
  total_countCDR3.L <- list(TRA = NULL, TRB = NULL)

  for (model in stat_names) {
    es <- motif$stat[[model]]

    for (ch in c("TRA", "TRB")) {
      total_countV[[ch]]      <- accumulate(total_countV[[ch]], es$countV[[ch]])
      total_countJ[[ch]]      <- accumulate(total_countJ[[ch]], es$countJ[[ch]])
      total_countL[[ch]]      <- accumulate(total_countL[[ch]], es$countL[[ch]])
      total_countVJ[[ch]]     <- accumulate_table(total_countVJ[[ch]], es$countVJ[[ch]])
      total_countVJ.L[[ch]]   <- accumulate_table_list(total_countVJ.L[[ch]], es$countVJ.L[[ch]])
      total_countCDR3.L[[ch]] <- accumulate_mat_list(total_countCDR3.L[[ch]], es$countCDR3.L[[ch]])
    }
  }

  list(total_countV = total_countV, total_countJ = total_countJ,
       countL       = total_countL, countVJ      = total_countVJ,
       countVJ.L    = total_countVJ.L, countCDR3.L = total_countCDR3.L)
}

# -- Full figure (VJ + length + CDR3) -----------------------------------------

build_and_save_full <- function(stats_input, baseline_model,
                                 rep_countV, rep_countJ, rep_sdV, rep_sdJ,
                                 plddt_vj, plddt_length, plddt_cdr3,
                                 label_input, label_baseline,
                                 out_file, use_plddt = TRUE,
                                 comp.baseline = TRUE, renormVJ = TRUE,
                                 n.input = NULL,
                                 plddt_limits = c(0.5, 1), plddt_midpoint = 0.75) {

  make_vj_info <- function(gene) make_info(gene, input_name = label_input,
                                                 baseline_name = label_baseline)
  make_ch_info <- function(chain) make_chain_info(chain, input_name = label_input,
                                                        baseline_name = label_baseline)

  # --- V/J plots ---
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
    plain_label        <- if (!is.null(n.input)) paste0(label_input, " (", n.input, ")") else label_input
    make_vj_info_plain <- function(gene) make_info(gene, input_name = plain_label,
                                                        baseline_name = label_baseline)
    pVA <- internal_plotVJ(count.es = stats_input$total_countV[["TRA"]], count.rep = rep_countV[["TRA"]],
                            sd.rep = rep_sdV[["TRA"]], info = make_vj_info_plain("TRAV"),
                            comp.baseline = comp.baseline, print.size = FALSE)
    pJA <- internal_plotVJ(count.es = stats_input$total_countJ[["TRA"]], count.rep = rep_countJ[["TRA"]],
                            sd.rep = rep_sdJ[["TRA"]], info = make_vj_info_plain("TRAJ"),
                            comp.baseline = comp.baseline, print.size = FALSE)
    pVB <- internal_plotVJ(count.es = stats_input$total_countV[["TRB"]], count.rep = rep_countV[["TRB"]],
                            sd.rep = rep_sdV[["TRB"]], info = make_vj_info_plain("TRBV"),
                            comp.baseline = comp.baseline, print.size = FALSE)
    pJB <- internal_plotVJ(count.es = stats_input$total_countJ[["TRB"]], count.rep = rep_countJ[["TRB"]],
                            sd.rep = rep_sdJ[["TRB"]], info = make_vj_info_plain("TRBJ"),
                            comp.baseline = comp.baseline, print.size = FALSE)
  }

  vj_plots <- list(TRA = list(V = pVA, J = pJA), TRB = list(V = pVB, J = pJB))

  # --- Length + CDR3 panels per chain ---
  pg.all <- list()

  for (ch in c("TRA", "TRB")) {
    chain_label <- if (ch == "TRA") "A" else "B"
    info_ch     <- make_ch_info(chain_label)

    # length distribution
    bl_len <- get_baseline_length_for_chain(stats_input, baseline_model, ch, info_ch, renormVJ)
    pLD <- plotLD(
      countL.es      = stats_input$countL[[ch]],
      countL.rep     = bl_len$countL.rep,
      plddtL.es      = if (use_plddt) plddt_length[[ch]] else NULL,
      info           = bl_len$info,
      sd.rep         = bl_len$sd.rep,
      plddt_limits   = plddt_limits,
      plddt_midpoint = plddt_midpoint
    )

    # CDR3 motif
    bl_cdr3 <- get_baseline_cdr3_for_chain(stats_input, baseline_model, ch, info_ch, renormVJ)
    CDR3 <- plotCDR3(
      countL.es      = stats_input$countL[[ch]],
      countL.rep     = baseline_model$countL[[ch]],
      countCDR3.es   = stats_input$countCDR3.L[[ch]],
      countCDR3.rep  = bl_cdr3$countCDR3.rep,
      plddtCDR3.L.es = if (use_plddt) plddt_cdr3[[ch]] else NULL,
      info           = bl_cdr3$info
    )

    cdr3_panel <- ggarrange(CDR3$ES_max, CDR3$Baseline_max, nrow = 2)

    pg.all[[ch]] <- ggarrange(
      vj_plots[[ch]]$V, vj_plots[[ch]]$J,
      pLD, cdr3_panel,
      ncol = 2, nrow = 2
    )
  }

  fig <- ggarrange(pg.all[["TRA"]], pg.all[["TRB"]], ncol = 2)
  ggsave(out_file, plot = fig, width = 20, height = 10, dpi = 300)
  invisible(fig)
}

# =============================================================================
# Main
# =============================================================================

df       <- read.table(motif_file, header = TRUE, sep = "\t")
df_plddt <- read.csv(plddt_file)

plddt_vj    <- compute_plddt_vj(df_plddt)
stats_input <- build_stats(df)

plddt_length <- list(
  TRA = make_plddt_by_length(df_plddt, "CDR3A"),
  TRB = make_plddt_by_length(df_plddt, "CDR3B")
)
plddt_cdr3 <- list(
  TRA = make_plddt_cdr3_matrices(df_plddt, seq_col = "CDR3A_seq", plddt_col = "CDR3A"),
  TRB = make_plddt_cdr3_matrices(df_plddt, seq_col = "CDR3B_seq", plddt_col = "CDR3B")
)

if (use_internal_baseline) {
  baseline_model <- MixTCRviz::baseline_HomoSapiens
  rep_countV     <- baseline_model$countV
  rep_countJ     <- baseline_model$countJ
  rep_sdV        <- baseline_model$sdV
  rep_sdJ        <- baseline_model$sdJ
} else {
  df_baseline    <- read.table(baseline_file, header = TRUE, sep = "\t")
  baseline_model <- MixTCRviz::build_stat(df_baseline, comp.VJL = as.integer(renormVJ))
  rep_countV     <- baseline_model$countV
  rep_countJ     <- baseline_model$countJ
  rep_sdV        <- NULL
  rep_sdJ        <- NULL
}

build_and_save_full(stats_input, baseline_model,
                    rep_countV, rep_countJ, rep_sdV, rep_sdJ,
                    plddt_vj, plddt_length, plddt_cdr3,
                    label_input, label_baseline,
                    out_full_plddt, use_plddt = TRUE, comp.baseline = use_internal_baseline,
                    renormVJ = renormVJ, n.input = nrow(df),
                    plddt_limits = plddt_limits, plddt_midpoint = plddt_midpoint)

build_and_save_full(stats_input, baseline_model,
                    rep_countV, rep_countJ, rep_sdV, rep_sdJ,
                    plddt_vj, plddt_length, plddt_cdr3,
                    label_input, label_baseline,
                    out_full_plain, use_plddt = FALSE, comp.baseline = use_internal_baseline,
                    renormVJ = renormVJ, n.input = nrow(df))

build_and_save_full(stats_input, baseline_model,
                    rep_countV, rep_countJ, rep_sdV, rep_sdJ,
                    plddt_vj, plddt_length, plddt_cdr3,
                    label_input, label_baseline,
                    out_full_plddt_novj, use_plddt = TRUE, comp.baseline = use_internal_baseline,
                    renormVJ = FALSE, n.input = nrow(df),
                    plddt_limits = plddt_limits, plddt_midpoint = plddt_midpoint)

build_and_save_full(stats_input, baseline_model,
                    rep_countV, rep_countJ, rep_sdV, rep_sdJ,
                    plddt_vj, plddt_length, plddt_cdr3,
                    label_input, label_baseline,
                    out_full_plain_novj, use_plddt = FALSE, comp.baseline = use_internal_baseline,
                    renormVJ = FALSE, n.input = nrow(df))

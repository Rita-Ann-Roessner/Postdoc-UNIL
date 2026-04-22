library(MixTCRviz)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Save internal MixTCRviz functions before source() shadows them
internal_plotVJ <- MixTCRviz:::plotVJ
internal_plotLD <- MixTCRviz:::plotLD

source("plotting_functions.R")

# =============================================================================
# Config — change these for each new dataset / system
# =============================================================================
motif_file            <- "B0801_RAKFKQLL_motif.txt"   # pre-filtered model TCRs
baseline_file         <- NULL                            # baseline TCRs (set to NULL when use_internal_baseline=TRUE)
plddt_file            <- "B0801_RAKFKQLL_cdr_plddts.csv"            # per-residue pLDDT values
renormVJ              <- TRUE
use_internal_baseline <- TRUE   # TRUE: use MixTCRviz::baseline_HomoSapiens
                                 # FALSE: use baseline_file as baseline
out_plddt_motif       <- "B0801_RAKFKQLL_motif_plddt.png"               # motif with pLDDT coloring
out_plain_motif       <- "B0801_RAKFKQLL_motif_plain.png"               # standard motif
# =============================================================================

# -- Utility functions --------------------------------------------------------

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

aa.list <- c("A","C","D","E","F","G","H","I","K","L",
             "M","N","P","Q","R","S","T","V","W","Y")
N.aa <- length(aa.list)

make_plddt_cdr3_matrices <- function(df, seq_col, plddt_col, aa_order = aa.list) {
  keep <- !is.na(df[[seq_col]]) & df[[seq_col]] != "" &
    !is.na(df[[plddt_col]]) & df[[plddt_col]] != ""
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) return(list())

  seqs   <- as.character(df[[seq_col]])
  plddts <- lapply(df[[plddt_col]], parse_num_list)
  lens   <- nchar(seqs)
  ok     <- mapply(function(s, p) nchar(s) == length(p), seqs, plddts)
  seqs   <- seqs[ok]; plddts <- plddts[ok]; lens <- lens[ok]
  if (length(seqs) == 0) return(list())

  out <- list()
  for (L in sort(unique(lens))) {
    idx      <- which(lens == L)
    seqs_L   <- seqs[idx]
    plddts_L <- plddts[idx]
    sum_mat  <- matrix(0, nrow = length(aa_order), ncol = L,
                       dimnames = list(aa_order, paste0("pos", seq_len(L))))
    n_mat    <- sum_mat
    for (i in seq_along(seqs_L)) {
      aa_vec <- strsplit(seqs_L[i], "", fixed = TRUE)[[1]]
      p_vec  <- plddts_L[[i]]
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

# -- Analysis functions -------------------------------------------------------

load_data <- function(motif_file, baseline_file = NULL) {
  model    <- read.table(motif_file, header = TRUE)
  baseline <- if (!is.null(baseline_file) && nzchar(baseline_file))
    read.table(baseline_file, header = TRUE) else NULL
  list(model = model, baseline = baseline)
}

# Returns: es, baseline.model, comp.baseline
# - use_internal_baseline=TRUE:  compare against MixTCRviz::baseline_HomoSapiens
# - use_internal_baseline=FALSE: compare against all TCRs in baseline_data
build_stats <- function(model, baseline_data, renormVJ = TRUE,
                        use_internal_baseline = FALSE) {
  if (use_internal_baseline) {
    motif          <- MixTCRviz(input1 = model, plot = FALSE)
    es             <- motif$stat[[1]]
    baseline.model <- MixTCRviz::baseline_HomoSapiens
    comp.baseline  <- TRUE
  } else {
    motif          <- MixTCRviz(input1 = model, input2 = baseline_data,
                                 renormVJ = renormVJ, plot = FALSE)
    es             <- motif$stat[[1]]
    baseline.model <- MixTCRviz::build_stat(baseline_data, comp.VJL = renormVJ)
    comp.baseline  <- FALSE
  }
  list(es = es, baseline.model = baseline.model, comp.baseline = comp.baseline)
}

compute_plddt_vj <- function(es, df_plddt) {
  make_plddt <- function(gene_col, cols) {
    sapply(split(df_plddt, df_plddt[[gene_col]]), function(d) {
      mean(unlist(lapply(cols, function(col) unlist(lapply(d[[col]], parse_num_list)))),
           na.rm = TRUE)
    })
  }
  if (is.null(es$plddtV)) es$plddtV <- list()
  if (is.null(es$plddtJ)) es$plddtJ <- list()
  es$plddtV[["TRA"]] <- make_plddt("TRAV", c("CDR1A", "CDR2A", "CDR3AV"))
  es$plddtV[["TRB"]] <- make_plddt("TRBV", c("CDR1B", "CDR2B", "CDR3BV"))
  es$plddtJ[["TRA"]] <- make_plddt("TRAJ", c("CDR3AJ"))
  es$plddtJ[["TRB"]] <- make_plddt("TRBJ", c("CDR3BJ"))
  es
}

compute_plddt_cdr3 <- function(es, df_plddt) {
  if (is.null(es$plddtCDR3.L)) es$plddtCDR3.L <- list()
  es$plddtCDR3.L[["TRA"]] <- make_plddt_cdr3_matrices(df_plddt, "CDR3A_seq", "CDR3A")
  es$plddtCDR3.L[["TRB"]] <- make_plddt_cdr3_matrices(df_plddt, "CDR3B_seq", "CDR3B")
  es
}

get_baseline_length_for_chain <- function(es, baseline.model, chain_key, info,
                                          renormVJ = TRUE) {
  info_out <- info
  if (!renormVJ) {
    return(list(countL.rep = baseline.model$countL[[chain_key]],
                sd.rep = baseline.model$sdL[[chain_key]], info = info_out))
  }
  has_input_vj     <- length(es$countVJ[[chain_key]]) > 0
  has_baseline_lvj <- !is.null(baseline.model$countL.VJ[[chain_key]])
  if (has_input_vj && has_baseline_lvj) {
    info_out["baseline.name"] <- paste(info_out["baseline.name"], "P(VJ)", sep = " | ")
    bs <- weighted_countL(baseline.model$countL.VJ[[chain_key]], es$countVJ[[chain_key]])
  } else {
    if (!has_input_vj)
      message("No P(VJ) in input1 for chain ", chain_key, " — using renormVJ=FALSE.")
    if (!has_baseline_lvj)
      message("No P(L|VJ) in baseline for chain ", chain_key, " — using renormVJ=FALSE.")
    bs <- baseline.model$countL[[chain_key]]
  }
  list(countL.rep = bs, sd.rep = baseline.model$sdL[[chain_key]], info = info_out)
}

get_baseline_cdr3_for_chain <- function(es, baseline.model, chain_key, info,
                                        renormVJ = TRUE) {
  info_out <- info
  if (!renormVJ) {
    return(list(countCDR3.rep = baseline.model$countCDR3.L[[chain_key]], info = info_out))
  }
  has_input_vj     <- max(sapply(es$countVJ.L[[chain_key]], length)) > 0
  has_baseline_vjl <- !is.null(baseline.model$countCDR3.VJL[[chain_key]])
  if (has_input_vj && has_baseline_vjl) {
    info_out["baseline.name"] <- paste(info_out["baseline.name"], "P(VJ)", sep = " | ")
    bs <- weighted_countCDR3(baseline.model$countCDR3.VJL[[chain_key]],
                             es$countVJ.L[[chain_key]])
  } else {
    if (!has_input_vj)
      message("No P(VJ|L) in input1 for chain ", chain_key, " — using renormVJ=FALSE.")
    if (!has_baseline_vjl)
      message("No P(CDR3|VJL) in baseline for chain ", chain_key, " — using renormVJ=FALSE.")
    bs <- baseline.model$countCDR3.L[[chain_key]]
  }
  list(countCDR3.rep = bs, info = info_out)
}

# use_plddt=TRUE:  custom plotVJ from plotting_functions.R (pLDDT-colored points)
# use_plddt=FALSE: MixTCRviz::plotVJ (standard, with sd from baseline.model)
build_vj_plots <- function(es, baseline.model, comp.baseline, use_plddt = TRUE) {
  if (use_plddt) {
    list(
      pVA = plotVJ(es$countV[["TRA"]], baseline.model$countV[["TRA"]],
                   plddt.es = es$plddtV[["TRA"]], sd.rep = baseline.model$sdV[["TRA"]],
                   info = make_info("TRAV"), species = es$species,
                   show.plddt.legend = TRUE),
      pVB = plotVJ(es$countV[["TRB"]], baseline.model$countV[["TRB"]],
                   plddt.es = es$plddtV[["TRB"]], sd.rep = baseline.model$sdV[["TRB"]],
                   info = make_info("TRBV"), species = es$species),
      pJA = plotVJ(es$countJ[["TRA"]], baseline.model$countJ[["TRA"]],
                   plddt.es = es$plddtJ[["TRA"]], sd.rep = baseline.model$sdJ[["TRA"]],
                   info = make_info("TRAJ"), species = es$species),
      pJB = plotVJ(es$countJ[["TRB"]], baseline.model$countJ[["TRB"]],
                   plddt.es = es$plddtJ[["TRB"]], sd.rep = baseline.model$sdJ[["TRB"]],
                   info = make_info("TRBJ"), species = es$species)
    )
  } else {
    list(
      pVA = internal_plotVJ(count.es = es$countV[["TRA"]],
                             count.rep = baseline.model$countV[["TRA"]],
                             sd.rep = baseline.model$sdV[["TRA"]],
                             info = make_info("TRAV"), species = es$species,
                             comp.baseline = comp.baseline),
      pVB = internal_plotVJ(count.es = es$countV[["TRB"]],
                             count.rep = baseline.model$countV[["TRB"]],
                             sd.rep = baseline.model$sdV[["TRB"]],
                             info = make_info("TRBV"), species = es$species,
                             comp.baseline = comp.baseline),
      pJA = internal_plotVJ(count.es = es$countJ[["TRA"]],
                             count.rep = baseline.model$countJ[["TRA"]],
                             sd.rep = baseline.model$sdJ[["TRA"]],
                             info = make_info("TRAJ"), species = es$species,
                             comp.baseline = comp.baseline),
      pJB = internal_plotVJ(count.es = es$countJ[["TRB"]],
                             count.rep = baseline.model$countJ[["TRB"]],
                             sd.rep = baseline.model$sdJ[["TRB"]],
                             info = make_info("TRBJ"), species = es$species,
                             comp.baseline = comp.baseline)
    )
  }
}

# Always uses MixTCRviz::plotLD (no pLDDT on length plots)
build_length_plots <- function(es, baseline.model, renormVJ, comp.baseline) {
  plots <- list()
  for (ch in c("TRA", "TRB")) {
    info_ch <- make_chain_info(sub("TR", "", ch))
    bl      <- get_baseline_length_for_chain(es, baseline.model, ch, info_ch, renormVJ)
    plots[[ch]] <- internal_plotLD(
      countL.es     = es$countL[[ch]],
      countL.rep    = bl$countL.rep,
      info          = bl$info,
      sd.rep        = bl$sd.rep,
      comp.baseline = comp.baseline
    )
  }
  plots
}

build_cdr3_panels <- function(es, baseline.model, renormVJ, comp.baseline,
                              use_plddt = TRUE) {
  panels <- list()
  for (ch in c("TRA", "TRB")) {
    info_ch <- make_chain_info(sub("TR", "", ch))
    bl      <- get_baseline_cdr3_for_chain(es, baseline.model, ch, info_ch, renormVJ)
    cdr3    <- plotCDR3(
      countL.es      = es$countL[[ch]],
      countL.rep     = baseline.model$countL[[ch]],
      countCDR3.es   = es$countCDR3.L[[ch]],
      countCDR3.rep  = bl$countCDR3.rep,
      plddtCDR3.L.es = if (use_plddt) es$plddtCDR3.L[[ch]] else NULL,
      info           = bl$info,
      comp.baseline  = comp.baseline
    )
    panels[[ch]] <- ggpubr::ggarrange(cdr3$ES_max, cdr3$Baseline_max, nrow = 2)
  }
  panels
}

assemble_figure <- function(vj, length_plots, cdr3_panels) {
  pg <- list()
  for (ch in c("TRA", "TRB")) {
    pV <- if (ch == "TRA") vj$pVA else vj$pVB
    pJ <- if (ch == "TRA") vj$pJA else vj$pJB
    pg[[ch]] <- ggpubr::ggarrange(pV, pJ, length_plots[[ch]], cdr3_panels[[ch]],
                                   ncol = 2, nrow = 2)
  }
  ggpubr::ggarrange(pg[["TRA"]], pg[["TRB"]], ncol = 2)
}

# =============================================================================
# Main
# =============================================================================

dat    <- load_data(motif_file, baseline_file)
stats  <- build_stats(dat$model, dat$baseline, renormVJ, use_internal_baseline)
es             <- stats$es
baseline.model <- stats$baseline.model
comp.baseline  <- stats$comp.baseline

df_plddt <- read.csv(plddt_file)

es <- compute_plddt_vj(es, df_plddt)
es <- compute_plddt_cdr3(es, df_plddt)

length_plots <- build_length_plots(es, baseline.model, renormVJ, comp.baseline)

vj_plddt <- build_vj_plots(es, baseline.model, comp.baseline, use_plddt = TRUE)
vj_plain <- build_vj_plots(es, baseline.model, comp.baseline, use_plddt = FALSE)

cdr3_plddt <- build_cdr3_panels(es, baseline.model, renormVJ, comp.baseline, use_plddt = TRUE)
cdr3_plain <- build_cdr3_panels(es, baseline.model, renormVJ, comp.baseline, use_plddt = FALSE)

fig_plddt <- assemble_figure(vj_plddt, length_plots, cdr3_plddt)
fig_plain <- assemble_figure(vj_plain, length_plots, cdr3_plain)

ggsave(out_plddt_motif, fig_plddt, width = 20, height = 10, dpi = 300)
ggsave(out_plain_motif, fig_plain, width = 20, height = 10, dpi = 300)

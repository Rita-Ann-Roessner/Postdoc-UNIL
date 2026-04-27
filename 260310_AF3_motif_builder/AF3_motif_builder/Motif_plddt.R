library(MixTCRviz)
library(ggplot2)
library(dplyr)
library(ggpubr)

source("plotting_functions.R")

### Helper Functions ###

make_info <- function(gene, input_name = "Input", baseline_name = "Baseline", model = "Model_default") {
  x <- c(gene, input_name, baseline_name, model)
  names(x) <- c("gene", "input1.name", "baseline.name", "model")
  x
}

make_chain_info <- function(chain, input_name = "Input", baseline_name = "Baseline", model = "Model_default") {
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

make_plddt <- function(df_plddt, gene_col, cols) {
  sapply(split(df_plddt, df_plddt[[gene_col]]), function(d) {
    mean(unlist(lapply(cols, function(col) unlist(lapply(d[[col]], parse_num_list)))), na.rm = TRUE)
  })
}

aa.list <- c("A","C","D","E","F","G","H","I","K","L",
             "M","N","P","Q","R","S","T","V","W","Y")
N.aa <- length(aa.list)

make_plddt_cdr3_matrices <- function(df, seq_col, plddt_col, aa_order = aa.list) {
  keep <- !is.na(df[[seq_col]]) & df[[seq_col]] != "" &
    !is.na(df[[plddt_col]]) & df[[plddt_col]] != ""
  df <- df[keep, , drop = FALSE]

  if (nrow(df) == 0) return(list())

  seqs <- as.character(df[[seq_col]])
  plddts <- lapply(df[[plddt_col]], parse_num_list)
  lens <- nchar(seqs)

  ok <- mapply(function(s, p) nchar(s) == length(p), seqs, plddts)
  seqs <- seqs[ok]
  plddts <- plddts[ok]
  lens <- lens[ok]

  if (length(seqs) == 0) return(list())

  out <- list()

  for (L in sort(unique(lens))) {
    idx <- which(lens == L)
    seqs_L <- seqs[idx]
    plddts_L <- plddts[idx]

    sum_mat <- matrix(0, nrow = length(aa_order), ncol = L,
                      dimnames = list(aa_order, paste0("pos", seq_len(L))))
    n_mat <- matrix(0, nrow = length(aa_order), ncol = L,
                    dimnames = list(aa_order, paste0("pos", seq_len(L))))

    for (i in seq_along(seqs_L)) {
      aa_vec <- strsplit(seqs_L[i], "", fixed = TRUE)[[1]]
      p_vec <- plddts_L[[i]]

      for (pos in seq_len(L)) {
        aa <- aa_vec[pos]
        if (!is.na(aa) && aa %in% aa_order && !is.na(p_vec[pos])) {
          sum_mat[aa, pos] <- sum_mat[aa, pos] + p_vec[pos]
          n_mat[aa, pos] <- n_mat[aa, pos] + 1
        }
      }
    }

    avg_mat <- sum_mat
    avg_mat[n_mat > 0] <- sum_mat[n_mat > 0] / n_mat[n_mat > 0]
    avg_mat[n_mat == 0] <- 0

    out[[paste0("L_", L)]] <- avg_mat
  }

  out
}

make_plddt_by_length <- function(df_plddt, col, prefix = "L_") {
  vals_list <- lapply(df_plddt[[col]], parse_num_list)

  cdr3_len <- sapply(vals_list, length)
  cdr3_mean <- sapply(vals_list, function(x) {
    if (length(x) == 0) return(NA_real_)
    mean(x, na.rm = TRUE)
  })

  keep <- cdr3_len > 0 & !is.na(cdr3_mean)
  cdr3_len <- cdr3_len[keep]
  cdr3_mean <- cdr3_mean[keep]

  out_raw <- tapply(cdr3_mean, cdr3_len, mean, na.rm = TRUE)
  out <- as.numeric(out_raw)
  names(out) <- paste0(prefix, names(out_raw))
  out
}

get_baseline_length_for_chain <- function(es, baseline.model, chain_key, info, renormVJ = TRUE) {
  info_out <- info

  if (!renormVJ) {
    return(list(
      countL.rep = baseline.model$countL[[chain_key]],
      sd.rep = baseline.model$sdL[[chain_key]],
      info = info_out
    ))
  }

  has_input_vj <- length(es$countVJ[[chain_key]]) > 0
  has_baseline_lvj <- !is.null(baseline.model$countL.VJ[[chain_key]])

  if (has_input_vj && has_baseline_lvj) {
    info_out["baseline.name"] <- paste(info_out["baseline.name"], "P(VJ)", sep = " | ")
    bs <- weighted_countL(
      baseline.model$countL.VJ[[chain_key]],
      es$countVJ[[chain_key]]
    )
  } else {
    if (!has_input_vj) {
      message("No P(VJ) information in input1 to compute baseline CDR3 ", chain_key,
              " length distribution | P(VJ). Using renormVJ = FALSE behavior.")
    }
    if (!has_baseline_lvj) {
      message("No P(L|VJ) information in baseline/input2 to compute baseline CDR3 ", chain_key,
              " length distribution | P(VJ). Using renormVJ = FALSE behavior.")
    }
    bs <- baseline.model$countL[[chain_key]]
  }

  list(
    countL.rep = bs,
    sd.rep = baseline.model$sdL[[chain_key]],
    info = info_out
  )
}

get_baseline_cdr3_for_chain <- function(es, baseline.model, chain_key, info, renormVJ = TRUE) {
  info_out <- info

  if (!renormVJ) {
    return(list(
      countCDR3.rep = baseline.model$countCDR3.L[[chain_key]],
      info = info_out
    ))
  }

  has_input_vj <- max(sapply(es$countVJ.L[[chain_key]], length)) > 0
  has_baseline_vjl <- !is.null(baseline.model$countCDR3.VJL[[chain_key]])

  if (has_input_vj && has_baseline_vjl) {
    info_out["baseline.name"] <- paste(info_out["baseline.name"], "P(VJ)", sep = " | ")
    bs <- weighted_countCDR3(
      baseline.model$countCDR3.VJL[[chain_key]],
      es$countVJ.L[[chain_key]]
    )
  } else {
    if (!has_input_vj) {
      message("No P(VJ|L) information in input1 to compute baseline CDR3 ", chain_key,
              " motif | P(VJ). Using renormVJ = FALSE behavior.")
    }
    if (!has_baseline_vjl) {
      message("No P(CDR3|VJL) information in baseline/input2 to compute baseline CDR3 ", chain_key,
              " motif | P(VJ). Using renormVJ = FALSE behavior.")
    }
    bs <- baseline.model$countCDR3.L[[chain_key]]
  }

  list(
    countCDR3.rep = bs,
    info = info_out
  )
}


run_pipeline_for_peptide <- function(af3_scores_file, af3_input_file, cdr_plddts_file,
                                     output_plot_file, threshold = 0.8, renormVJ = TRUE) {
  df <- read.table(af3_scores_file, header = TRUE)[, c("id", "AF3_iptm_pair_mean")]
  anno <- read.table(af3_input_file, header = TRUE)

  df <- merge(df, anno, by = "id")

  baseline <- df
  model <- df[df$AF3_iptm_pair_mean > threshold, ]

  motif <- MixTCRviz(input1 = model, input2 = baseline, renormVJ = renormVJ)
  es <- motif$stat[[1]]
  baseline.stat <- MixTCRviz::build_stat(
    baseline,
    chain = "AB",
    species = es$species,
    comp.VJL = 0
  )
  baseline.model <- MixTCRviz::build_stat(baseline, comp.VJL = renormVJ)

  df_plddt <- read.csv(cdr_plddts_file)

  if (is.null(es$plddt)) es$plddt <- list()
  if (is.null(es$plddtV)) es$plddtV <- list()
  if (is.null(es$plddtJ)) es$plddtJ <- list()
  if (is.null(es$plddtL)) es$plddtL <- list()
  if (is.null(es$plddtCDR3.L)) es$plddtCDR3.L <- list()

  es$plddtV[["TRA"]] <- make_plddt(df_plddt, "TRAV", c("CDR1A", "CDR2A", "CDR3AV"))
  es$plddtV[["TRB"]] <- make_plddt(df_plddt, "TRBV", c("CDR1B", "CDR2B", "CDR3BV"))
  es$plddtJ[["TRA"]] <- make_plddt(df_plddt, "TRAJ", c("CDR3AJ"))
  es$plddtJ[["TRB"]] <- make_plddt(df_plddt, "TRBJ", c("CDR3BJ"))

  pVA <- plotVJ(es$countV[["TRA"]], baseline.stat$countV[["TRA"]],
                plddt.es = es$plddtV[["TRA"]],
                info = make_info("TRAV"), species = es$species,
                show.plddt.legend = TRUE)

  pVB <- plotVJ(es$countV[["TRB"]], baseline.stat$countV[["TRB"]],
                plddt.es = es$plddtV[["TRB"]],
                info = make_info("TRBV"), species = es$species)

  pJA <- plotVJ(es$countJ[["TRA"]], baseline.stat$countJ[["TRA"]],
                plddt.es = es$plddtJ[["TRA"]],
                info = make_info("TRAJ"), species = es$species)

  pJB <- plotVJ(es$countJ[["TRB"]], baseline.stat$countJ[["TRB"]],
                plddt.es = es$plddtJ[["TRB"]],
                info = make_info("TRBJ"), species = es$species)

  infoA <- make_chain_info("A")
  infoB <- make_chain_info("B")

  es$plddtL[["TRA"]] <- make_plddt_by_length(df_plddt, "CDR3A")
  es$plddtL[["TRB"]] <- make_plddt_by_length(df_plddt, "CDR3B")

  length_plots <- list()
  length_info <- list(TRA = infoA, TRB = infoB)

  for (ch in names(length_info)) {
    baseline_length <- get_baseline_length_for_chain(
      es = es,
      baseline.model = baseline.model,
      chain_key = ch,
      info = length_info[[ch]],
      renormVJ = renormVJ
    )

    length_plots[[ch]] <- plotLD(
      countL.es = es$countL[[ch]],
      countL.rep = baseline_length$countL.rep,
      plddtL.es = es$plddtL[[ch]],
      info = baseline_length$info,
      sd.rep = baseline_length$sd.rep
    )
  }

  pLDA <- length_plots[["TRA"]]
  pLDB <- length_plots[["TRB"]]

  es$plddtCDR3.L[["TRA"]] <- make_plddt_cdr3_matrices(
    df = df_plddt, seq_col = "CDR3A_seq", plddt_col = "CDR3A"
  )
  es$plddtCDR3.L[["TRB"]] <- make_plddt_cdr3_matrices(
    df = df_plddt, seq_col = "CDR3B_seq", plddt_col = "CDR3B"
  )

  chain_map <- list(
    TRA = list(seq_col = "CDR3A_seq", plddt_col = "CDR3A", info = infoA),
    TRB = list(seq_col = "CDR3B_seq", plddt_col = "CDR3B", info = infoB)
  )

  CDR3 <- list()
  cdr3_panels <- list()

  for (ch in names(chain_map)) {
    info_ch <- chain_map[[ch]]$info

    baseline_cdr3 <- get_baseline_cdr3_for_chain(
      es = es,
      baseline.model = baseline.model,
      chain_key = ch,
      info = info_ch,
      renormVJ = renormVJ
    )

    CDR3[[ch]] <- plotCDR3(
      countL.es = es$countL[[ch]],
      countL.rep = baseline.stat$countL[[ch]],
      countCDR3.es = es$countCDR3.L[[ch]],
      countCDR3.rep = baseline_cdr3$countCDR3.rep,
      plddtCDR3.L.es = es$plddtCDR3.L[[ch]],
      info = baseline_cdr3$info
    )

    cdr3_panels[[ch]] <- ggpubr::ggarrange(
      CDR3[[ch]]$ES_max,
      CDR3[[ch]]$Baseline_max,
      nrow = 2
    )
  }

  pg.all <- list()
  pg.all[["TRA"]] <- ggpubr::ggarrange(pVA, pJA, pLDA, cdr3_panels[["TRA"]], ncol = 2, nrow = 2)
  pg.all[["TRB"]] <- ggpubr::ggarrange(pVB, pJB, pLDB, cdr3_panels[["TRB"]], ncol = 2, nrow = 2)

  fig <- ggpubr::ggarrange(pg.all[["TRA"]], pg.all[["TRB"]], ncol = 2)

  grid::grid.newpage()
  grid::grid.draw(fig)

  ggsave(
    output_plot_file,
    fig,
    width = 20,
    height = 10,
    dpi = 300
  )

  invisible(topdir)
}


### Main ###

base_output_dir <- ".."
threshold <- 0.8 # AF3_iptm_pair_mean 
renormVJ  <- TRUE

dico <- list(
  "LLWNGPMAV" = "HLA_A0201",
  "ELAGIGILTV" = "HLA_A0201",
  "GILGFVFTL"  = "HLA_A0201"
)

for (peptide in names(dico)) {
  topdir <- file.path(base_output_dir, paste0("step2_", peptide))

  af3_scores_file  <- file.path(topdir, paste0("chainA_B_random_pair_", peptide, "_output.txt"))
  af3_input_file   <- file.path(topdir, paste0("chainA_B_random_pair_", peptide, "_input.txt"))
  cdr_plddts_file  <- file.path(topdir, "cdr_plddts.csv")
  output_plot_file <- file.path(topdir, "final_motif.png")

  message("Running pipeline for ", peptide, " / ", dico[[peptide]])
  run_pipeline_for_peptide(
    af3_scores_file  = af3_scores_file,
    af3_input_file   = af3_input_file,
    cdr_plddts_file  = cdr_plddts_file,
    output_plot_file = output_plot_file,
    threshold        = threshold,
    renormVJ         = renormVJ
  )
}
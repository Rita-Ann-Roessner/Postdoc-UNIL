# Extracting specificity encoded in the CDR3s by computing Jensen-Shannon-Divergence between Baseline and Input motif.
# Rita Ann Roessner (ritaann.roessner@unil.ch)
# 2026-04-27

library(MixTCRviz)
library(ggplot2)
library(dplyr)
library(ggpubr)

# =============================================================================
# Input
# =============================================================================
motif_file <- "../test.txt"


# =============================================================================
# Math background
# =============================================================================
#
# We replicate the exact MixTCRviz weblogo pipeline and compute JSD between
# the two logos (input vs VJ-weighted baseline) for CDR3A and CDR3B.
#
# --- Two probability matrices ---
#
#   P  (input logo)    : 20 AAs × L positions, each column sums to 1
#                        = column-normalised raw count matrix from input
#   Q  (baseline logo) : 20 AAs × L positions, each column sums to 1
#                        = VJ-weighted baseline (see below)
#
# --- VJ-weighted baseline (what MixTCRviz actually plots) ---
#
# Instead of using the marginal baseline CDR3 motif, MixTCRviz re-weights the
# baseline CDR3 distribution using the V-J usage *of the input*. For each
# (V, J) pair observed in the input at length L:
#
#   w(V,J) = count_input(V,J,L) / Σ_{v,j} count_input(v,j,L)    [input P(VJ|L)]
#
#   Q_base(V,J,L) = col-normalised baseline CDR3 matrix for that specific (V,J,L)
#
#   Q = col-normalise[ Σ_{V,J}  w(V,J) · Q_base(V,J,L) ]
#
# This controls for the fact that V and J choice strongly determines CDR3
# shape (V anchors the N-term, J anchors the C-term). By projecting the
# baseline onto the input's V-J usage, any remaining divergence (JSD > 0)
# reflects genuine antigen-driven AA selection, not just VJ bias.
#
# --- Length selection ---
#
# lmax = argmax_L count_input(L),  restricted to lengths present in baseline
# (same rule as MixTCRviz: dominant length in the input)
#
# --- JSD at one position p ---
#
#   M_p  = (P_p + Q_p) / 2          [mixture; M_p(aa) ≥ max(P,Q)/2 > 0 always]
#
#   JSD_p = 0.5 · KL(P_p ‖ M_p)  +  0.5 · KL(Q_p ‖ M_p)
#
#   KL(P ‖ M) = Σ_aa  P(aa) · log₂[ P(aa) / M(aa) ]
#
# Because M is the average of P and Q, M(aa) > 0 wherever either has mass,
# so log₂(P/M) is always finite.  Convention: 0·log₂(0) = 0.
# JSD ∈ [0, 1] bits (log₂):  0 = identical, 1 = fully non-overlapping.
#
# --- Scalar per chain ---
#
#   JSD_chain = mean over positions  JSD_p


# =============================================================================
# Helper functions
# =============================================================================

# Column-wise normalisation: raw count matrix → probability matrix (cols sum to 1).
normalize_cols <- function(mat) {
  s <- colSums(mat)
  s[s == 0] <- 1
  sweep(mat, 2, s, "/")
}

kl_div <- function(p, q) {
  idx <- p > 0
  sum(p[idx] * log2(p[idx] / q[idx]))
}

jsd_vec <- function(p, q) {
  m <- (p + q) / 2
  0.5 * kl_div(p, m) + 0.5 * kl_div(q, m)
}

# JSD scalar between two already-normalised 20×L matrices.
jsd_motif <- function(P, Q) {
  mean(vapply(seq_len(ncol(P)), function(j) jsd_vec(P[, j], Q[, j]), numeric(1)))
}

# Most common CDR3 length in the input that is also present in the baseline.
dominant_length <- function(countL_input, countL_baseline) {
  L_both <- intersect(names(countL_input), names(countL_baseline))
  if (length(L_both) == 0) return(NA_character_)
  names(which.max(countL_input[L_both]))    # e.g. "L_12"
}


# =============================================================================
# Main
# =============================================================================

df    <- read.table(motif_file, header = TRUE, sep = "\t")
motif <- MixTCRviz(input1 = df, plot.cdr3.norm=1)
print(motif$plot)

baseline <- MixTCRviz::baseline_HomoSapiens

models <- names(motif$stat)
chains <- c(CDR3A = "TRA", CDR3B = "TRB")

results <- lapply(models, function(model) {
  es  <- motif$stat[[model]]
  row <- list(model = model)
  jsd_vals <- c()

  for (cdr3_name in names(chains)) {
    chain <- chains[[cdr3_name]]

    # --- lmax: most common input length present in baseline ---
    lmax_c <- dominant_length(es$countL[[chain]], baseline$countL[[chain]])
    row[[paste0("best_L_", cdr3_name)]] <- lmax_c

    if (is.na(lmax_c)) {
      warning(sprintf("Model %s | %s: no overlapping CDR3 lengths.", model, chain))
      row[[paste0("JSD_", cdr3_name)]] <- NA
      next
    }

    # --- Input logo: column-normalise raw count matrix ---
    inp_mat <- es$countCDR3.L[[chain]][[lmax_c]]
    if (is.null(inp_mat)) {
      warning(sprintf("Model %s | %s | %s: missing input count matrix.", model, chain, lmax_c))
      row[[paste0("JSD_", cdr3_name)]] <- NA
      next
    }
    P <- normalize_cols(inp_mat)

    # --- Baseline logo: VJ-weighted baseline (replicates MixTCRviz | P(VJ)) ---
    # weighted_countCDR3 uses input's P(VJ|L) to re-weight the baseline CDR3
    # distributions per (V,J,L) compartment, then column-normalises the result.
    bs_all <- MixTCRviz::weighted_countCDR3(
      countCDR3.VJL.baseline = baseline$countCDR3.VJL[[chain]],
      countVJ.L.es            = es$countVJ.L[[chain]]
    )
    Q <- bs_all[[lmax_c]]

    if (is.null(Q) || all(Q == 0)) {
      warning(sprintf("Model %s | %s | %s: empty VJ-weighted baseline matrix.", model, chain, lmax_c))
      row[[paste0("JSD_", cdr3_name)]] <- NA
      next
    }

    jsd <- jsd_motif(P, Q)
    row[[paste0("JSD_", cdr3_name)]] <- jsd
    jsd_vals <- c(jsd_vals, jsd)
  }

  row[["JSD_avg"]] <- if (length(jsd_vals) == 2) mean(jsd_vals) else NA
  as.data.frame(row, stringsAsFactors = FALSE)
})

result_table <- bind_rows(results)
print(result_table)

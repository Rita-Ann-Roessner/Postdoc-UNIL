# Definition of various functions used to prepare input files for computations
# with AF3 or Boltz.
# Some of these functions could however also be useful in other contexts.


#' Getting full V/J sequences from their names.
#'
#' Function to get the part of the full V/J gene sequence that is outside of
#' the CDR3.
#'
#' This function also returns the number of AA that were deleted from the
#' given V/J segment in given TCR sequence.
#'
#' We can either give the chain and a tcrTab that has the columns
#' (TRAV, TRAJ, cdr3_TRA, ... needed based on the asked gene and chain), and
#' we'll transform the gene names found there from e.g. TRAV1-1 to their
#' corresponding full gene sequence, or we can directly give the full gene
#' sequences and we'll extract the part that is before/after the CDR3.
#'
#' @param gene Either "V" or "J".
#' @param chain Tells if it is the *A* or *B* TCR chain. We either need to give
#'  *chain* and *tcrTab* (and possibly *CDR3*) inputs or *geneSeq* and *CDR3*.
#' @param CDR3 (optional) sequence of the CDR3 (otherwise it is taken from
#'    the tcrTab). It can be a vector for all TCRs or a single value if the same
#'    CDR3 is used with all TCRs.
#' @param geneSeq Instead of giving the *tcrTab* we can give directly the
#'    sequence(s) of the genes present for our TCRs (to still determine the
#'    part that was outside of the CDR3 from it).
#' @param nAA_trimMax is the maximum number of AAs that should be removed from
#'    the gene sequence - if more than that, it means that the trimmed part is
#'    likely not really corresponding to CDR3 part but to a bigger part.
#' @param species Indicates from which species the TCRs are (HomoSapiens or
#'    MusMusculus).
#'
#' @return A list (seq, nAA_trimmed) indicating the sequence of the
#'    TRAV/TRAJ/TRBV/TRBJ that is before/after the CDR3 for each TCRs, as well
#'    as the number of AA that were removed during the deletion step from the
#'    TRAV/...
#'
#' @export
get_TRseq <- function(gene, chain=NULL, tcrTab=NULL,
  geneSeq=NULL, CDR3=NULL, nAA_trimMax=15, species="HomoSapiens"){
  if ((is.null(tcrTab) && is.null(geneSeq)) || (!is.null(tcrTab) && !is.null(geneSeq))){
    stop("Should give either tcrTab or geneSeq, but not both.")
  }
  if (!gene %in% c("V", "J")){
    stop("gene should be 'V' or 'J', was ", gene, ".")
  }
  if ((!is.null(chain)) && (!chain %in% c("A", "B"))){
    stop("When given, chain should be either 'A' or 'B'.")
  }
  geneSeq_noCDR3 <- NULL
  # Will be defined below if using a tcrTab but not when giving a geneSeq as input
  if (!is.null(tcrTab)){
    if (is.null(chain)){
      stop("Need to give the 'chain' when using 'tcrTab' as input.")
    }
    if (gene == "V"){
      geneTab <- MixTCRviz::cdr123[[species]][[paste0("TR", chain)]][
        tcrTab[[paste0("TR", chain, gene)]], ]
      # cdr123... is a data.frame with row.names corresponding to TRAV/TRBV
      # segment names, column 'full' giving the full sequence of the TRAV/TRBV
      # and column CDR3 giving the CDR3 AA that are part of the TRxV.
      CDR3pattern <- paste0("(", geneTab$CDR3, ")$")
    } else { # it's J gene
      geneTab <- MixTCRviz::Jseq[[species]][[paste0("TR", chain)]][
        tcrTab[[paste0("TR", chain, gene)]], ]
      CDR3pattern <- paste0("^(", geneTab$CDR3, ")")
    }
    CDR3pattern <- gsub("\\*", "\\\\*", CDR3pattern)
    # Escape the * characters that can be present in some gene CDR3 patterns
    # (these are usually nonFunctional genes that won't really be processed
    # much later on, but otherwise it was making a bug).
    geneSeq <- gsub("g", "", geneTab$full)
    # Remove the gaps that are present in the full sequences (only present for V
    # but won't change the sequence for J so fine to always do it in case some
    # J would also have gaps in some alignments).
    geneSeq_noCDR3 <- stringr::str_replace(geneSeq, pattern=CDR3pattern,
      replacement="")


    if (is.null(CDR3)){
      CDR3 <- tcrTab[[paste0("cdr3_TR", chain)]]
      # Get this CDR3 sequence from the tab only if wasn't given in input
    } else {
      warning("The tcrTab is given in input and CDR3 sequences as well, ",
        "using thus the input CDR3 seq and not the ones from the tcrTab, make ",
        "sure this is expected.")
    }
  }

  if (length(CDR3) != length(geneSeq)){
    if (length(CDR3) == 1){
      CDR3 <- rep(CDR3, times=length(geneSeq))
    } else if (length(geneSeq) == 1){
      geneSeq <- rep(geneSeq, times=length(CDR3))
      geneSeq_noCDR3 <- rep(geneSeq_noCDR3, times=length(CDR3))
    } else {
      stop("geneSeq and CDR3 sequences should be vectors of the same size or ",
        "of size 1.")
    }
  }
  nAA <- nchar(geneSeq)
  if (is.null(geneSeq_noCDR3)){
    # Need to determine this value as input wasn't given as a tcrTab but directly
    # as the full geneSeq.
    geneSeq_noCDR3 <- rep(NA, times=length(geneSeq))
    for (i in 5:1){
      # We'll try to find the begin (or end for J) of the CDR3 in the geneSeq
      # sequence, with as many AAs from the CDR3 as possible. That's why we start
      # by keeping many AAs from CDR3 and decrease at each step.
      # When tcrTab was given as in
      cInds <- is.na(geneSeq_noCDR3)
      # Sequences that couldn't be determined yet.
      if (gene == "V"){
        pat <- paste0("(.*)", stringr::str_sub(CDR3[cInds], end=i))
      } else if (gene == "J"){
        pat <- paste0(stringr::str_sub(CDR3[cInds], start=-i), "(.*)")
      }
      tSeq <- stringr::str_extract(geneSeq[cInds], pattern=pat, group=1)
      # str_extract will keep the AAs from the geneSeq sequence that were matching
      # the 1st parenthesis group ('(.*)'), i.e. the AAs that were before (or
      # after for J) the CDR3 AAs present in the pattern 'pat'.
      tooShort <- is.na(tSeq) | (nAA[cInds] - nchar(tSeq) > nAA_trimMax)
      # Check if some trimmed sequences are too short (the .* will keep the
      # longest matching sequence, but maybe by chance that was for example a
      # 'CAV' in the middle of the Va and the CDR3 part was 'CAE'). In this case
      # we'll need to see if shorter 'pat' would find a better fit.
      tSeq[tooShort] <- NA
      geneSeq_noCDR3[cInds] <- tSeq
      # And assign the currently found trimmed sequences (could be only NA values).
    }
  }
  if (anyNA(geneSeq_noCDR3)){
    cInd <- which(is.na(geneSeq_noCDR3))[1]
    warning("Couldn't find the part of the sequence outside of CDR3 for ",
      sum(is.na(geneSeq_noCDR3)), " sequences. E.g., for ", geneSeq[cInd], " with ",
      CDR3[cInd], " from gene ", gene, ".")
  }

  # And we determine the number of AA that have been deleted from the V or J.
  # For this, we need to append the CDR3 and V/J part outside of CDR3, and then
  # remove as few as possible AAs from the full V/J sequence.
  if (gene == "V"){
    matchedSeq <- paste0(geneSeq_noCDR3, CDR3)
    fromLeft <- T
  } else {
    matchedSeq <- paste0(CDR3, geneSeq_noCDR3)
    fromLeft <- F
  }
  nAA_max <- max(nchar(matchedSeq))
  matchedSeq_mat <- stringr::str_pad(matchedSeq, width=nAA_max,
    side=ifelse(fromLeft, "right", "left"), use_width=F) %>%
    stringr::str_split(pattern="", simplify=T)
  geneSeq_mat <- stringr::str_pad(geneSeq, width=nAA_max,
    side=ifelse(fromLeft, "right", "left"), use_width=F) %>%
    stringr::str_split(pattern="", simplify=T)

  diffAA <- geneSeq_mat != matchedSeq_mat
  if (!fromLeft){
    diffAA <- diffAA[,ncol(diffAA):1]
  }
  nAA_match <- apply(diffAA, MARGIN=1, FUN=match, x=T) - 1
  # Check the first AA mismatched per tcr (either from the start if V gene
  # or from the end if J gene).
  nAA_trimmed <- nAA - nAA_match
  nAA_trimmed[is.na(nAA_match)] <- 0
  # Cases where nAA_match is.na corresponds to full equality between geneSeq and
  # matchedSeq, meaning we shouldn't trim any AA from this tcr's gene sequence.
  return(list(seq=geneSeq_noCDR3, nAA_trimmed=nAA_trimmed))
}


#' Exports .yaml files of the TCR / peptide / MHC sequences in the format needed
#' as input by Boltz.
#'
#' @param cTab The table of the TCRs to export.
#' @param outCols Names of the columns from cTab to export (these names will be
#'    used as 'id's in the yaml file).
#' @param nameCol Tells which column to use to get the IDs used to make the file
#'    names.
#' @param outFolder Folder in which to export the yaml files.
#' @param batchNr When given, this will be used to make a subolder of
#'    the outFolder telling data is from given this batchNr.
#'
#' @return nothing.
#'
#' @export
export_boltz_yaml <- function(cTab, outCols, nameCol, outFolder, batchNr=NULL){
  if (!all(c(outCols, nameCol) %in% colnames(cTab))){
    stop("Some columns from outCols or nameCol don't exist in cTab")
  }
  if (!is.null(batchNr)){
    outFolder <- paste0(outFolder, "/batch", batchNr, "/")
  }
  if (!dir.exists(outFolder)){
    dir.create(outFolder, recursive=T)
  }
  purrr::pwalk(cTab, .f=function(...){
    x <- list(...)
    oFile <- paste0(outFolder, "/", x[[nameCol]], ".yaml")
    cat("sequences:\n", file=oFile)
    for (cCol in outCols){
      cat("  - protein:",
        "\n      id: ", gsub("peptide", "Pep", cCol),
        "\n      sequence: ", x[[cCol]], "\n",
        sep="", file=oFile, append=T)
      # We can't use \t to delimit columns in yaml files, need to use multiple
      # spaces instead. Note that I'm replacing peptide" by Pep because boltz
      # had issue when the id was peptide... (probably thinking we were giving
      # the id 'pept' and then starting another 'id' with just 'e').
    }
  })
}

#' Exports .json files of the TCR / peptide / MHC sequences in the format needed
#' as input by AF3.
#'
#' @param cTab The table of the TCRs to export.
#' @param outCols Names of the columns from cTab to export (these names will be
#'    used as 'id's in the yaml file).
#' @param nameCol Tells which column to use to get the IDs used to make the file
#'    names.
#' @param outFolder Folder in which to export the yaml files.
#' @param batchNr When given, this will be used to make a subolder of
#'    the outFolder telling data is from given this batchNr.
#' @param preMSA_folder When given, we'll take the pre-aligned TCR segments to
#'    build msaPaired and msaUnpaired for AF3 inputs, speeding up computations
#'    (we'll search for files recursively in subfolders from this one to allow
#'    organizing files in subfolders). When this is given, the cTab should
#'    additionally contain the columns TRAV, n_TRAV, TRAJ, n_TRAJ, ... to get
#'    information from these; if some are missing we won't add the given pre-msa
#'    values (could e.g. only give TRAV, or TRBV+TRBJ).
#' @param templateFree indicates if we should add an empty template (so that AF3
#'     won't search for structural templates if these weren't already present),
#'     or if allow getting templates (recommended).
#' @param noPepMSA Boolean. When T, we make that no MSA is performed by AF3
#'     for the peptide chain (i.e., it'll only use the peptide sequence, but
#'     the unpairedMsa and pairedMsa and templates will be left as *empty strings*
#'     (to not make this alignment)).
#' @param modelSeeds A vector of integer seeds that we want to use if we want to
#'     run with multiple seeds.
#'
#' @return nothing.
#'
#' @export

.remap_tcr_templates <- function(templates, seg_len, n_trim, full_len, isV){
  if (is.null(templates) || length(templates) == 0) return(NULL)

  out <- list()

  for (tpl in templates){
    qidx <- unlist(tpl$queryIndices)
    tidx <- unlist(tpl$templateIndices)

    if (isV){
      keep <- qidx <= (seg_len - n_trim)
      qidx_new <- qidx[keep]
      tidx_new <- tidx[keep]
    } else {
      keep <- qidx > n_trim
      qidx_new <- qidx[keep] + full_len - seg_len
      tidx_new <- tidx[keep]
    }

    if (length(qidx_new) == 0) next

    tpl$queryIndices <- as.list(qidx_new)
    tpl$templateIndices <- as.list(tidx_new)
    out[[length(out) + 1]] <- tpl
  }

  out
}

export_AF3_json <- function(cTab, outCols, nameCol, outFolder, batchNr=NULL,
  preMSA_folder=NULL, templateFree=F, noPepMSA=F, noTCRMSA=F, modelSeeds=1){
  if (("MHC_b" %in% outCols) && !"MHC_b" %in% colnames(cTab)){
    warning("MHC_b column absent from the cTab to export for AF3 json.",
      "\n Please make sure this is expected.")
    outCols <- setdiff(outCols, "MHC_b")
  }
  if (!all(c(outCols, nameCol) %in% colnames(cTab))){
    stop("Some columns from outCols or nameCol don't exist in cTab")
  }
  if (!is.null(batchNr)){
    outFolder <- paste0(outFolder, "/batch", batchNr, "/")
  }
  if (!dir.exists(outFolder)){
    dir.create(outFolder, recursive=T)
  }
  preMSA <- list()
  if (!is.null(preMSA_folder)){
    avail_MSA <- list.files(preMSA_folder, pattern=".json$", recursive=T)
    # Reading the pre-aligned MSA files needed for the current TCRs.
    for (cCol in outCols){
      if (cCol %in% c("TCRa", "TCRb")){
        cChain <- toupper(gsub("TCR", "", cCol))
        subCols <- paste0("TR", cChain, c("V", "J"))
        if (!all(c(subCols, paste0("n_", subCols)) %in% colnames(cTab))){
          warning("Some subCols columns seem missing from cTab to be able to ",
            "add all the preMSA information")
        }
      } else if (grepl("^MHC_", cCol)){
        subCols <- gsub("MHC_", "MHC_allele_", cCol)
        if (!subCols %in% colnames(cTab)){
          warning(subCols, " column is missing from cTab, can't add related ",
            "preMSA values.")
        }
      } else {
        subCols <- cCol
      }
      for (cSubCol in subCols){
        x <- unique(cTab[[cSubCol]])
        for (cx_o in x){
          cx <- tolower(gsub("\\.", "_", make.names(gsub("\\*", "__", cx_o))))
          # Transform to same naming convention than used for the files, replacing
          # the special characters, etc. (and AF3 outputs file names in lower
          # letters only).
          if (grepl("__.*$", cx)){
            warning("For the moment we don't consider the alleles, just the ",
              "main gene sequences for prealigned MSAs (would need to run ",
              "msa separately for each allele to get corresponding files if ",
              "we want these as well).")
            cx <- gsub("__.*$", "", cx)
          }
          if ((cCol == "peptide") && noPepMSA){
            cfile <- grep(paste0("/emptyPeptide_no_msa.json"), avail_MSA, value=T)
            # Special case, not reading preMSA from current peptide but one
            # telling to not make any msa (i.e., with field ' "unpairedMsa": "" '
            # for example).
          } else {
            cfile <- grep(paste0("/", cx, "_data.json"), avail_MSA, value=T)
          }
          if (length(cfile) > 1){
            stop("There are multiple preMSA files matching ", cx, " data.")
          } else if (length(cfile) == 0){
            warning("There wasn't any preMSA file for ", cx)
          } else {
            y <- jsonlite::read_json(paste0(preMSA_folder, "/", cfile))
            y <- y$sequences[[1]]$protein
            if (cCol %in% c("TCRa", "TCRb")){
              y$unpairedMsa <- gsub("^>query\\n[[:alpha:]]*\\n", "", y$unpairedMsa)
              y$pairedMsa <- gsub("^>query\\n[[:alpha:]]*\\n", "", y$pairedMsa)
              # Remove the first "line" from these MSA of TCRs as we'll need
              # to replace these by the full-length query sequence.
            } else if ((cCol == "peptide") && noPepMSA){
              y$sequence <- cx_o
              # To have the correct sequence as we used a 'universal' file
              # not specific to the current peptide.
            }
            preMSA[[cCol]][[cx_o]] <- y
          }
        }
      }
    }
  }
  purrr::pwalk(cTab, .f=function(...){
    x <- list(...)
    out <- list(name=x[[nameCol]], modelSeeds=as.list(modelSeeds), dialect="alphafold3",
      version=3)
    for (cCol in outCols){
      cID <- gsub("[^[:upper:]]", "", toupper(cCol))
      # Need to have upper-case letters only for the id, otherwise AF3 doesn't
      # work. So we first make sure small letters are replaced by capital ones
      # and then also remove symbols (like '_' there were possibly present).
      pStruct <- list(list(protein=list(id=cID, sequence=x[[cCol]])))
      if (!is.null(preMSA[[cCol]])){
        # Will use corresponding alignments when available.
        if (cCol %in% c("TCRa", "TCRb")){
          cChain <- toupper(gsub("TCR", "", cCol))
          subCols <- paste0("TR", cChain, c("V", "J"))
          tcrSeq <- x[[cCol]]
          nAA <- nchar(tcrSeq)

          for (cSubCol in subCols){
            cMSA <- preMSA[[cCol]][[x[[cSubCol]]]]
            if (!is.null(cMSA) && !is.null(x[[paste0("n_", cSubCol)]])){
              n_trim <- x[[paste0("n_", cSubCol)]]
              seg_len <- nchar(cMSA$sequence)
              n_add <- nAA + n_trim - seg_len

              if (grepl("TR.V", cSubCol)){
                pattern <- paste0("([[:lower:]]*[-[:upper:]]){", n_trim, "}",
                  "(\\n(>|$))")
                repl <- paste0(stringr::str_dup("-", times=n_add), "\\2")
                isV <- TRUE
              } else if (grepl("TR.J", cSubCol)){
                pattern <- paste0("(\\n)(?!(>|$))([-[:upper:]][[:lower:]]*){",
                  n_trim, "}")
                repl <- paste0("\\1", stringr::str_dup("-", times=n_add))
                isV <- FALSE
              } else {
                stop("cSubCol doesn't seem to match TR.V/J: ", cSubCol)
              }
              if (noTCRMSA) {
                pStruct[[1]]$protein$unpairedMsa <- ""
                pStruct[[1]]$protein$pairedMsa   <- ""
              } else {
                if (!is.null(cMSA$unpairedMsa)){
                  cMSA$unpairedMsa <- gsub(pattern, repl, x=cMSA$unpairedMsa, perl=TRUE)
                  pStruct[[1]]$protein$unpairedMsa <-
                    paste0(pStruct[[1]]$protein$unpairedMsa, cMSA$unpairedMsa)
                }

                if (!is.null(cMSA$pairedMsa)){
                  cMSA$pairedMsa <- gsub(pattern, repl, x=cMSA$pairedMsa, perl=TRUE)
                  pStruct[[1]]$protein$pairedMsa <-
                    paste0(pStruct[[1]]$protein$pairedMsa, cMSA$pairedMsa)
                }
              }

              ## NEW: carry over templates too
              cTemplates <- .remap_tcr_templates(
                templates = cMSA$templates,
                seg_len = seg_len,
                n_trim = n_trim,
                full_len = nAA,
                isV = isV
              )

              if (!is.null(cTemplates)){
                pStruct[[1]]$protein$templates <-
                  c(pStruct[[1]]$protein$templates, cTemplates)
              }
            }
          }

          if (!noTCRMSA &&
            (!is.null(pStruct[[1]]$protein$unpairedMsa) ||
            !is.null(pStruct[[1]]$protein$pairedMsa))) {
          qu_str <- paste0(">query\n", tcrSeq, "\n")
          pStruct[[1]]$protein$unpairedMsa <- paste0(
            qu_str, pStruct[[1]]$protein$unpairedMsa
          )
          pStruct[[1]]$protein$pairedMsa <- paste0(
            qu_str, pStruct[[1]]$protein$pairedMsa
          )
        }
        } else {
          cCol_i <- cCol
          if (grepl("^MHC_", cCol)){
            cCol_i <- gsub("MHC_", "MHC_allele_", cCol)
          }
          if (!is.null(preMSA[[cCol]][[x[[cCol_i]]]])){
            pStruct[[1]]$protein <- preMSA[[cCol]][[x[[cCol_i]]]]
            # Can directly use all pre-MSA, templates, sequence, ... for
            # the peptide and MHC matches.
          }
        }
      }
      if (templateFree){
        if (is.null(pStruct[[1]]$protein$templates)){
          pStruct[[1]]$protein$templates <- list()
          # Put to an empty templates so that AF3 doesn't search for it if
          # it isn't already present from preMSA files (when templateFree=T).
        }
      }
      if (!purrr::is_empty(pStruct[[1]]$protein$templates)){
        # Doing a little trick so that the queryIndices and templateIndices (if
        # present), will be on a single line: otherwise write_json with pretty
        # was putting each index on a separate row.
        for (i in seq_along(pStruct[[1]]$protein$templates)){
          pStruct[[1]]$protein$templates[[i]]$queryIndices <-
            unlist(pStruct[[1]]$protein$templates[[i]]$queryIndices)
          pStruct[[1]]$protein$templates[[i]]$templateIndices <-
            unlist(pStruct[[1]]$protein$templates[[i]]$templateIndices)
        }
      }

      out$sequences <- c(out$sequences, pStruct)
    }
    oFile <- paste0(outFolder, "/", x[[nameCol]], ".json")
    jsonlite::write_json(out, path=oFile, pretty=T, auto_unbox=T)
  })
}


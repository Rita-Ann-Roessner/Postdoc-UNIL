# Script used to build the input files used by AF3 or Boltz, based on some table
# of TCR sequences.

# Input arguments 1 to 2 are compulsory and 3 to 8 are optional:
#   1. In general this is the full path to a tsv/csv file containing
#      a table of TCRs. If it has an 'id' column, this will be used to name the
#      output files (if this is a numeric column, we'll add 'tcr' in front of
#      the numbers).
#      MODIFICATION: Also a single TCR chain can be supplied.
#     This first argument, can however also be "preMSA" to build the alignment
#      files for additional MHC alleles; it will then change the input options
#      meaning (see below for explanations after the standard ones).
#   2. Output folder where to save the exported files.
#   3. exportType, indicating for which predictor we want to output the TCRs
#     (AF3 (default), Boltz, or AF3_Boltz to have both (can be lower or upper
#     letters)).
#   4. Sequence of the peptide.
#   5. Name of the first MHC allele (alpha chain for MHC-II).
#   6. Name of the 2nd "MHC allele" (B2M for MHC-I or beta chain for MHC-II).
#     It can also be 'none' if we don't want to include the B2M in the
#     structural predictions.
#   7. Species (HomoSapiens or MusMusculus).
#   8. An index of the seed to use (default to 1) - it can be multiple values
#     separated by coma to use multiple seeds (e.g. 1,17,32891).

# The inputs 4-7 can be omitted if they are given in the input tcr tab file
# (columns need to be named "peptide" (or "epitope"), "MHC_allele_a" (or "MHC"),
# "MHC_allele_b" and "species" respectively). If some columns are to be taken
# from the tcr tab and some given here (e.g. we want to use seeds), we can
# replace their values by "input".
# If MHC_allele_b is not given and MHC is of class I, we'll use B2M (only for
# human for the moment, would need to update if useful for mice as well).
# If species is not given it'll be determined based on the MHC_allele_a.
# Note that all TCRs need to be from the same species, if not, split the input
# in two files.

# When input argument 1 is equal to preMSA, the other input arguments have
# the following meaning:
#   2. Output folder where to save the exported files.
#   3. and following: names of the alleles for which we want to build input
#      files to build their alignment files with AF3.

# Preparing environment and getting input values ---------------------------------
if (Sys.info()["sysname"] == "Linux"){
  .libPaths("/work/FAC/FBM/LLB/dgfeller/epitope_pred/R_library/4.2")
  # To be able to load installed libraries on the cluster.
  source(paste0("/work/FAC/FBM/LLB/dgfeller/epitope_pred/MixTCR_internal/",
    "Struct_preds/Struct_preds_fun.R"))
  # File defining some functions used here.
}

library(tidyverse)
# Load this set of packages that are used in multiple places, allowing also
# the use of ' %>% '.

# preMSA_folder <- paste0(dataRootFolder, "/TCR/Struct_preds/preMSA/",
#   "AF3_models/")
preMSA_folder <- paste0("/work/FAC/FBM/LLB/dgfeller/epitope_pred/jracle/preMSA/",
  "AF3_models/")
# Will add the 'species' subFolder path below after having determined it.

allelesTabFile <- paste0(dirname(preMSA_folder), "/Alleles_tab.txt")
MHCallelesTab <- read_tsv(allelesTabFile, show_col_types=F)
# Importing the MHCallelesTab containing the AA sequences of the MHCs (and B2M).


batchSize <- 800
# To output sequences in batches to run boltz/AF3 separately on each batch. For
# AF3, it seems that 800 is a good compromise between the waiting time in the
# queue and the running time.

noPepMSA <- T
# Set to T, so that we won't use MSAs for the peptide part (for AF3 predictions),
# and make that AF3 doesn't search for this MSA (in general it finds few
# sequences against which to align the peptide and it shouldn't thus have a big
# impact on the results (this avoid the need of pre-aligning each peptide that
# we want to study)).

if (!interactive()){
  args <- commandArgs(trailingOnly=T)
} else {
  if (!file.exists("Struct_preds_fun.R")){
    stop("When working in interactive mode, current working directory ",
      "should contain the file 'Struct_preds_fun.R', which isn't the case.")
  }
  source("Struct_preds_fun.R")

  # args <- c("test_functions/Struct_preds/TCRtab_example.txt",
  #   "test_functions/Struct_preds/outFolder/")
  tx <- ""
  i <- 1
  args <- c()
  while (T) {
    tx <- readline(paste("Please enter input argument", i, "(or 'q' to finish): "))
    if (tx != "q"){
      args <- c(args, tx)
      i <- i+1
    } else {
      break
    }
  }
}

jrUtils::cat_withTime("Starting preparing inputs for structural predictors.",
  "Input arguments were:\n  -", paste(args, collapse="\n  - "), "\n")

tcrFile <- args[1]
outPath <- args[2]

if (tcrFile != "preMSA"){

# Getting input values in the case of a standard tcrFile in the input ---------------
peptide <- MHC_allele_a <- MHC_allele_b <- species <- NULL
if (length(args) >= 3){
  exportType <- args[3]
} else {
  exportType <- "AF3"
}
if ((length(args) >= 4) && ("input" != args[4])){
  peptide <- args[4]
}
if ((length(args) >= 5) && ("input" != args[5])){
  MHC_allele_a <- args[5]
}
if ((length(args) >= 6) && ("input" != args[6])){
  MHC_allele_b <- args[6]
}
if ((length(args) >= 7) && ("input" != args[7])){
  species <- args[7]
}
if (length(args) >= 8){
  modelSeeds <- as.numeric(str_split_1(args[8], pattern=","))
} else {
  modelSeeds <- 1
}


# Reading the input tab and transforming/getting some variables based on this -----------
exportBoltz <- grepl("boltz", tolower(exportType))
exportAF3 <- grepl("af3", tolower(exportType))
# To export the sequences in the format needed by Boltz and/or AF3.

if (grepl(".csv", tcrFile)){
  tcrTab <- read_csv(tcrFile)
} else {
  tcrTab <- read_tsv(tcrFile)
}

# Possibly add the id column (either from first column if it didn't have a name
# or creating this id column if absent). This will be used to name the various
# files created.
if (colnames(tcrTab)[1] == "...1"){
  colnames(tcrTab)[1] <- "id"
}
if (! "id" %in% colnames(tcrTab)){
  tcrTab <- add_column(tcrTab, .before=1, id=seq_len(nrow(tcrTab)))
}
if (is.numeric(tcrTab$id)){
  # Add tcr in front to not have filenames starting with numbers (and make
  # it'll be numbers with 4 digits (or more))
  tcrTab$id <- paste0("tcr", sprintf("%04.0f", tcrTab$id))
}

# Checking if the MHC/peptide/species are given or obtaining from the tcrTab
# or from inferred values.
if (is.null(MHC_allele_a)){
  cCol <- grep("^(MHC|MHC_allele_a)$", colnames(tcrTab))
  if (length(cCol) > 1){
    warning("Columns MHC and MHC_allele_a are both present or appear multiple ",
      "times in tcrTab, will use first one.")
    cCol <- cCol[1]
  }
  if (length(cCol) == 1){
    MHC_allele_a <- tcrTab[[cCol]]
  } else {
    stop("MHC_allele_a wasn't given in input and isn't found in tcrTab ",
      "(or the columns MHC and MHC_allele_a are both present)")
  }
}
MHC_allele_a <- paste0(ifelse(grepl("^(A|B|C)", MHC_allele_a), "HLA_", ""),
  MHC_allele_a)
# Allele needs to be of the form HLA_A0201, not just A0201 (this simple check
# will work in standard cases, but we don't check for mouse or other types).

if (is.null(MHC_allele_b)){
  if ("MHC_allele_b" %in% colnames(tcrTab)){
    MHC_allele_b <- tcrTab[["MHC_allele_b"]]
  } else {
    if (all(grepl("^HLA_(A|B|C)", MHC_allele_a))){
      MHC_allele_b <- "B2M"
    } else {
      stop("MHC_allele_b wasn't given in input, isn't found in tcrTab ",
        "and it seems that not all MHC_allele_a are from human HLA-I, so ",
        "we can't infer this MHC_allele_b value.")
    }
  }
}

if (is.null(peptide)){
  cCol <- grep("^(peptide|epitope)$", colnames(tcrTab))
  if (length(cCol) > 1){
    warning("Columns peptide and epitope are both present or appear multiple ",
      "times in tcrTab, will use first one.")
    cCol <- cCol[1]
  }
  if (length(cCol) == 1){
    peptide <- tcrTab[[cCol]]
  } else {
    stop("peptide wasn't given in input and isn't found in tcrTab ",
      "(or the columns peptide and epitope are both present)")
  }
}

if (is.null(species)){
  if ("species" %in% colnames(tcrTab)){
    species <- unique(tcrTab[["species"]])
    if (length(species) != 1){
      stop("Species was taken from the input tcrTab. We need to have only ",
        "one species, which wasn't the case.")
    }
  } else {
    if (all(grepl("^HLA", MHC_allele_a))){
      species <- "HomoSapiens"
    } else if (all(grepl("^H2", MHC_allele_a))){
      species <- "MusMusculus"
    } else {
      stop("species wasn't given in input, isn't found in tcrTab ",
        "and it seems that not all MHC_a are from the same species, which ",
        "is needed.")
    }
  }
}
preMSA_folder <- paste0(preMSA_folder, species, "/")

outCols <- c("TCRa", "TCRb", "peptide", "MHC_a", "MHC_b")
# Columns that will be outputted from the tcrTab.

tcrTab$MHC_allele_a <- MHC_allele_a
tcrTab$MHC_allele_b <- MHC_allele_b
tcrTab$MHC_a <- MHCallelesTab$seq[match(MHC_allele_a, MHCallelesTab$Allele)]
if (all(tolower(MHC_allele_b) == "none")){
  outCols <- setdiff(outCols, "MHC_b") # Not outputting "MHC_b"
} else {
  tcrTab$MHC_b <- MHCallelesTab$seq[match(MHC_allele_b, MHCallelesTab$Allele)]
}
tcrTab$peptide <- peptide

# Transforming the TCR sequences with V/J names to full AA sequence --------

# Check if we are running in single chain mode.
alpha_cols <- c("TRAV", "TRAJ", "cdr3_TRA")
beta_cols  <- c("TRBV", "TRBJ", "cdr3_TRB")

# Check alpha chain validity
alpha_valid <- 
  all(alpha_cols %in% colnames(tcrTab)) &&
  !all(is.na(tcrTab[, alpha_cols, drop = FALSE]))

if (!alpha_valid){
outCols <- setdiff(outCols, "TCRa")} # Not outputting "TCRa"

# Check beta chain validity
beta_valid <- 
  all(beta_cols %in% colnames(tcrTab)) &&
  !all(is.na(tcrTab[, beta_cols, drop = FALSE]))

if (!beta_valid){
outCols <- setdiff(outCols, "TCRb")} # Not outputting "TCRb"

if (sum(c(alpha_valid, beta_valid)) == 1){
  message('Running in single chain mode.')}

# Initialize logical vector for bad sequences
badIds <- rep(FALSE, nrow(tcrTab))

if (alpha_valid) {
  TRAV <- get_TRseq("V", "A", tcrTab, species=species)
  TRAJ <- get_TRseq("J", "A", tcrTab, species=species)

  tcrTab <- bind_cols(tcrTab, list(
  TCRa=paste0(TRAV$seq, tcrTab$cdr3_TRA, TRAJ$seq),
  n_TRAV=TRAV$nAA_trimmed, n_TRAJ=TRAJ$nAA_trimmed))

  bad_alpha <- (!grepl(paste0("^[", paste(c("-", MHC2R::aa_letters), collapse=""), "]*$"), tcrTab$TCRa)) | is.na(tcrTab$cdr3_TRA) | is.na(TRAV$seq) | is.na(TRAJ$seq)
  badIds <- badIds | bad_alpha
}

if (beta_valid) {
  TRBV <- get_TRseq("V", "B", tcrTab, species=species)
  TRBJ <- get_TRseq("J", "B", tcrTab, species=species)
  tcrTab <- bind_cols(tcrTab, list(
  TCRb=paste0(TRBV$seq, tcrTab$cdr3_TRB, TRBJ$seq),
  n_TRBV=TRBV$nAA_trimmed, n_TRBJ=TRBJ$nAA_trimmed))

  bad_beta <- (!grepl(paste0("^[", paste(c("-", MHC2R::aa_letters), collapse=""), "]*$"), tcrTab$TCRb)) | is.na(tcrTab$cdr3_TRB) | is.na(TRBV$seq) | is.na(TRBJ$seq)
  badIds <- badIds | bad_beta
}

# Warning and remove bad sequences
badIds <- which(badIds)
if (length(badIds) > 0) {
  warning(length(badIds), " TCRs contain bad AAs or missing gene/CDR3 info (e.g. ",
    paste(tcrTab$id[badIds[1:min(3, length(badIds))]], collapse=", "), "). Removing these from the exported TCRs.")
  tcrTab <- tcrTab[-badIds,]
}

tcrTab <- distinct(tcrTab, pick(all_of(outCols)), .keep_all=T)

# Exporting the json or yaml files ----------------------------------------
nBatches <- ceiling(nrow(tcrTab)/batchSize)
walk(seq_len(nBatches), .f=function(i){
  cInds <- ((i-1)*batchSize+1) : min(i*batchSize, nrow(tcrTab))
  tcrTab <- tcrTab[cInds,]
  if (exportBoltz){
    export_boltz_yaml(tcrTab, nameCol="id", outFolder=paste0(outPath, "/Boltz_input"),
      batchNr=i, outCols=outCols)
  }
  if (exportAF3){
    export_AF3_json(tcrTab, nameCol="id", outFolder=paste0(outPath, "/AF3_input"),
      batchNr=i, outCols=outCols, preMSA_folder=preMSA_folder, noPepMSA=noPepMSA,
      modelSeeds=modelSeeds)
  }
})

# Exporting the updated tcrTab, including the id column and some other columns -------
# that were possibly added, and without the duplicates or sequences that had
# issues.
tcrTab <- select(tcrTab, -any_of(c("MHC_a", "MHC_b", "TCRa",
  "TCRb", "n_TRAV", "n_TRAJ", "n_TRBV", "n_TRBJ")))

# Don't keep these columns with long sequences of AAs because they aren't really
# used further and could be reobtained from here.
oFile_root <- gsub("\\.(txt|tsv|csv)$", "", basename(tcrFile))
oFile <- jrUtils::f_new_file(paste0(outPath, "/", oFile_root, "_input.txt"))
write_tsv(tcrTab, file=oFile)


} else {
# Instead of standard tcrFile, we want to build files for MHCs pre-alignment ---------------
  if (length(args) < 3){
    stop("When using the 'preMSA' option, there should be at lease 3 input ",
      "arguments.")
  }
  MHC_s <- args[-(1:2)]
  MHCtab <- filter(MHCallelesTab, Allele %in% c(MHC_s, paste0("HLA_", MHC_s))) %>%
    rename(MHC="seq")
  if (nrow(MHCtab) != length(MHC_s)){
    stop("It seems that at least some MHCs asked in input for preMSA are absent ",
      "from the allelesTab from ", allelesTabFile)
  }

  export_AF3_json(MHCtab, outCols="MHC", nameCol="Allele",
    outFolder=outPath)


}


# End, outputting some text -----------------------------------------------
jrUtils::cat_withTime("Finished preparing input files for structual predictors.\n",
  "Following warning messages appeared during the analysis\n  - ",
  paste0(names(jrUtils::warnings_u()), collapse="\n  - "), "\n", sep="")


# EOF ---------------------------------------------------------------------



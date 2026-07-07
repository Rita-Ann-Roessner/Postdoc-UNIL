library(MixTCRviz)
library(pdftools)

peptides <- c("ELAGIGILTV", "GILGFVFTL", "GLCTLVAML", "LLWNGPMAV")

for (peptide in peptides) {

  binders_file <- file.path(peptide, sprintf("A0201_%s.csv", peptide))
  out_dir      <- file.path(peptide, "motif_plain_tmp")

  if (!file.exists(binders_file)) {
    message(sprintf("Skipping %s — file not found: %s", peptide, binders_file))
    next
  }

  binders <- read.csv(binders_file, stringsAsFactors = FALSE)
  message(sprintf("Building reference motif for %s (%d TCRs)", peptide, nrow(binders)))

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  MixTCRviz(
    input1     = binders,
    output.path = out_dir,
    plot       = TRUE,
    renormVJ   = TRUE,
    set.title  = sprintf("A0201_%s", peptide),
    verbose    = 0
  )

  # MixTCRviz writes a PDF; convert the first page to PNG
  pdf_file <- list.files(out_dir, pattern = "\\.pdf$", full.names = TRUE)[1]
  png_out  <- file.path(peptide, "motif_plain.png")

  if (!is.na(pdf_file) && file.exists(pdf_file)) {
    # Use pdftools if available, otherwise fall back to bitmap()
    pdftools::pdf_convert(pdf_file, format = "png", pages = 1,
                          filenames = png_out, dpi = 150)
    if (file.exists(png_out)) {
      message(sprintf("  Saved: %s", png_out))
    }
  } else {
    message(sprintf("  No PDF found in %s for %s", out_dir, peptide))
  }
}

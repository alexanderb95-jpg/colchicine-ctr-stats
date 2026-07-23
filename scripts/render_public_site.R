#!/usr/bin/env Rscript
# Render slim GitHub Pages site (methods, results, figures only).
root <- if (requireNamespace("here", quietly = TRUE)) {
  here::here()
} else {
  normalizePath(".", winslash = "/")
}
setwd(root)

if (!file.exists("CTR_Submission_Results_Protocol2101175.Rmd")) {
  stop("Run from repository root (CTR_Submission_Results_Protocol2101175.Rmd not found).")
}

rmarkdown::render(
  "CTR_Submission_Results_Protocol2101175.Rmd",
  output_file = "index.html",
  params = list(public_site = TRUE),
  quiet = FALSE
)

message("Wrote index.html (public_site = TRUE).")

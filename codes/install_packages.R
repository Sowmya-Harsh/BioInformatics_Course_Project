# HOW TO RUN ON IRIS HPC:
# module load env/release/2023b
# module load lang/R-bundle-CRAN/2024.06-foss-2023b
# Rscript install_packages.R

rm(list = ls())
.libPaths(c("~/Rlibs", .libPaths()))

# CRAN packages
pkgs_cran <- c(
  "Seurat", "tidyverse", "dplyr", "ggplot2",
  "viridis", "viridisLite",
  "pheatmap", "RColorBrewer",
  "ggrepel", "patchwork",
  "openxlsx", "future", "scales",
  "stringr", "tidyr", "tibble"
)

for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org", lib = "~/Rlibs")
  }
}

# Bioconductor (for org.Hs.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = "~/Rlibs")
}
BiocManager::install(
  c("org.Hs.eg.db", "enrichplot"),  #clusterprofiler not available
  lib = "~/Rlibs"
)

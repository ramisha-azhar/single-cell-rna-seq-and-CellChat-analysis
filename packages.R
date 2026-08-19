if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "BiocNeighbors"))

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("sqjin/CellChat")
devtools::install_github('immunogenomics/presto')


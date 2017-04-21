install.packages("devtools")

source("https://bioconductor.org/biocLite.R")  # Install Bioconductor
biocLite(c("Biostrings","iterators","ade4"))                                     # Install required packages

devtools::install_github("crarlus/paprbag") # install paprbag

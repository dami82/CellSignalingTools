\name{CellSignalingTools-package}
\alias{CellSignalingTools-package}
\alias{CellSignalingTools}
\docType{package}
\title{
Analyze Fluoresence Microscopy Data and Scientific Images for Molecular Biology Applications
}
\description{
Analyze Fluoresence Microscopy data and calculate nuclear to cytoplasmic signal ratios based on signal profiles or stacked immunofluorescence TIFF images. 
}
\author{
Damiano Fantini <damiano.fantini@gmail.com>
}
\references{
http://www.biotechworld.it/bioinf/2016/03/09/analyzing-fluoresence-microscopy-data-with-r/
}
\examples{
####-------- Example 01: IF data analysis ----------------
####
data("leio_cells_dapi")
data("leio_cells_green")
dapi_mask <- get_dapi_range(leio_cells_dapi, 0.05)
d_msk <- dapi_mask[800:1100, 800:1100]
g_img <- leio_cells_green[800:1100, 800:1100]
nc_signal <- apply_mask(g_img, d_msk, bckgrnd = 0.05)
nc_ratio_tiff(nc_signal)
plot_signal(nc_signal)
####
####-------- Example 02: signal profile analysis ---------
####
data("leio_basal_profile")
data("leio_pq100_profile")
items <- list(leio_basal_profile, leio_pq100_profile)
conditions <- c('CTR','PQ100')
results <- lapply(1:2, (function(n){
  nc_ratio_profile(items[[n]], 
  "baseline",
  "mean", 
  conditions[n], 
  "PTEN", 
  0.065)
})) 
}

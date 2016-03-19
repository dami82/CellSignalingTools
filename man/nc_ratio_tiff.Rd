\name{nc_ratio_tiff}
\alias{nc_ratio_tiff}
\title{
Determine Nucleus to Cytoplasm Signal Ratio
}
\description{
Determine the ratio of Nuclear to non-nuclear Signal intensity from processed stacked TIFF images. It also plots a graphic representation of nuclear and non-nuclear signal. 
}
\usage{
nc_ratio_tiff(signal_data, signal_cutoff = 0)
}
\arguments{
  \item{signal_data}{
is a list containing 3 numeric matrices of identical dimensions. Typically, this is the result of a apply_mask() call.
}
  \item{signal_cutoff}{
a number in the range 0 to 1.}
}
\value{
It returns a number corresponding to the ratio between nuclear and non-nuclear signal.
}
\references{
http://www.biotechworld.it/bioinf/2016/03/09/analyzing-fluoresence-microscopy-data-with-r/
}
\author{
Damiano Fantini
}
\note{
It is advisable not to set a signal cutoff.
}
\examples{
data("leio_cells_dapi")
data("leio_cells_green")
dapi_mask <- get_dapi_range(leio_cells_dapi, 0.05)
d_msk <- dapi_mask[500:1100, 700:1100]
g_img <- leio_cells_green[500:1100, 700:1100] 
nc_signal <- apply_mask(g_img, d_msk, bckgrnd = 0.05)
nc_ratio_tiff(nc_signal)
}

\name{plot_signal}
\alias{plot_signal}
\title{
Generates an Graphic Represntation of the Nuclear and Extranuclear Signal
}
\description{
Generates a four-panel graphics that summarizes the results of the nuclear to cytoplasm signal distribution analysis. Panels include an image of the input signal, an image of the nuclear signal, an image of the non-nuclear signale and a faux-color image showing both nuclear and non-nuclear signal in the same image.
}
\usage{
plot_signal(signal_map)
}
\arguments{
  \item{signal_map}{
  is a list containing 3 numeric matrices of identical dimensions. Typically, this is the result of a apply_mask() call.
}
}
\references{
http://www.biotechworld.it/bioinf/2016/03/09/analyzing-fluoresence-microscopy-data-with-r/
}
\author{
Damiano Fantini
}
\examples{
data("leio_cells_dapi")
data("leio_cells_green")
dapi_mask <- get_dapi_range(leio_cells_dapi, 0.05)

d_mask <- dapi_mask[500:1100, 700:1100]
leio_img <- leio_cells_green[500:1100, 700:1100]
nc_signal <- apply_mask(leio_img, d_mask, bckgrnd = 0.05)
plot_signal(nc_signal)
}

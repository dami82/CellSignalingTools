\name{get_dapi_range}
\alias{get_dapi_range}
\title{
Define a DAPI Mask
}
\description{
Use a numeric matrix from a DAPI TIFF image to create a DAPI Mask. It identifies high signal DAPI areas (corresponding to nuclei) and converts the image to a boolean matrix. }
\usage{
get_dapi_range(dapi_measures, cutoff = "auto")
}
\arguments{
  \item{dapi_measures}{
is a numeric matrix corresponding to an imported 8-bit TIFF image.
}
  \item{cutoff}{
is a numeric value in the range 0 to 1 and is used for the backgroung correction.
}
}
\value{
It returns a boolean matrix.
}
\references{
http://www.biotechworld.it/bioinf/2016/03/09/analyzing-fluoresence-microscopy-data-with-r/
}
\author{
Damiano Fantini <damiano.fantini@gmail.com>
}
\examples{
data("leio_cells_dapi")
dapi_mask <- get_dapi_range(leio_cells_dapi, 0.05)
my_cols <- c("gray75", "gray45")
image(dapi_mask[400:1200,500:1200], col = my_cols)
}

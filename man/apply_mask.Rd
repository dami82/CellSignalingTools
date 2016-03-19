\name{apply_mask}
\alias{apply_mask}
\title{
Apply a DAPI Mask over a Stacked Green Channel Image
}
\description{
Apply a DAPI mask over a specific antibody signal matrix (retrieved from an immunofluorescence TIFF image). Split the signal corresponding to a specific antibody (usually the green or the red channel) into nuclear and extra-nuclear areas based on DAPI Mask. Images in the Green/Red Channel and in the DAPI channel have to be perfectly stacked. 
}
\usage{
apply_mask(green_image, dapi_mask, bckgrnd = 0)
}
\arguments{
  \item{green_image}{
is a numeric matrix corresponding to signal intensity from an immunofluorescence TIFF image
}
  \item{dapi_mask}{
is a boolean matrix having the same dimensions as green_image. Typically, this is the result of a get_dapi_range() call.
}
  \item{bckgrnd}{
is a number in the range 0 to 0.9 and is used for background subtraction of the scaled (0 to 1) signals of the green_image
}
}
\value{
Returns a list of three numeric matrices having identical dimensions.
  \item{nucl }{Signal intensities corresponding to nuclear areas of the green_image }
  \item{cyto }{Signal intensities corresponding to extra-nuclear areas of the green_image }
  \item{input }{Scaled and background-corrected input matrix (green_image) }
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
leio_img <- leio_cells_green[500:1100, 700:1100]
d_mask <- dapi_mask[500:1100, 700:1100]
nc_signal <- apply_mask(leio_img, d_mask, bckgrnd = 0.05)
plot_signal(nc_signal)
}
\name{remove_dust}
\alias{remove_dust}
\title{
Clean and Adjust Signal of a DAPI Mask 
}
\description{
Clean and Adjust Signal of a DAPI Mask. Voids in large signal areas are filled while small signal areas (aspecific signal) are removed. 
}
\usage{
remove_dust(signal_map, nuclei_dims)
}
\arguments{
  \item{signal_map}{
is a boolean matrix. Typically, this is the result of a get_dapi_range() call. 
}
  \item{nuclei_dims}{
is an integer vector containing two elements. Typically, this is the result of a get_nucl_dims() call.
}
}
\value{
The function returns a boolean matrix
}
\references{
http://www.biotechworld.it/bioinf/2016/03/09/analyzing-fluoresence-microscopy-data-with-r/
}
\author{
Damiano Fantini
}
\examples{
data("leio_cells_dapi")
dapi_mask <- get_dapi_range(leio_cells_dapi)
n_dims <- get_nucl_dims(dapi_mask)
clean_mask <- remove_dust(dapi_mask[300:800,400:800], n_dims)
image(clean_mask, col = c("gray75", "white"))
}

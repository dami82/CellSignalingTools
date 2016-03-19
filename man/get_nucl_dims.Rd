\name{get_nucl_dims}
\alias{get_nucl_dims}

\title{
Determine Maximum Width and Height of Nuclei in a DAPI Mask  
}
\description{
Determine the maximum width and the maximum height (in pixels) of DAPI areas identified as nuclei. 
}
\usage{
get_nucl_dims(dapi_range)
}
\arguments{
  \item{dapi_range}{
is a boolean matrix. Typically, this is the result of a get_dapi_range() call. 
}
}
\value{
It returns an integer vector containing two elements
}
\references{
http://www.biotechworld.it/bioinf/2016/03/09/analyzing-fluoresence-microscopy-data-with-r/
}
\author{
Damiano Fantini
}
\examples{
data("leio_cells_dapi")
dapi_mask <- get_dapi_range(leio_cells_dapi, 0.05)
get_nucl_dims(dapi_mask)
}

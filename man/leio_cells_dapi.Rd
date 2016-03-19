\name{leio_cells_dapi}
\alias{leio_cells_dapi}
\docType{data}
\title{
Leiomyoma Cells Immunofluorescence Signal Intensities in the DAPI Channel
}
\description{
Numeric matrix of signal intensity corresponding to the DAPI channel of an immunofluorescence experiment performed on Leiomyoma cells. The dataset was retrieved from a stack of two images that also includes the leio_cells_green dataset.
}
\usage{data("leio_cells_dapi")
}
\details{
Leiomyoma cells were analyzed by immunofluorescence. Fixed cells were stained with an anti-PTEN-Alexa488 antibody and counterstained with DAPI. 8-bit grayscale images in the green and the DAPI channels were acquired. Stacked TIFF images corresponding to both channels were imported via the "tiff" package and saved as numeric matrices.
}
\source{
Vidimar V et al, 2016? (in preparation)
}
\references{
http://www.biotechworld.it/bioinf/2016/03/09/analyzing-fluoresence-microscopy-data-with-r/
}
\examples{
data(leio_cells_dapi)
image(leio_cells_dapi[400:900,500:1000])
}

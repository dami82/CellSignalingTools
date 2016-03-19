\name{leio_basal_profile}
\alias{leio_basal_profile}
\docType{data}
\title{
Immunofluorescence Signal Profiles of control Leiomyoma Cells
}
\description{
Position and Signal Intensity Values of a single cell from an Immunofluorescence experiment performed on control Leiomyoma cells. Cells were immunostained with an anti-PTEN-Alexa488 antibody and then counterstained with DAPI. Stacked 8-bit images corresponding to the blue (B1) and the green (G1) channels were analyzed using the "line select" and the "plot profile" tools of ImageJ. 

}
\usage{data("leio_basal_profile")}
\format{
  A data frame with 470 observations on the following 3 variables.
  \describe{
    \item{\code{X1}}{a numeric vector; contains position information along x-axis}
    \item{\code{G1}}{a numeric vector; contains values of signal intensities from the green channel}
    \item{\code{B1}}{a numeric vector; contains values of signal intensities from the DAPI channel}
  }
}
\details{
Leiomyoma Cells were mock treated for 6 hours and then fixed with PFA. Fixed samples were stained with an anti-PTEN primary antibody and then counterstained with DAPI. For further information, please check the reference cited below. This dataset is part of an experiment including a Paraquat-treated sample as well ("leio_pq100_profile").
}
\source{
Vidimar V et al, 2016? (in preparation)
}
\references{
http://www.biotechworld.it/bioinf/2016/03/09/analyzing-fluoresence-microscopy-data-with-r/
}
\examples{
data(leio_basal_profile)
with(leio_basal_profile, plot(B1~X1, pch = 16, cex = 0.75))
with(leio_basal_profile, lines(G1~X1, col = "green", lwd = 3))
}

\name{scale_01}
\alias{scale_01}
\title{
Scale Matrices to a Zero to One Range
}
\description{
Scale Vectors or Matrices to a Zero to One Range. 
}
\usage{
scale_01(num_matrix)
}
\arguments{
  \item{num_matrix}{
is a numeric matrix or a numeric vector  
}
}
\value{
a numeric matrix or a numeric vector
}
\references{
http://www.biotechworld.it/bioinf/2016/03/09/analyzing-fluoresence-microscopy-data-with-r/
}
\author{
Damiano Fantini
}
\examples{
scale_01(c(-2,1,1,2,4,10))
}
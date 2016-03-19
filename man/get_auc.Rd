\name{get_auc}
\alias{get_auc}
\title{
Determine Area Under the Curve
}
\description{
Determine the area under the curve given a numeric matrix containing position values (x-axis) and signal intensities (y-axis).  
}
\usage{
get_auc(signal_matrix)
}
\arguments{
  \item{signal_matrix}{
is a numeric matrix containing two columns.  
}
}
\details{
The numeric matrix passed as argument must contain only two columns. The first column defines the position on the x-axis. The second column provides the signal intensity values (y-axis). The base is set at 0. 
}
\value{
Returns a numeric value corresponding to the area under the curve
}
\references{
http://www.biotechworld.it/bioinf/2016/03/09/analyzing-fluoresence-microscopy-data-with-r/
}
\author{
Damiano Fantini
}
\examples{
my_mat <- matrix(c(0,0.25,2,2.9,3,0,3,7,5,0), ncol = 2)
plot(my_mat, ylim = c(0,10), type = "n")
lines(my_mat, lwd = 3)
get_auc(my_mat)
}

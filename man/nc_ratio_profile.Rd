\name{nc_ratio_profile}
\alias{nc_ratio_profile}
\title{
Nuclear to Cytoplasmic Signal Ratio Analysis on Stacked Signal Profiles
}
\description{
Analyze Stacked Signal Profiles (calculated using plotProfile tool of ImageJ or similar software) and determine the nuclear to non-nuclear signal ratio. A plot displaying the signal intensity in the green and the DAPI channels is also displayed.
}
\usage{
nc_ratio_profile(profile_tab, 
bckgr_subtr = "baseline", 
mode = "auc", 
condition = "", 
antibody = "Ab - I", 
show_labels = TRUE, 
noise = 0.05)
}
\arguments{
  \item{profile_tab}{
is a data frame containing a X1 column (position along x), an one or more pairs of Gi, Bi columns (signal intensity in the Green channel and Blue channel along x). If profile_tab is not correctly formatted, the functionn will not work correctly.
}
  \item{bckgr_subtr}{
is one of the strings in c("baseline, "confocal""). It is the type of background subtraction that will be employed
}
  \item{mode}{
is one of the strings in c("auc", "mean"). It determines whether the nuclear/cyto signal ratio is calculated based on average signal or area-under-curve.
}
  \item{condition}{
is an optional string that will be used as title of the plot 
}
  \item{antibody}{
is a string that will be used as a label for the signal in the green channel 
}
  \item{show_labels}{
is a boolean and determines whether displaying extra information on the profile plots
}
  \item{noise}{
is a number in the range 0 to 1 and define the signal threshold to be used for background correction
}
}
\value{
list of numeric vectors. Each element of the list corresponds to one of the cells analyzed. Each numeric vector includes values for nuclear signal, non-nuclear signal and ratio nucl / non-nucl.
}
\references{
http://www.biotechworld.it/bioinf/2016/03/09/analyzing-fluoresence-microscopy-data-with-r/
}
\author{
Damiano Fantini
}
\examples{
data("leio_basal_profile")
data("leio_pq100_profile")

items <- list(leio_basal_profile, leio_pq100_profile)
conditions <- c('CTR','PQ100')
results <- lapply(1:2, (function(n){
  nc_ratio_profile(items[[n]], 
  bckgr_subtr = "baseline", 
  mode = "mean", 
  condition = conditions[n], 
  antibody = "PTEN", 
  show_labels = TRUE, 
  noise = 0.05)
})) 
}

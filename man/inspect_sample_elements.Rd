% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaprotr_inspect_sample_elements.R
\name{inspect_sample_elements}
\alias{inspect_sample_elements}
\title{inspect_sample_elements}
\usage{
inspect_sample_elements(spectral_count_object, force = FALSE)
}
\arguments{
\item{spectral_count_object}{List defined as "spectral_count_object" 
containing dataframes with abundance expressed as spectral counts. 
The spectral data can be organized by peptides, subgroups, 
groups or taxonomic levels.}

\item{force}{Logic value set at FALSE by default in order to ask 
permission to create a pdf file in the workstation of the user.}
}
\value{
Barplots (pdf) ilustrating the common spectral elements 
   (peptides, subgroups, groups, taxonomic elements) per sample
   in a "spectral_count_object".
}
\description{
Displays a graph that indicates the number of common elements from a 
"spectral_count_object" (peptides, subgroups, groups or taxonomic entities) 
per sample. This function is useful to distinguish heterogeneity between 
samples in an experimental design.
}
\examples{
\dontshow{.old_wd <- setwd(tempdir())}

data(fecal_waters)
inspect_sample_elements(fecal_waters)

\dontshow{setwd(.old_wd)}
}

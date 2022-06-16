\name{GE600}
\alias{GE600}
\alias{case.dat}
\alias{ctrl.dat}

\docType{data}
\title{
ISLET deconvolution example raw input data
}
\description{
\code{GE600} contains the raw example datasets for ISLET. It has the gene expression values for 50 genes among 200 subjects, with 3 temporal measures per subject. Among these 200 subjects, 100 are cases (i.e. cancer) and 100 are controls (i.e. healthy). Data were separated by case/control status, into two \code{SummarizedExperiment} objects. These raw example datasets will need to be converted by the \code{dataprep} function, and then they will be ready for downstream deconvolution (function \code{islet.solve}) and/or differentially expressed gene testing (function \code{islet.test}).
}
\usage{data(GE600)}
\format{
  Two \code{SummarizedExperiment} objects containing the following elements, respectively for each object:
  \describe{
    \item{\code{counts}}{A gene expression value data frame of 50 genes by 300 sample. 100 cancer/control subjects times 3 repeated measurement each subject.}
    \item{\code{colData}}{Sample meta-data. First column is the subject ID, shows the relationship between the 300 samples IDs and their 100 subject IDs. The remaining 6 columns (i.e. column 2-7) are the cell type proportions of all samples by their 6 cell types.}
  }
}

\examples{
data(GE600)
ls()
## [1] "case.dat" "ctrl.dat"

#show case.dat as an example
case.dat
#An object of class "SummarizedExperiment"

dim(assays(case.dat)$counts)
##[1]  50 300

dim(colData(case.dat))
##[1] 300   7

#Then, we can proceed to data preparation step, function 'dataprep' for ISLET.
##The rest of the deconvolution/csDE analysis will then follow.
}
\keyword{datasets}

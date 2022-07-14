\name{GE600age}
\alias{GE600age}
\alias{case.age.dat}
\alias{ctrl.age.dat}

\docType{data}
\title{
ISLET example datasets for slope variable testing in csDE
}
\description{
\code{GE600age} contains the example input datasets for ISLET's slope testing function. It has the gene expression values for 50 genes among 200 subjects, with 3 temporal measures per subject. Among these 200 subjects, 100 are cases (i.e. cancer) and 100 are controls (i.e. healthy). Temporal measures are at different age, with an \code{age} variable stored in the metadata. This is the main variable-of-interest in downstream testing. Data were separated by case/control status, into two \code{SummarizedExperiment} objects. These example datasets will need to be converted by the \code{dataprep} function, and then they will be ready for downstream deconvolution (function \code{islet.solve}) and/or differentially expressed gene testing (function \code{islet.test.slope}).


}
\usage{data(GE600age)}
\format{
  Two \code{SummarizedExperiment} objects containing the following elements, respectively for each object:
  \describe{
    \item{\code{counts}}{A gene expression value data frame of 50 genes by 300 sample. 100 cancer/control subjects times 3 repeated measurement each subject.}
    \item{\code{colData}}{Sample meta-data. First column is the subject ID, shows the relationship between the 300 samples IDs and their 100 subject IDs. The second column is the age variable for each sample, which is the main variable in downstream testing. The remaining 6 columns (i.e. column 3-8) are the cell type proportions of all samples by their 6 cell types.}
  }
}

\examples{
data(GE600age)
ls()
## [1] "case.age.dat" "ctrl.age.dat"

#show case.dat as an example
case.age.dat
#An object of class "SummarizedExperiment"

dim(assays(case.age.dat)$counts)
##[1]  50 300

dim(colData(case.age.dat))
##[1] 300   8

#Then, we can proceed to data preparation step, function 'dataprep' for ISLET.
##The rest of the csDE testing on age(slope) effect will then follow.
}
\keyword{datasets}

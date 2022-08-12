\name{GE600age}
\alias{GE600age}
\alias{GE600age_se}


\docType{data}
\title{
ISLET example datasets for slope variable testing in csDE
}
\description{
\code{GE600age} contains the example input datasets for ISLET's slope testing function. It has the gene expression values for 50 genes among 200 subjects, with 3 temporal measures per subject. Among these 200 subjects, 100 are cases (i.e. cancer) and 100 are controls (i.e. healthy). Temporal measures are at different age, with an \code{age} variable stored in the metadata. This is the main variable-of-interest in downstream testing. Data were combined by case/control status, into one \code{SummarizedExperiment} object. These example datasets will need to be converted by the \code{dataPrep} function, and then they will be ready for downstream deconvolution (function \code{isletSolve}) and/or differentially expressed gene testing (function \code{isletTest}).


}
\usage{data(GE600age)}
\format{
  One \code{SummarizedExperiment} object containing the following elements:
  \describe{
    \item{\code{counts}}{A gene expression value data frame of 50 genes by 600 sample. 200 cancer/control subjects times 3 repeated measurement each subject.}
    \item{\code{colData}}{Sample meta-data. The first column is the case/ctrl group status, the second column is the subject ID, shows the relationship between the 600 samples IDs and their 200 subject IDs. The third column is the age variable for each sample, which is the main variable in downstream testing. The remaining 6 columns (i.e. column 4-9) are the cell type proportions of all samples by their 6 cell types.}
  }
}

\examples{
data(GE600age)
ls()
## [1] "GE600age_se"

#show GE600age_se details
GE600age_se
#An object of class "SummarizedExperiment"

dim(assays(GE600age_se)$counts)
##[1]  50 600

dim(colData(GE600age_se))
##[1] 600   9

#Then, we can proceed to data preparation step, function 'dataPrep' for ISLET.
##The rest of the csDE testing on age(slope) effect will then follow.
}
\keyword{datasets}

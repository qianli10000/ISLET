\name{dataPrep}
\alias{dataPrep}
\title{Data preparation step utility function before using ISLET}
\description{
This function, \code{dataPrep}, is a necessary step to make your data ready and acceptable to ISLET. It takes one \code{SummarizedExperiment} objects listed below as the input, and produce a list ready to be feed into ISLET downstream deconvolution (function \code{isletSolve}) and/or cell type-specific differentially expressed gene testing (function \code{isletTest}).
}
\usage{
dataPrep(dat_se)
}

\arguments{
  \item{dat_se}{
A \code{SummarizedExperiment} object containing both case group and control. It should contain the observed feature by sample matrix, group label (case or ctrl) for each sample, subject IDs for each sample, and cell type proportions for each sample. The feature by sample matrix is stored in the assay, while the remaining components are stored in the colData slot. In the colData slot, the group label should be stored in the first column, subject IDs should be stroed in the second column, and the cell type proportions take the remaining columns. Subject IDs need to be sorted before using the \code{dataPrep} function.
}
}

\details{
This is the initial step for using ISLET, to prepare your data input ready for downstream deconvolution (function \code{isletSolve}) and/or differentially expressed gene testing (function \code{isletTest}). The input data must follow requirements listed above.
}

\value{
\code{dataPrep} returns a S4 object, containing elements ready to serve as the input for downstream deconvolution (function \code{isletSolve}) and/or differentially expressed gene testing (function \code{isletTest}).
}

\author{
Hao Feng <hxf155@case.edu>
}

\examples{
data(GE600)
ls()
## [1] "GE600_se"

study123input <- dataPrep(dat_se=GE600_se)
}

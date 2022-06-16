\name{dataprep}
\alias{dataprep}
\title{Data preparation step utility function before using ISLET}
\description{
This function, \code{dataprep}, is a necessary step to make your data ready and acceptable to ISLET. It takes two \code{SummarizedExperiment} objects listed below as the input, and produce a list ready to be feed into ISLET downstream deconvolution (function \code{islet.solve}) and/or differentially expressed gene testing (function \code{islet.test}).
}
\usage{
dataprep(case_dat_se, ctrl_dat_se)
}

\arguments{
  \item{case_dat_se}{
A \code{SummarizedExperiment} object for case group. It should contain the observed feature by sample matrix, subject IDs for each sample, and cell type proportions for each sample. The feature by sample matrix is stored in the assay, while the remaining components are stored in the colData slot. In the colData slot, the subject IDs take the first column, and the cell type proportions take the remaining columns.
}
  \item{ctrl_dat_se}{Similar to \code{case_dat_se}, but for control group. Requirements are the same.
}
}

\details{
This is the initial step for using ISLET, to prepare your data input ready for downstream deconvolution (function \code{islet.solve}) and/or differentially expressed gene testing (function \code{islet.test}). The input data matrices must follow requirements listed above, and samples/subjects must be ordered and match across matrices.
}

\value{
\code{dataprep} returns a list of elements ready to serve as the input for downstream deconvolution (function \code{islet.solve}) and/or differentially expressed gene testing (function \code{islet.test}).
}

\author{
Hao Feng <hxf155@case.edu>
}

\examples{
data(GE600)
ls()
## [1] "case.dat" "ctrl.dat"

study123input = dataprep(case_dat_se = case.dat, ctrl_dat_se = ctrl.dat)
}

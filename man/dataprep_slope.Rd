\name{dataprep_slope}
\alias{dataprep_slope}
\title{Data preparation function specifically for slope effect testing}
\description{
This function, \code{dataprep_slope}, is a necessary step to make your data ready and acceptable to ISLET slope effect testing. It takes two \code{SummarizedExperiment} objects listed below as the input, and produce a list ready to be feed into ISLET downstream deconvolution (function \code{islet.solve}) and/or differentially expressed gene testing (function \code{islet.test}).
}
\usage{
dataprep_slope(case_dat_se, ctrl_dat_se)
}

\arguments{
  \item{case_dat_se}{
A \code{SummarizedExperiment} object for case group. It should contain the observed feature by sample matrix, subject IDs for each sample, a slope variable for each sample, and cell type proportions for each sample. The feature by sample matrix is stored in the assay, while the remaining components are stored in the colData slot. In the colData slot, the subject IDs take the first column, the slope variable takes the second column, and the cell type proportions take the remaining columns.
}
  \item{ctrl_dat_se}{Similar to \code{case_dat_se}, but for control group. Requirements are the same.
}
}

\details{
This is the initial step for using ISLET, to prepare your data input ready for downstream cell-type-specific differentially expressed gene testing (function \code{islet.test}) with respect to the slope variable. The input data matrices must follow requirements listed above, and samples/subjects must be ordered and match across matrices.
}

\value{
\code{dataprep_slope} returns a list of elements ready to serve as the input for cell-type-specific differentially expressed gene testing (function \code{islet.test}).
}

\author{
Hao Feng <hxf155@case.edu>
}


\examples{
data(GE600age)
ls()
## [1] "case.age.dat" "ctrl.age.dat"

#(1) Data preparation
study456input = dataprep_slope(case_dat_se = case.age.dat, ctrl_dat_se = ctrl.age.dat)

#(2) Test for slope effect(i.e. age) difference in csDE testing
result.test = islet.test(input = study456input)

}

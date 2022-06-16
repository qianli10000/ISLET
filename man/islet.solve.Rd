\name{islet.solve}
\alias{islet.solve}
\title{
Solving individual-specific and cell-type-specific reference panels
}
\description{
This function, \code{islet.solve}, is a core function of ISLET. It takes the output from data preparation function \code{dataprep}, and solve for individual-specific and cell-type-specific reference panels. It has parallel computing implemented to speed up the EM algorithm application.
}
\usage{
islet.solve(input)
}

\arguments{
  \item{input}{
The list object output from data preparation function \code{dataprep}.
}
}
\details{
The deconvolution result output is a list of length \emph{G}, where \emph{G} is the number of features (i.e. genes). For each of the \emph{G} elements, it is a list of length 7. The description of the 7 components is listed below.
}
\value{
\item{case.m}{A vector of length \emph{K}, where \emph{K} is the number of cell types. It contains the mean cell type profiles (i.e. mean cell type gene expression values for \emph{K} cell types), for case group, on a specific feature. }
\item{ctrl.m}{A vector of length \emph{K}, where \emph{K} is the number of cell types. It contains the mean cell type profiles (i.e. mean cell type gene expression values for \emph{K} cell types), for control group, on a specific feature.}
\item{case.indv}{A subject by cell type matrix containing all the feature values (i.e. gene expression values), for case group. It is one of the main products the individual-specific and cell-type-specific solve algorithm. Values are the summation of fixed effect (\code{case.m}) and individual random effect (not shown explicitly).}
\item{ctrl.indv}{A subject by cell type matrix containing all the feature values (i.e. gene expression values), for control group. It is one of the main products the individual-specific and cell-type-specific solve algorithm. Values are the summation of fixed effect (\code{ctrl.m}) and individual random effect (not shown explicitly).}
\item{var.k}{A vector of length \emph{K}, where \emph{K} is the number of cell types. It contains the estimated feature value variance for the \emph{K} cell types, on a specific feature.}
\item{var.0}{A scalar. The estimated residual variance, on a specific feature.}
\item{LLK}{A scalar. The log-likelihood from the current model. It can be useful for testing purpose such as Likelihood Ratio Test.}
}

\author{
Hao Feng <hxf155@case.edu>
}


\examples{
data(GE600)
ls()
## [1] "case.dat" "ctrl.dat"

#(1) Data preparation
study123input = dataprep(case_dat_se = case.dat, ctrl_dat_se = ctrl.dat)

#(2) Individual-specific and cell-type-specific deconvolution
result = islet.solve(input = study123input)
}

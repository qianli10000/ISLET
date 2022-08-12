\name{isletSolve}
\alias{isletSolve}
\title{
Solving individual-specific and cell-type-specific reference panels
}
\description{
This function, \code{isletSolve}, is a core function of ISLET. It takes the output from data preparation function \code{dataPrep}, and solve for individual-specific and cell-type-specific reference panels. It has parallel computing implemented to speed up the EM algorithm application.
}
\usage{
isletSolve(input, ncores = min(detectCores()-1, 15) )
}

\arguments{
  \item{input}{
The list object output from data preparation function \code{dataPrep}.
}
  \item{ncores}{
The number of computing cores requested for parallel computing to speed up the estimation. If not specified, the default will learn from the computing environment.
}

}
\details{
For both case group and control group, the deconvolution result is a list of length \emph{K}, where \emph{K} is the number of cell types. For each of the \emph{K} elements, it is a matrix of dimension \emph{G} by \emph{N}. It stores the deconvoluted feature (\emph{G}) by subject (\emph{N}) values, for each of the \emph{K} elements.
}
\value{
\item{case.ind.ref}{A list of length \emph{K}, where \emph{K} is the number of cell types. For each of the \emph{K} elements in this list, it is a feature by subject matrix containing all the feature values (i.e. gene expression values), for case group. It is one of the main products the individual-specific and cell-type-specific solve algorithm. }
\item{ctrl.ind.ref}{A list of length \emph{K}, where \emph{K} is the number of cell types. For each of the \emph{K} elements in this list, it is a feature by subject matrix containing all the feature values (i.e. gene expression values), for control group. It is one of the main products the individual-specific and cell-type-specific solve algorithm. }
\item{mLLK}{A scalar. The log-likelihood from the current model. It can be useful for testing purpose such as Likelihood Ratio Test.}
}

\author{
Hao Feng <hxf155@case.edu>
}


\examples{
data(GE600)
ls()
## [1] "GE600_se"

#(1) Data preparation
study123input <- dataPrep(dat_se=GE600_se)

#(2) Individual-specific and cell-type-specific deconvolution
result <- isletSolve(input=study123input)
}

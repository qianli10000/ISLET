\name{ctrlEst}
\alias{ctrlEst}

\title{ISLET reference estimation results extraction}
\description{
This function, \code{ctrlEst}, extracts the estimated reference panels from the control group. It takes one \code{outputSol} objects as the input, and produce a list of matrices containing the estimated reference panels. The length of the list is equal to the number of cell types. For each matrix, it contains the gene by subject reference panel for that specific cell type.
}
\usage{
ctrlEst(res.sol)
}

\arguments{
  \item{res.sol}{
A  \code{outputSol} objects from \code{isletSolve} step.
}
}

\details{
This is the accessor function help to extract the estimated reference panels from \code{isletSolve} step.
}

\value{
A list of matrices containing the estimated reference panels. The length of the list is equal to the number of cell types. For each matrix, it contains the gene by subject reference panel for that specific cell type.
}

\author{
Hao Feng <hxf155@case.edu>
}

\examples{
data(GE600)

study123input <- dataPrep(dat_se=GE600_se)
res.sol <- isletSolve(input=study123input)
caseVal <- ctrlEst(res.sol)

}

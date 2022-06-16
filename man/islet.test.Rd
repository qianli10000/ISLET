\name{islet.test}
\alias{islet.test}
\title{
Testing for cell-type-specific Differential Expression (csDE) genes
}
\description{
This function, \code{islet.test}, can take the output from data preparation function \code{dataprep}, and test for csDE genes. It uses Likelihood Ratio Test (LRT), iterating all cell types. The output is a matrix of p-values from LRT. It has parallel computing implemented to speed up the EM algorithm application.
}
\usage{
islet.test(input)
}

\arguments{
  \item{input}{
The list object output from data preparation function \code{dataprep}.
}
}
\details{
This function implement a LRT, and run individually for each cell type, and then aggregate the results together into a matrix.
}
\value{
 \item{}{A p-value matrix, in the dimension of feature by cell type. Each element is the LRT p-value, by contrasting case group and control group, for one feature in one cell type.}
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

#(2) [optional for csDE genes testing] Individual-specific and cell-type-specific deconvolution
result.solve = islet.solve(input = study123input)

#(3) Test for csDE genes
result.test = islet.test(input = study123input)
}

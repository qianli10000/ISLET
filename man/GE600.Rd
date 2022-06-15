\name{GE600}
\alias{GE600}
\alias{case_obs}
\alias{ctrl_obs}
\alias{case_prop}
\alias{ctrl_prop}
\alias{case_subject}
\alias{ctrl_subject}


\docType{data}
\title{
ISLET deconvolution example raw input data
}
\description{
\code{GE600} contains the raw example datasets for ISLET. It has the gene expression values for 50 genes among 200 subjects, with 3 temporal measures per subject. These raw example datasets will need to be converted by the \code{dataprep} function, and then they will be ready for downstream deconvolution (function \code{islet.solve}) and/or differentially expressed gene testing (function \code{islet.test}).
}
\usage{data(GE600)}
\format{
  6 dataframes/matrices of the following elements:
  \describe{
    \item{\code{case_obs}}{A gene expression value data frame of 50 genes by 300 cancer sample, for the cancer group.}
    \item{\code{ctrl_obs}}{A gene expression value data frame of 50 genes by 300 control sample, for the control group.}
    \item{\code{case_prop}}{A cell proportion matrix of 300 cancer samples by 6 cell types, for the cancer group.}
    \item{\code{ctrl_prop}}{A cell proportion matrix of 300 control samples by 6 cell types, for the control group.}
    \item{\code{case_subject}}{A two-column data frame showing the relationship between the 300 cancer samples and their 100 subjects, for the cancer group. The first column is the \emph{sample_id}. The second column is the \emph{subject_id}. Each subject has 3 repeated/temporal measurements.}
    \item{\code{ctrl_subject}}{A two-column data frame showing the relationship between the 300 control samples and their 100 subjects, for the control group. The first column is the \emph{sample_id}. The second column is the \emph{subject_id}. Each subject has 3 repeated/temporal measurements.}

  }
}

\examples{
data(GE600)
ls()
## "case_obs"     "case_prop"    "case_subject" "ctrl_obs"     "ctrl_prop"    "ctrl_subject"


#Then, we can proceed to data preparation step for ISLET.
#(1) Data preparation
study123input = dataprep(case_obs = case_obs, ctrl_obs = ctrl_obs,
         case_prop = case_prop, ctrl_prop = ctrl_prop,
         case_subject = case_subject, ctrl_subject = ctrl_subject)
##The rest of the analysis will follow.
}
\keyword{datasets}

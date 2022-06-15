\name{dataprep}
\alias{dataprep}
\title{Data preparation step utility function before using ISLET}
\description{
This function, \code{dataprep}, is a necessary step to make your data ready and acceptable to ISLET. It takes 6 elements listed below as the input, and produce a list ready to be feed into ISLET downstream deconvolution (function \code{islet.solve}) and/or differentially expressed gene testing (function \code{islet.test}).
}
\usage{
dataprep(case_obs, ctrl_obs,
         case_prop, ctrl_prop,
         case_subject, ctrl_subject)
}

\arguments{
  \item{case_obs}{
A feature by sample data frame (i.e. gene by sample), for case group (i.e. cancer samples), to deconvolute.
}
  \item{ctrl_obs}{
A feature by sample data frame (i.e. gene by sample), for control group (i.e. healthy samples), to deconvolute. The features (i.e. rows) must match features in \code{case_obs}.
}
  \item{case_prop}{
A sample by cell type matrix containing cell type proportions, for case group (i.e. cancer samples). The samples (i.e. rows) must match samples (i.e. columns) in \code{case_obs}.
}
  \item{ctrl_prop}{
A sample by cell type matrix containing cell type proportions, for control group (i.e. healthy samples). The samples (i.e. rows) must match samples (i.e. columns) in \code{ctrl_obs}.
}
  \item{case_subject}{
A two-column data frame showing the relationship between samples and subjects, for case group (i.e. cancer samples). The first column is the \emph{sample_id}, which must match the samples (i.e. columns) in \code{case_obs}. The second column is the \emph{subject_id}, which shows the subject ID associated with each sample. Naturally, \emph{subject_id} can be repetitive while \emph{sample_id} cannot.
}
  \item{ctrl_subject}{
A two-column data frame showing the relationship between samples and subjects, for control group (i.e. healthy samples). The first column is the \emph{sample_id}, which must match the samples (i.e. columns) in \code{case_obs}. The second column is the \emph{subject_id}, which shows the subject ID associated with each sample. The subject IDs across case and control groups must be unique.
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
## "case_obs"     "case_prop"    "case_subject" "ctrl_obs"     "ctrl_prop"    "ctrl_subject"

study123input = dataprep(case_obs = case_obs, ctrl_obs = ctrl_obs,
         case_prop = case_prop, ctrl_prop = ctrl_prop,
         case_subject = case_subject, ctrl_subject = ctrl_subject)
}

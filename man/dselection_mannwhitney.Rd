% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dselection_MannWhitney.R
\name{dselection_MannWhitney}
\alias{dselection_MannWhitney}
\title{dselection_MannWhitney}
\usage{
dselection_MannWhitney(S_Y, S_X, S, alpha)
}
\arguments{
\item{S_Y}{score vector for the test set}

\item{S_X}{score vector for the calibration set}

\item{S}{set of the selected indices}

\item{alpha}{significance level}
}
\value{
An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for
the number of true discoveries in closed testing procedure using
Wilcoxon-Mann-Whitney local test applied to conformal \emph{p}-values.
The selection set, i.e., the set of hypothesis
indices that we are interested in can be chosen.
}
\description{
dselection_MannWhitney
}
\examples{
set.seed(321)
Sxy = sample(x=1:1000, size=100)
Sx = sample(Sxy, size=70)
Sy = setdiff(Sxy, Sx)
dselection_MannWhitney(S_Y=Sy, S_X=Sx, alpha = 0.1, S=1:10)




}

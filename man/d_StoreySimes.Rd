% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/d_StoreySimes.R
\name{d_StoreySimes}
\alias{d_StoreySimes}
\title{d_StoreySimes}
\usage{
d_StoreySimes(S_Y, S_X, alpha = 0.1, lambda = 0.5)
}
\arguments{
\item{S_Y}{: score vector of test observations}

\item{S_X}{: score vector of calibration observations}

\item{alpha}{: significance level of the local test. Default value is set equal to 0.1}

\item{lambda}{: parameter involved in the computation of Storey estimator. Default value is set equal to 0.5}
}
\value{
An integer which is the \eqn{(1 − \alpha)}-confidence lower bound for the
number of true discoveries in closed testing procedure using Simes local
test with Storey's estimator for the proportion of true null hypotheses applied to conformal \emph{p}-values.
The selection set, i.e. the set of hypothesis indices that we are
interested in is \eqn{[m]=:{1,...,m}} by default.
Then, Storey estimator is computed as
\deqn{\hat\pi_0 = \frac{1+\sum_{i=1}^m \mathbb{1}\{p_i>\lambda\}}{m(1-\lambda)}}
where \eqn{\lambda\in(0,1)} and \eqn{p_i} is the \emph{p}-value related
to the null hypothesis \eqn{H_i, \hspace{2mm} i=1,\ldots,m}.
}
\description{
It returns the lower bound for the number of true discoveries in closed testing procedure
using Simes local test with Storey estimator for the proportion of true
null hypotheses.
}
\examples{
set.seed(321)
Sxy = sample(x=1:1000, size=100)
Sx = sample(Sxy, size=70)
Sy = setdiff(Sxy, Sx)
d_StoreySimes(S_Y=Sy, S_X=Sx)



}

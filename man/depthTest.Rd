% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/depthTest.R
\name{depthTest}
\alias{depthTest}
\title{Test for a difference in distribution}
\usage{
depthTest(sample1, sample2, method, depthFun = NULL)
}
\arguments{
\item{sample1}{A matrix or data frame containing the first sample}

\item{sample2}{A matrix or data frame containing the second sample}

\item{method}{The depth method to use (one of "halfspace", "mahalanobis", "pvb", "lp", "pd", "lcd", or "custom")}

\item{depthFun}{A depth function which takes two datasets and returns the depths}
}
\value{
The p-value and test statistic for a test of difference in distribution
}
\description{
Test for a difference in distribution
}
\examples{
library(mvtnorm)
samp1 <- rmvnorm(100, mean = rep(0, 10))
samp2 <- rmvnorm(100, mean = rep(1, 10))
depthTest(samp1, samp2, "mahalanobis")
}
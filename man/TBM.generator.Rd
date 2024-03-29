% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HLloyd.R
\name{TBM.generator}
\alias{TBM.generator}
\title{Generate a tensor block model under permutation and noises.}
\usage{
TBM.generator(p, r, delta = 1, sigma = 1)
}
\arguments{
\item{p}{A vector of tensor dimensions.}

\item{r}{A vector of cluster numbers in each mode.}

\item{delta}{The slice separation of core tensor.}

\item{sigma}{Observational noise level.}
}
\value{
Observational tensor; clutering labels; signal tensor.
}
\description{
Generate a tensor block model under permutation and noises.
}

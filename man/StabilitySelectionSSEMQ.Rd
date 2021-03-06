% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solver.R
\name{StabilitySelectionSSEMQ}
\alias{StabilitySelectionSSEMQ}
\title{StabilitySelectionSSEMQ}
\usage{
StabilitySelectionSSEMQ(
  X,
  Y,
  B,
  F,
  Sk,
  sigma2,
  lambda,
  rho,
  Wb,
  Wf,
  p,
  R = 100,
  threshold = 1e-04
)
}
\arguments{
\item{X}{eQTL matrix}

\item{Y}{Gene expression matrix}

\item{B}{initialized GRN-matrix}

\item{F}{initialized eQTL effect matrix}

\item{Sk}{eQTL index of genes}

\item{sigma2}{initialized noise variance from ridge regression}

\item{lambda}{Optimized selected hyperparameter of lasso term in SSEM for B from opt.SSEMiPALM}

\item{rho}{Optimized selected Hyperparameter of lasso term in SSEM for F from opt.SSEMiPALM}

\item{Wb}{weight matrix of B for adaptive lasso terms. Default as 1 / B(init) from eidge regression.}

\item{Wf}{weight matrix of F for adaptive lasso terms. Default as 1 / F(init) from ridge regression.}

\item{p}{number of genes}

\item{R}{number of replicates in bootstrap of stability selection}

\item{threshold}{threshold of estimates}
}
\value{
list Matrices for B, F, (I-B)^{-1}F, entries with probability of nonzeros within R replicates
}
\description{
StabilitySelectionSSEMQ
}

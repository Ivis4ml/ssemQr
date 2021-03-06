% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generator.R
\name{randomeQTLdata}
\alias{randomeQTLdata}
\title{randomeQTLdata}
\usage{
randomeQTLdata(
  n,
  p,
  k,
  sparse = 1,
  sqtl = 0.5,
  sigma2 = 0.01,
  intercept = 5,
  type = c("DG", "ER"),
  dag = TRUE,
  coefs = c(0.2, 0.4),
  esize = c(0.5, 1),
  overlap = c("none", "weak", "strong"),
  qdist = 100,
  ncore = 0.1,
  span = FALSE,
  noqtl = FALSE,
  rmv = 0.1
)
}
\arguments{
\item{n}{number of observations}

\item{p}{number of candidate genes}

\item{k}{number of candidate eQTLs (nonzero)}

\item{sparse}{expected number of nonzero edges of each gene}

\item{sqtl}{sparsity of cis-eQTL}

\item{esize}{effect size of eQTL}

\item{overlap}{SNP effects are overlapped or not. Default as FALSE}

\item{qdist}{average distance between eQTLs}

\item{span}{span SNP average or not}

\item{noqtl}{several genes have no cis-eQTL. Default FALSE}

\item{rmv}{ratio of genes have no cis-eQTL}
}
\description{
generate simulated network structure for eQTLSEM algorithm, with defined
number of genes and number of candidate eQTLs.
require("eQTLsem")
}
\examples{
data = randomeQTLdata(n = 100, p = 10, k = 30, type = "DG")
}

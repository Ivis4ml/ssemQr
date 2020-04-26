##' generate simulated network structure for eQTLSEM algorithm, with defined
##' number of genes and number of candidate eQTLs.
##' require("eQTLsem")
##' @title randomeQTLdata
##' @param n number of observations
##' @param p number of candidate genes
##' @param k number of candidate eQTLs (nonzero)
##' @param sparse expected number of nonzero edges of each gene
##' @param sqtl sparsity of cis-eQTL
##' @param esize effect size of eQTL
##' @param overlap SNP effects are overlapped or not. Default as FALSE
##' @param qdist average distance between eQTLs
##' @param span span SNP average or not
##' @param noqtl several genes have no cis-eQTL. Default FALSE
##' @param rmv ratio of genes have no cis-eQTL
##' @export
##' @importFrom igraph graph_from_adjacency_matrix is.dag sample_pa as_adjacency_matrix
##' @importFrom mvtnorm rmvnorm
##' @importFrom qtl sim.map sim.cross pull.geno find.markerpos sim.geno makeqtl
##' @examples
##' data = randomeQTLdata(n = 100, p = 10, k = 30, type = "DG")
randomeQTLdata = function(n, p, k, sparse = 1, sqtl = 0.5, sigma2 = 0.01, intercept = 5.0,
                          type = c("DG", "ER"), dag = TRUE, coefs = c(0.2, 0.4), esize = c(0.5, 1.0),
                          overlap = c("none", "weak", "strong"), qdist = 100, ncore = 0.1,
                          span = FALSE, noqtl = FALSE, rmv = 0.1) {
  coef.min = coefs[1]
  coef.max = coefs[2]
  type     = match.arg(type)
  overlap  = match.arg(overlap)
  DG = function() {
    B = matrix(0, nrow = p, ncol = p)
    d = p * p
    # number of edges
    ne = rbinom(1, d, sparse / (p - 1))
    niter = 0
    # candidate edge indices
    ce = setdiff(seq(1, d), which(diag(p) != 0))
    while (sum(B) <= ne & niter < 2 * d) {
      ix = sample(ce, 1)
      B[ix] = TRUE
      if (dag) {
        g = graph_from_adjacency_matrix(B)
        B[ix] = is.dag(g)
      }
      niter = niter + 1
    }
    ## simulated edge indices
    ei = which(B != 0)
    B[ei] = runif(length(ei), min = coef.min, max = coef.max) * sample(c(-1, 1), length(ei), replace = T)
    ## check matrix (invertible)
    Det = abs(det(diag(p) - B))
    if(Det > 1e-2 || Matrix::norm(B, "f") < 1) {
      Matrix::Matrix(B, sparse = T)
    } else {
      NULL
    }
  }
  ER = function() {
    B = as_adjacency_matrix(sample_pa(p, power = sparse), names = FALSE, sparse = FALSE)
    ord = sample(p)
    B = B[ord, ord]
    ## simulated edge indices
    ei = which(B != 0)
    B[ei] = runif(length(ei), min = coef.min, max = coef.max) * sample(c(-1, 1), length(ei), replace = T)
    ## check matrix (invertible)
    Det = abs(det(diag(p) - B))
    if(Det > 1e-2 || Matrix::norm(B, "f") < 1) {
      Matrix::Matrix(B, sparse = T)
    } else {
      NULL
    }
  }
  ## eQTL mapping configuration
  ## simulate F2 family by qpgraph. Therefore, for n genes, we simulate candidate eQTls
  ## for simplicity, we put all eQtls on same chromosome, with even spacing.
  eQTLs = function(qlen = 100) {
    nq = k / p
    # for each genes, candidate eQtls
    cq = nq / sqtl
    # simulated genetic map
    map = sim.map(len = qlen * p, n.mar = cq * p, eq.spacing = TRUE, include.x = FALSE)
    # simulate f2 cross genotype
    f2  = sim.cross(map, type = "f2", n.ind = n)
    X = t(pull.geno(f2)) - 1
    F = matrix(0, nrow = p, ncol = cq * p)
    Sk = lapply(1:p, function(i) {
      s = seq((i-1) * cq, i * cq - 1) + 1
      ## v = sample(seq(1, cq - nq + 1), 1)
      if (!span) {
        F[i, s[1:nq]] <<- runif(nq, min = esize[1], max = esize[2])
      } else {
        F[i, s[round(seq(1, cq, length.out = nq))]] <<- runif(nq, min = esize[1], max = esize[2])
      }
      s
    })
    list(F = F, X = X, Sk = Sk, Cross = f2, map = map)
  }
  nQTLs = function(qlen = 100) {
    nq = k / p
    # for each genes, candidate eQtls
    cq = nq / sqtl
    # simulated genetic map
    map = sim.map(len = qlen * p, n.mar = cq * p, eq.spacing = TRUE, include.x = FALSE)
    # simulate f2 cross genotype
    f2  = sim.cross(map, type = "f2", n.ind = n)
    X = t(pull.geno(f2)) - 1
    F = matrix(0, nrow = p, ncol = cq * p)
    Sk = lapply(1:p, function(i) {
      s = seq((i-1) * nq, (i-1) * nq + cq - 1) + 1
      ## v = sample(seq(1, cq - nq + 1), 1)
      if (!span) {
        F[i, s[1:nq]] <<- runif(nq, min = esize[1], max = esize[2])
      } else {
        F[i, s[round(seq(1, cq, length.out = nq))]] <<- runif(nq, min = esize[1], max = esize[2])
      }
      s
    })
    list(F = F, X = X, Sk = Sk, Cross = f2, map = map)
  }
  sQTLs = function(qlen = 100) {
    nq = k / p
    # for each genes, candidate eQtls
    cq = nq / sqtl
    # simulated genetic map
    map = sim.map(len = qlen * p, n.mar = cq * p, eq.spacing = TRUE, include.x = FALSE)
    # simulate f2 cross genotype
    f2  = sim.cross(map, type = "f2", n.ind = n)
    X = t(pull.geno(f2)) - 1
    F = matrix(0, nrow = p, ncol = cq * p)
    Sk = lapply(1:p, function(i) {
      s = seq(i, i + cq - 1)
      ## v = sample(seq(1, cq - nq + 1), 1)
      if (!span) {
        F[i, s[1:nq]] <<- runif(nq, min = esize[1], max = esize[2])
      } else {
        F[i, s[round(seq(1, cq, length.out = nq))]] <<- runif(nq, min = esize[1], max = esize[2])
      }
      s
    })
    list(F = F, X = X, Sk = Sk, Cross = f2, map = map)
  }
  B = NULL
  while (is.null(B)) {
    if (type == "DG") {
      B = DG()
    } else {
      B = ER()
    }
  }
  QTL = if (overlap == "none") {
    eQTLs(qlen = qdist)
  } else if (overlap == "weak") {
    nQTLs(qlen = qdist)
  } else {
    sQTLs(qlen = qdist)
  }

  F = QTL$F
  if (noqtl) {
    rmix = sample(p, rbinom(1, p, rmv))
    F[rmix, ] = 0
  }
  X = QTL$X
  Sk = QTL$Sk
  E = sqrt(sigma2) * t(rmvnorm(n, mean = rep(0, p), sigma = diag(p)))
  U = tcrossprod(rnorm(p, 0, intercept), rep(1, n))
  Y = solve(diag(p) - B) %*% (F %*% X + U + E)
  rownames(X) = sapply(1:nrow(X), function(i){sprintf("rs_%05d", i)})
  rownames(Y) = sapply(1:nrow(Y), function(i){sprintf("g_%05d", i)})
  V = F %*% X
  SNR = sqrt(var(as.vector(V)) / var(as.vector(E)))
  Cross = QTL$Cross
  Cross$pheno = as.data.frame(t(Y))
  list(
    Data = list(X = X, Y = Y, Sk = Sk),
    Vars = list(B = B, F = Matrix(F, sparse = T), Mu = U, n = n, p = p, k = k, q = ncol(F), SNR = SNR),
    eQTL = list(Cross  = Cross, map = QTL$map)
  )
}

##' generate simulated network structure for eQTLSEM algorithm, with fixed GRN and eQTL mapping
##' require("eQTLsem")
##' @title eQTLdataExt
##' @param B Gene regulatory network matrix
##' @param F eQTL mapping matrix
##' @param sigma2 noise of observations
##' @export
eQTLdataExt = function(B, F, Sk, QTL, n, k, sigma2 = 0.01, intercept = 5.0) {
  p = nrow(B)
  f2  = sim.cross(QTL$map, type = "f2", n.ind = n)
  X = t(pull.geno(f2)) - 1
  E = sqrt(sigma2) * t(rmvnorm(n, mean = rep(0, p), sigma = diag(p)))
  U = tcrossprod(rnorm(p, 0, intercept), rep(1, n))
  Y = solve(diag(p) - B) %*% (as.matrix(F) %*% X + U + E)
  rownames(X) = sapply(1:nrow(X), function(i){sprintf("rs_%05d", i)})
  rownames(Y) = sapply(1:nrow(Y), function(i){sprintf("g_%05d", i)})
  V = F %*% X
  SNR = sqrt(var(as.vector(V)) / var(as.vector(E)))
  Cross = f2
  Cross$pheno = as.data.frame(t(Y))
  list(
    Data = list(X = X, Y = Y, Sk = Sk),
    Vars = list(B = B, F = F, Mu = U, n = n, p = p, k = k, q = ncol(F), SNR = SNR),
    eQTL = list(Cross  = Cross, map = QTL$map)
  )
}



##' Generate simulated network structure for eQTLSEM algorithm, with defined
##' number of genes and number of cis-eQTLs and trans-eQTLs. To be different
##' with `randomeQTLdata`, `randQtlsdata`'s eQTL contains all Qtls with
##' require("eQTLsem")
##' @title randQtlsdata
##' @param n number of observations
##' @param p number of candidate genes
##' @param k number of candidate eQTLs (cis-eQTL + trans-eQTL)
##' @param sparse expected number of nonzero edges of each gene
##' @param s expected number of eQTLs
##' @param esize effect size of eQTL
##' @export
randQtlsdata = function(n, p, k, sparse = 1, s = 3, sigma2 = 0.01, intercept = 5.0,
                        type = c("DG", "ER"), dag = TRUE, coefs = c(0.2, 0.4), esize = c(0.5, 1.0),qdist = 100, ncore = 0.1) {
  coef.min = coefs[1]
  coef.max = coefs[2]
  type     = match.arg(type)
  DG = function() {
    B = matrix(0, nrow = p, ncol = p)
    d = p * p
    # number of edges
    ne = rbinom(1, d, sparse / (p - 1))
    niter = 0
    # candidate edge indices
    ce = setdiff(seq(1, d), which(diag(p) != 0))
    while (sum(B) <= ne & niter < 2 * d) {
      ix = sample(ce, 1)
      B[ix] = TRUE
      if (dag) {
        g = graph_from_adjacency_matrix(B)
        B[ix] = is.dag(g)
      }
      niter = niter + 1
    }
    ## simulated edge indices
    ei = which(B != 0)
    B[ei] = runif(length(ei), min = coef.min, max = coef.max) * sample(c(-1, 1), length(ei), replace = T)
    ## check matrix (invertible)
    Det = abs(det(diag(p) - B))
    if(Det > 1e-2 || Matrix::norm(B, "f") < 1) {
      Matrix::Matrix(B, sparse = T)
    } else {
      NULL
    }
  }
  ER = function() {
    B = as_adjacency_matrix(sample_pa(p, power = sparse), names = FALSE, sparse = FALSE)
    ord = sample(p)
    B = B[ord, ord]
    ## simulated edge indices
    ei = which(B != 0)
    B[ei] = runif(length(ei), min = coef.min, max = coef.max) * sample(c(-1, 1), length(ei), replace = T)
    ## check matrix (invertible)
    Det = abs(det(diag(p) - B))
    if(Det > 1e-2 || Matrix::norm(B, "f") < 1) {
      Matrix::Matrix(B, sparse = T)
    } else {
      NULL
    }
  }
  ## eQTL mapping configuration
  ## simulate F2 family by qpgraph. Therefore, for n genes, we simulate candidate eQTls
  ## for simplicity, we put all eQtls on same chromosome, with even spacing.
  eQTLs = function(qlen = 100) {
    # simulated genetic map
    map = sim.map(len = qlen * p, n.mar = k, eq.spacing = TRUE, include.x = FALSE)
    # simulate f2 cross genotype
    f2  = sim.cross(map, type = "f2", n.ind = n)
    X = t(pull.geno(f2)) - 1
    F = matrix(0, nrow = p, ncol = k)
    Sk = lapply(1:p, function(i) {
      seq(1, k)
    })
    d = p * s
    ne = sample(seq(1, p * k), d)
    F[ne] = runif(d, min = esize[1], max = esize[2])
    list(F = F, X = X, Sk = Sk, Cross = f2, map = map)
  }
  B = NULL
  while (is.null(B)) {
    if (type == "DG") {
      B = DG()
    } else {
      B = ER()
    }
  }
  QTL = eQTLs(qlen = qdist)
  F = QTL$F
  X = QTL$X
  Sk = QTL$Sk
  E = sqrt(sigma2) * t(rmvnorm(n, mean = rep(0, p), sigma = diag(p)))
  U = tcrossprod(rnorm(p, 0, intercept), rep(1, n))
  Y = solve(diag(p) - B) %*% (F %*% X + U + E)
  rownames(X) = sapply(1:nrow(X), function(i){sprintf("rs_%05d", i)})
  rownames(Y) = sapply(1:nrow(Y), function(i){sprintf("g_%05d", i)})
  V = F %*% X
  SNR = sqrt(var(as.vector(V)) / var(as.vector(E)))
  Cross = QTL$Cross
  Cross$pheno = as.data.frame(t(Y))
  list(
    Data = list(X = X, Y = Y, Sk = Sk),
    Vars = list(B = B, F = Matrix(F, sparse = T), Mu = U, n = n, p = p, k = k, q = ncol(F), SNR = SNR),
    eQTL = list(Cross  = Cross, map = QTL$map)
  )
}


##' generate simulated network structure for eQTLSEM algorithm, with defined
##' number of genes and number of candidate eQTLs.
##' require("eQTLsem")
##' @title randomeQTLdata4H
##' @param n number of observations
##' @param p number of candidate genes
##' @param k number of candidate eQTLs (nonzero)
##' @param sparse expected number of nonzero edges of each gene
##' @param sqtl sparsity of cis-eQTL
##' @param esize effect size of eQTL
##' @param overlap SNP effects are overlapped or not. Default as FALSE
##' @param qdist average distance between eQTLs
##' @param span span SNP average or not
##' @param noqtl several genes have no cis-eQTL. Default FALSE
##' @param rmv ratio of genes have no cis-eQTL
##' @export
##' @importFrom igraph graph_from_adjacency_matrix is.dag sample_pa as_adjacency_matrix
##' @importFrom mvtnorm rmvnorm
##' @importFrom qtl sim.map sim.cross pull.geno find.markerpos sim.geno makeqtl
##' @examples
##' data = randomeQTLdata4H(n = 100, p = 10, k = 30, type = "DG")
randomeQTLdata4H = function(n, p, k, sparse = 1, sqtl = 0.5, sigma2 = 0.01, intercept = 5.0,
                            type = c("DG", "ER"), dag = TRUE, coefs = c(0.2, 0.4), esize = c(0.5, 1.0),
                            overlap = c("none", "weak", "strong"), qdist = 100, ncore = 0.1,
                            span = FALSE, noqtl = FALSE, rmv = 0.1) {
  coef.min = coefs[1]
  coef.max = coefs[2]
  type     = match.arg(type)
  overlap  = match.arg(overlap)
  DG = function() {
    B = matrix(0, nrow = p, ncol = p)
    d = p * p
    # number of edges
    ne = rbinom(1, d, sparse / (p - 1))
    niter = 0
    # candidate edge indices
    ce = setdiff(seq(1, d), which(diag(p) != 0))
    while (sum(B) <= ne & niter < 2 * d) {
      ix = sample(ce, 1)
      B[ix] = TRUE
      if (dag) {
        g = graph_from_adjacency_matrix(B)
        B[ix] = is.dag(g)
      }
      niter = niter + 1
    }
    ## simulated edge indices
    ei = which(B != 0)
    B[ei] = runif(length(ei), min = coef.min, max = coef.max) * sample(c(-1, 1), length(ei), replace = T)
    ## check matrix (invertible)
    Det = abs(det(diag(p) - B))
    if(Det > 1e-2 || Matrix::norm(B, "f") < 1) {
      Matrix::Matrix(B, sparse = T)
    } else {
      NULL
    }
  }
  ER = function() {
    B = as_adjacency_matrix(sample_pa(p, power = sparse), names = FALSE, sparse = FALSE)
    ord = sample(p)
    B = B[ord, ord]
    ## simulated edge indices
    ei = which(B != 0)
    B[ei] = runif(length(ei), min = coef.min, max = coef.max) * sample(c(-1, 1), length(ei), replace = T)
    ## check matrix (invertible)
    Det = abs(det(diag(p) - B))
    if(Det > 1e-2 || Matrix::norm(B, "f") < 1) {
      Matrix::Matrix(B, sparse = T)
    } else {
      NULL
    }
  }
  ## eQTL mapping configuration
  ## simulate F2 family by qpgraph. Therefore, for n genes, we simulate candidate eQTls
  ## for simplicity, we put all eQtls on same chromosome, with even spacing.
  eQTLs = function(qlen = 100) {
    nq = k / p
    # for each genes, candidate eQtls
    cq = nq / sqtl
    # simulated genetic map
    map = sim.map(len = qlen * p, n.mar = cq * p, eq.spacing = TRUE, include.x = FALSE)
    # simulate f2 cross genotype
    f2  = sim.cross(map, type = "f2", n.ind = n)
    X = t(pull.geno(f2)) - 1
    F = matrix(0, nrow = p, ncol = cq * p)
    Sk = lapply(1:p, function(i) {
      s = seq((i-1) * cq, i * cq - 1) + 1
      ## v = sample(seq(1, cq - nq + 1), 1)
      if (!span) {
        F[i, s[1:nq]] <<- runif(nq, min = esize[1], max = esize[2])
      } else {
        F[i, s[round(seq(1, cq, length.out = nq))]] <<- runif(nq, min = esize[1], max = esize[2])
      }
      s
    })
    list(F = F, X = X, Sk = Sk, Cross = f2, map = map)
  }
  nQTLs = function(qlen = 100) {
    nq = k / p
    # for each genes, candidate eQtls
    cq = nq / sqtl
    # simulated genetic map
    map = sim.map(len = qlen * p, n.mar = cq * p, eq.spacing = TRUE, include.x = FALSE)
    # simulate f2 cross genotype
    f2  = sim.cross(map, type = "f2", n.ind = n)
    X = t(pull.geno(f2)) - 1
    F = matrix(0, nrow = p, ncol = cq * p)
    Sk = lapply(1:p, function(i) {
      s = seq((i-1) * nq, (i-1) * nq + cq - 1) + 1
      ## v = sample(seq(1, cq - nq + 1), 1)
      if (!span) {
        F[i, s[1:nq]] <<- runif(nq, min = esize[1], max = esize[2])
      } else {
        F[i, s[round(seq(1, cq, length.out = nq))]] <<- runif(nq, min = esize[1], max = esize[2])
      }
      s
    })
    list(F = F, X = X, Sk = Sk, Cross = f2, map = map)
  }
  sQTLs = function(qlen = 100) {
    nq = k / p
    # for each genes, candidate eQtls
    cq = nq / sqtl
    # simulated genetic map
    map = sim.map(len = qlen * p, n.mar = cq * p, eq.spacing = TRUE, include.x = FALSE)
    # simulate f2 cross genotype
    f2  = sim.cross(map, type = "f2", n.ind = n)
    X = t(pull.geno(f2)) - 1
    F = matrix(0, nrow = p, ncol = cq * p)
    Sk = lapply(1:p, function(i) {
      s = seq(i, i + cq - 1)
      ## v = sample(seq(1, cq - nq + 1), 1)
      if (!span) {
        F[i, s[1:nq]] <<- runif(nq, min = esize[1], max = esize[2])
      } else {
        F[i, s[round(seq(1, cq, length.out = nq))]] <<- runif(nq, min = esize[1], max = esize[2])
      }
      s
    })
    list(F = F, X = X, Sk = Sk, Cross = f2, map = map)
  }
  B = NULL
  while (is.null(B)) {
    if (type == "DG") {
      B = DG()
    } else {
      B = ER()
    }
  }
  QTL = if (overlap == "none") {
    eQTLs(qlen = qdist)
  } else if (overlap == "weak") {
    nQTLs(qlen = qdist)
  } else {
    sQTLs(qlen = qdist)
  }

  F = QTL$F
  if (noqtl) {
    rmix = sample(p, rbinom(1, p, rmv))
    F[rmix, ] = 0
  }
  X = QTL$X
  ## Replace the simulated X by Human Sapiens's Sample Genotype
  data(G)
  ni = sample(seq(1, length(G)), 1)
  Geno = G[[ni]]
  nq = k / p
  # for each genes, candidate eQtls
  cq = nq / sqtl
  nix = sample(seq(1, ncol(Geno[[1]])), ncol(X))
  if (nrow(Geno[[1]]) < cq) {
    stop("`Pool of SNP to be sampled < expected number of SNP, Only support up to 15 SNP per genes!`")
  }
  niy = sample(seq(1, nrow(Geno[[1]])), cq)
  CandGene = lapply(sample(Geno, p), function(x){x[niy, nix]})
  X = do.call(rbind, CandGene)
  cat(sprintf("DIM of QTLX = (%d, %d); DIM of GENO = (%d, %d)\n", nrow(QTL$X), ncol(QTL$X), nrow(X), ncol(X)))
  Sk = QTL$Sk
  E = sqrt(sigma2) * t(rmvnorm(n, mean = rep(0, p), sigma = diag(p)))
  U = tcrossprod(rnorm(p, 0, intercept), rep(1, n))
  Y = solve(diag(p) - B) %*% (F %*% X + U + E)
  rownames(X) = sapply(1:nrow(X), function(i){sprintf("rs_%05d", i)})
  rownames(Y) = sapply(1:nrow(Y), function(i){sprintf("g_%05d", i)})
  V = F %*% X
  SNR = sqrt(var(as.vector(V)) / var(as.vector(E)))
  Cross = QTL$Cross
  Cross$pheno = as.data.frame(t(Y))
  list(
    Data = list(X = X, Y = Y, Sk = Sk),
    Vars = list(B = B, F = Matrix(F, sparse = T), Mu = U, n = n, p = p, k = k, q = ncol(F), SNR = SNR),
    eQTL = list(Cross  = Cross, map = QTL$map),
    ij = list(ni = ni, nix = nix, niy = niy)
  )
}

##' generate simulated network structure for eQTLSEM algorithm, with fixed GRN and eQTL mapping
##' require("eQTLsem")
##' @title eQTLdataExt4H
##' @param B Gene regulatory network matrix
##' @param F eQTL mapping matrix
##' @param sigma2 noise of observations
##' @export
eQTLdataExt4H = function(B, F, Sk, QTL, n, k, sigma2 = 0.01, intercept = 5.0, ij) {
  p = nrow(B)
  f2  = sim.cross(QTL$map, type = "f2", n.ind = n)
  QX = t(pull.geno(f2)) - 1
  ####
  ni  = ij$ni
  nix = ij$nix
  niy = ij$niy
  ## Replace the simulated X by Human Sapiens's Sample Genotype
  data(G)
  Geno = G[[ni]]
  nix = sample(seq(1, ncol(Geno[[1]])), ncol(QX))
  CandGene = lapply(sample(Geno, p), function(x){x[niy, nix]})
  X = do.call(rbind, CandGene)
  cat(sprintf("DIM of QTLX = (%d, %d); DIM of GENO = (%d, %d)\n", nrow(QX), ncol(QX), nrow(X), ncol(X)))
  ###
  E = sqrt(sigma2) * t(rmvnorm(n, mean = rep(0, p), sigma = diag(p)))
  U = tcrossprod(rnorm(p, 0, intercept), rep(1, n))
  Y = solve(diag(p) - B) %*% (as.matrix(F) %*% X + U + E)
  rownames(X) = sapply(1:nrow(X), function(i){sprintf("rs_%05d", i)})
  rownames(Y) = sapply(1:nrow(Y), function(i){sprintf("g_%05d", i)})
  V = F %*% X
  SNR = sqrt(var(as.vector(V)) / var(as.vector(E)))
  Cross = f2
  Cross$pheno = as.data.frame(t(Y))
  list(
    Data = list(X = X, Y = Y, Sk = Sk),
    Vars = list(B = B, F = F, Mu = U, n = n, p = p, k = k, q = ncol(F), SNR = SNR),
    eQTL = list(Cross  = Cross, map = QTL$map),
    ij = list(ni = ni, nix = nix, niy = niy)
  )
}




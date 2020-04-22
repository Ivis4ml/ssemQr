## Objective value of regression
##' @title obj.ridgeRegression
##' @param X eQTL matrix
##' @param Y gene expression matrix
##' @param fit regression fit result object
##' @param trans  if rows for sample, trans = TRUE, otherwise, trans = FALSE. Default FALSE
##' @return error squared norm of ||(I-B)Y - FX||_2^2
obj.ridgeRegression = function(X, Y, fit, trans = F) {
  if (!trans) {
    X = t(X)
    Y = t(Y)
  }
  error = .Call("ObjL2Regression", X, Y, fit, PACKAGE = "ssemQr")
  error
}

## Lambda max of ridge regression
##' @title lamax.ridgeRegression
##' @param X eQTL matrix
##' @param Y gene expression matrix
##' @return max lambda
lamax.ridgeRegression = function(X, Y, Sk, n, p, k) {
  X = t(X)
  Y = t(Y)
  lambda = .Call("L2lamax", X, Y, Sk, n, p, k, PACKAGE = "ssemQr")
  lambda
}


inert_opt = function(opts = c("continuous", "linear"), init = 0) {
  switch(
    opts,
    "continuous" = function(k) {
      init
    },
    "linear"  = function(k) {
      (k - 1) / (k + 2)
    }
  )
}

##' @title proc.centerSSEM
##' @param X eQTL matrix
##' @param Y gene expression matrix
##' @return  centered X and Y and mean vectors
proc.centerSSEM = function(X, Y) {
  mX = rowMeans(X)
  mY = rowMeans(Y)
  center = function(m) {
    apply(m, 1, function(x){ x - mean(x) })
  }
  X = center(X)
  Y = center(Y)
  list(X = t(X),
       Y = t(Y),
       mX = mX, mY = mY)
}


sigma2SSEM = function(X, Y, B, F, n, p) {
  error = (norm(Y - B %*% Y - F %*% X, "f"))**2
  error / (n * p)
}

loglikSSEM = function(B, F, Wb, Wf, lambda, rho, sigma2, Det, n, p) {
  loglik = -n/2 * log(Det**2)
  if (is.infinite(lambda)) {
    l1B = 0 * Wb * abs(B)
  } else {
    l1B = lambda * Wb * abs(B)
  }
  if (is.infinite(rho)) {
    l1F = 0 * Wf * abs(F)
  } else {
    l1F = rho * Wf * abs(F)
  }
  l1B[is.infinite(Wb)] = 0
  l1F[is.infinite(Wf)] = 0
  loglik + p * n / 2 * log(2 * pi * sigma2) + n * p + sum(l1B) + sum(l1F)
}

loglikENET = function(B, F, Wb, Wf, lambda, alpha, rho, sigma2, Det, n, p) {
  loglik = -n/2 * log(Det**2)
  if (is.infinite(lambda)) {
    l1B = 0 * Wb * abs(B)
  } else {
    l1B = lambda * Wb * abs(B)
  }
  if (is.infinite(rho)) {
    l2F = 0 * (1 - alpha) * sum(F**2)
    l1F = 0 * alpha * Wf * abs(F)
  } else {
    l2F = rho * (1 - alpha) * sum(F**2)
    l1F = rho * alpha * Wf * abs(F)
  }
  l1B[is.infinite(Wb)] = 0
  l1F[is.infinite(Wf)] = 0
  loglik + p * n / 2 * log(2 * pi * sigma2) + n * p + sum(l1B) + l2F + sum(l1F)
}

##' function generator function for network B
##' @title cwiseGradient4SSEM2B
##' @param n number of observations
##' @param c cofactor vector
##' @param Y Matrix of gene expression
##' @param R Residual matrix
##' @param Y2norm Column of YtY
##' @param sigma2 noise variance
##' @return function whose argument is column vector bi
cwiseGradient4SSEM2B = function(n, c, Y, R, Y2norm, sigma2) {
  function(x) {
    n * t(c) + (Y2norm * x - tcrossprod(R, Y)) / sigma2
  }
}

## speed-up lipschitz computation
## abs(z - c %*% b) >= abs(z) - abs(c %*% b) >= abs(z) - sqrt(c2 * b2)
cwiseLipschitz4SSEM2B = function(n, z, c2, b2, Y2norm, sigma2, ImB, i) {
  if (abs(z) - sqrt(c2 * b2) > 0) {
    n * c2 / (abs(z) - sqrt(c2 * b2))**2 + Y2norm / sigma2
  } else {
    cwiseLipschitzSML(n, c2, Y2norm, sigma2, ImB, i)[1]
    ## Inf
  }
}

cwiseLipschitzSML = function(n, c2, Y2norm, sigma2, ImB, i) {
  ## x^tWx + 2R^tx + C
  p = nrow(ImB)
  I = chol2inv(chol(crossprod(ImB[, -i])))
  D = det(crossprod(ImB[, -i]))
  W = diag(p - 1) - ImB[-i, -i] %*% tcrossprod(I, ImB[-i, -i])
  R = ImB[-i, -i] %*% tcrossprod(I, ImB[i, -i, drop = FALSE])
  C = 1 - ImB[i, -i, drop = FALSE] %*% tcrossprod(I, ImB[i, -i, drop = FALSE])
  e = 1e-6
  v = solve(W + diag(e, p - 1), -1 * R)
  L = n * c2 / (crossprod(v, I %*% v) + 2 * crossprod(R, v) + C) / (D + e) + Y2norm / sigma2
  if (L < 0) {
    L = Inf
  }
  L
}

## proximal operator of lasso
## prox_{c} c/2||x - u||_2^2 + lambda * w ||x||_1
proxLasso = function(u, w, lambda, c) {
  pmax(u - lambda * w / c, 0) + pmin(u + lambda * w / c, 0)
}

## proximal operator of elastic net
## proc_{c} c/2||x - u||_2^2 + (1 - alpha) * lambda * ||x||_2^2 + alpha * lambda * w ||x||_1
proxElasticNet = function(u, w, alpha, lambda, c) {
  (pmax(u - alpha * lambda * w / c, 0) + pmin(u + alpha * lambda * w / c, 0)) / (1 + 2 * (1 - alpha) * lambda / c)
}

##' function generator function for network B
##' @title rwiseGradient4SSEM2F
##' @param ImBi updated (I - B)[i,]
##' @param Y Matrix of gene expression
##' @param Xi Submatrix of eQTLs
##' @param sigma2 noise variance
##' @return function whose argument is column vector fi
rwiseGradient4SSEM2F = function(ImBi, Y, Xi, sigma2) {
  function(f) {
    -1 / sigma2 * tcrossprod(ImBi %*% Y - f %*% Xi, Xi)    ## 1 x sk
  }
}

bayesianInfocriterion = function(X, Y, B, F, mu, Det, sigma2, p, c) {
  n    = ncol(Y)
  q    = nrow(X)
  err  = norm(Y - B %*% Y - F %*% X - tcrossprod(mu, rep(1, n)), "f")**2
  logl = -n/2*log(Det**2) + n * p / 2 * log(sigma2) + err / 2 / sigma2
  df   = sum(B != 0) + sum(F != 0) + 1
  if (q > n & p > n) { # high dimensional case
    kappa = 0
    2 * logl + 6 * (1 + kappa) * df * log(n)  ## CL-BIC
  } else {
    2 * logl + df * log(n)
  }
}

threshSEM = function(n, p, q) {
  sqrt(log(p * (p - 1) + q) / n / p)
}


QiPALM = function(X, Y, fit) {
  centered = proc.centerSSEM(X, Y)
  X  = centered[["X"]]            # k x n
  Y  = centered[["Y"]]            # p x n
  p  = nrow(Y)
  n  = ncol(Y)
  abs(n * solve(diag(p) - fit$B) + (tcrossprod((diag(p) - fit$B) %*% Y - fit$F %*% X, Y) / fit$sigma2))
}


FPR = function(X, B, PREC = 0) {
  X = as.matrix(X)
  B = as.matrix(B)
  sum(abs(X[B == 0] > PREC)) / sum(B == 0)
}


##' get TPR FDR of B and F
##' @title getMetrics
##' @param data modeled data, contain real values
##' @export
getMetrics = function(data, fit, only.Qtl = FALSE) {
  threshold = 0
  if (only.Qtl) {
    TPRF = fssemR:::TPR(fit$F, data$Vars$F, threshold)
    FDRF = fssemR:::FDR(fit$F, data$Vars$F, threshold)
    list(F = list(TPR = TPRF, FDR = FDRF))
  } else {
    TPRB = fssemR:::TPR(fit$B, data$Vars$B, threshold)
    FDRB = fssemR:::FDR(fit$B, data$Vars$B, threshold)
    TPRF = fssemR:::TPR(fit$F, data$Vars$F, threshold)
    FDRF = fssemR:::FDR(fit$F, data$Vars$F, threshold)
    list(B = list(TPR = TPRB, FDR = FDRB),
         F = list(TPR = TPRF, FDR = FDRF))
  }
}

##' calculated precision-recall curve with specified quantile cutoff
##' @title calcPR
##' @param X  estimated matrix of GRN or eQTL-mapping
##' @param B  ground true value in simulation
##' @param by precision of curve
##' @return list of PRcurve
##' @export
calcPR = function(X, B, by = 5 * 1e-4) {
  X = as.vector(abs(X))
  B = as.vector(B != 0)
  probs = seq(0, 1, by = by)
  TB = sum(B)
  cutoff = c(Inf, sort(quantile(c(max(X), 0), probs = probs), decreasing = TRUE))
  PR = plyr::llply(cutoff, function(c) {
    c(sum(B[X >= c]) / TB, sum(B[X >= c]) / sum(X >= c))
  })
  PR = do.call(rbind, PR)
  colnames(PR) = c("recall", "precision")
  as.data.frame(PR)
}


##' Download GTeX borrow from Yarn
##' @title downloadGTEx
##' @return ExpressionSet set
##' @export
##' @examples
##' # obj <- downloadGTExv6p(type = 'genes', file = '/media/xinchou/Storage/gtexv6p.rds')
downloadGTExv6p = function (type = "genes", file = NULL, ...)
{
  phenoFile <- "https://storage.googleapis.com/gtex_analysis_v6p/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
  pheno2File <- "https://storage.googleapis.com/gtex_analysis_v6p/annotations/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"
  geneFile <- "https://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz"
  message("Downloading and reading files")
  pdFile <- tempfile("phenodat", fileext = ".txt")
  downloader::download(phenoFile, destfile = pdFile)
  pd <- readr::read_tsv(pdFile)
  pd <- as.matrix(pd)
  rownames(pd) <- pd[, "SAMPID"]
  ids <- sapply(strsplit(pd[, "SAMPID"], "-"), function(i) paste(i[1:2],
                                                                 collapse = "-"))
  pd2File <- tempfile("phenodat2", fileext = ".txt")
  downloader::download(pheno2File, destfile = pd2File)
  pd2 <- readr::read_tsv(pd2File)
  pd2 <- as.matrix(pd2)
  rownames(pd2) <- pd2[, "SUBJID"]
  pd2 <- pd2[which(rownames(pd2) %in% unique(ids)), ]
  pd2 <- pd2[match(ids, rownames(pd2)), ]
  rownames(pd2) <- colnames(counts)
  pdfinal <- Biobase::AnnotatedDataFrame(data.frame(cbind(pd, pd2)))
  if (type == "genes") {
    countsFile <- tempfile("counts", fileext = ".gz")
    downloader::download(geneFile, destfile = countsFile)
    cnts <- suppressWarnings(readr::read_tsv(geneFile, skip = 2))
    genes <- unlist(cnts[, 1])
    geneNames <- unlist(cnts[, 2])
    counts <- cnts[, -c(1:2)]
    counts <- as.matrix(counts)
    rownames(counts) <- genes
    for (i in 1:nrow(readr::problems(cnts))) {
      counts[readr::problems(cnts)$row[i], readr::problems(cnts)$col[i]] <- 1e+05
    }
    throwAway <- which(rowSums(counts) == 0)
    counts <- counts[-throwAway, ]
    genes <- sub("\\..*", "", rownames(counts))
    host <- "dec2013.archive.ensembl.org"
    biomart <- "ENSEMBL_MART_ENSEMBL"
    dataset <- "hsapiens_gene_ensembl"
    attributes <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name",
                    "start_position", "end_position", "gene_biotype")
  }
  message("Creating ExpressionSet")
  pdfinal <- pdfinal[match(colnames(counts), rownames(pdfinal)),
                     ]
  es <- Biobase::ExpressionSet(as.matrix(counts))
  Biobase::phenoData(es) <- pdfinal
  pData(es)["GTEX-YF7O-2326-101833-SM-5CVN9", "SMTS"] <- "Skin"
  pData(es)["GTEX-YEC3-1426-101806-SM-5PNXX", "SMTS"] <- "Stomach"
  #message("Annotating from biomaRt")
  #es <- yarn::annotateFromBiomart(obj = es, genes = genes, host = host,
  #                          biomart = biomart, dataset = dataset, attributes = attributes)
  message("Cleaning up files")
  unlink(pdFile)
  unlink(pd2File)
  unlink(countsFile)
  if (!is.null(file))
    saveRDS(es, file = file)
  return(es)
}



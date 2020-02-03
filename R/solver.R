##' @title ridgeRegression
##' @description Ridge regression for initialization
##' @param X     eQTL matrix.
##' @param Y     Gene expression matrix
##' @param Sk    eQTL index of genes
##' @param gamma  Hyperparameter for ridge regression
##' @param n      number of observations
##' @param p      number of genes
##' @param k      number of eQTLs
##' @param trans  if rows for sample, trans = TRUE, otherwise, trans = FALSE. Default FALSE
##' @param sparse Default is False
##' @return fit List of SEM model
##' \describe{
##' \item{Bs}{ coefficient matrices of gene regulatory networks}
##' \item{fs}{ eQTL's coefficients w.r.t each gene}
##' \item{Fs}{ coefficient matrices of eQTL-gene effect}
##' \item{mu}{ Bias vector}
##' \item{sigma2}{ estimate of covariance in SEM}
##' }
##' @export
##' @importFrom stats rbinom rnorm runif sd
##' @importFrom Matrix Matrix
ridgeRegression = function(X, Y, Sk, gamma, n, p, k, trans = FALSE, sparse = FALSE) {
  if (!trans) {
    X = t(X)
    Y = t(Y)
  }
  fit = .Call("L2Regression", X, Y, Sk, gamma, n, p, k, PACKAGE = "SEMQTLR")
  if (sparse) {
    fit$F = Matrix(fit$F, sparse = T)
  }
  fit
}


## cross-validation on ridge regression
##' @title cv.ridgeRegression
##' @param X      eQTL matrix
##' @param Y      Gene expression matrix
##' @param Sk      eQTL index of genes
##' @param ngamma  number of hyper-parameter in CV
##' @param nfold   CVfold number. Default 5/10
##' @param n       number of observations
##' @param p       number of genes
##' @param k       number of eQTLs
##' @return gamma_min optimal gamma to minimize cross-validation error
##' @export
##' @examples
##' gamma = cv.ridgeRegression(data$Data$X, data$Data$Y, data$Data$Sk, ngamma = 20, nfold = 5, 100, 10, 60)
cv.ridgeRegression = function(X, Y, Sk, ngamma = 20, nfold = 5, n, p, k) {
  cverr.mat = matrix(0, nrow = ngamma, ncol = nfold)
  cvfold    = sample(seq(1, nfold), size = n, replace = T)
  gamma_max = lamax.ridgeRegression(X, Y, Sk, n, p, k)
  gamma_factors = 10 ** (seq(0, -6, length.out = ngamma)) * gamma_max
  Xtrain = vector("list", nfold)
  Xtest  = vector("list", nfold)
  Ytrain = vector("list", nfold)
  Ytest  = vector("list", nfold)
  for (i in 1:nfold) {
    Xtrain[[i]] = X[, cvfold != i, drop = F]
    Xtest[[i]]  = X[, cvfold == i, drop = F]
    Ytrain[[i]] = Y[, cvfold != i, drop = F]
    Ytest[[i]]  = Y[, cvfold == i, drop = F]
  }
  igamma = 1
  for (gamma in gamma_factors) {
    for (i in 1:nfold) {
      ntrain = sum(cvfold != i)
      ntest  = sum(cvfold == i)
      fit = ridgeRegression(Xtrain[[i]], Ytrain[[i]], Sk, gamma, ntrain, p, k, trans = F)
      cverr.mat[igamma, i] = obj.ridgeRegression(Xtest[[i]], Ytest[[i]], fit) / ntest
    }
    igamma = igamma + 1
  }
  cvmean = rowMeans(cverr.mat)
  print(cvmean)
  gamma.min = gamma_factors[which.min(cvmean)]
  gamma.min
}

##' @title initRegression
##' @param X      eQTL matrix
##' @param Y      Gene expression matrix
##' @param Sk     eQTL index of genes
##' @param n      number of observations
##' @param p      number of genes
##' @param k      number of eQTLs
##' @export
initRegression = function(X, Y, Sk, n, p, k) {
  centered = proc.centerSSEM(X, Y)
  X  = centered[["X"]]
  Y  = centered[["Y"]]
  mX = centered[["mX"]]
  mY = centered[["mY"]]
  B = matrix(0, nrow = p, ncol = p)
  F = matrix(0, nrow = p, ncol = k)
  sigma2 = sigma2SSEM(X, Y, B, F, n, p)
  Wb = matrix(1, nrow = p, ncol = p)
  diag(Wb) = 0
  Wb = 1 / Wb
  Wf = matrix(0, nrow = p, ncol = k)
  lapply(1:p, function(i) {
    Wf[i, Sk[[i]]] <<- 1
  })
  Wf = 1 / Wf
  fit = list(B = B, F = t(F), sigma2 = sigma2, Wb = Wb, Wf = t(Wf))
}

##' Sparse Structural Equation Model (SEM) with iPALM for eQTL mapping and GRN inference
##' @title SSEMiPALM
##' @description Implementing iPALM algorithm to solve the sparse SEM for eQTL mapping and network inference
##' @param X  eQTL matrix
##' @param Y  Gene expression matrix
##' @param B  initialized GRN-matrix
##' @param F  initialized eQTL effect matrix
##' @param Sk  eQTL index of genes
##' @param sigma2 initialized noise variance from ridge regression
##' @param lambda Hyperparameter of lasso term in SSEM for B
##' @param rho    Hyperparameter of lasso term in SSEM for F
##' @param Wb weight matrix of B for adaptive lasso terms. Default as 1 / B(init) from eidge regression.
##' @param Wf weight matrix of F for adaptive lasso terms. Default as 1 / F(init) from ridge regression.
##' @param p  number of genes
##' @param maxit maximum iteration number. Default 100
##' @param inert inertial function for iPALM. Default as k-1/k+2
##' @param threshold convergence threshold. Default 1e-6
##' @param verbose Default TRUE
##' @param sparse Sparse Matrix or not
##' @param strict Converge strictly or not. Default False
##' @param trans  Fs matrix is transposed to k x p or not. If Fs from ridge regression, trans = TRUE, else, trans = FALSE
##' @param B2norm B2norm matrices upper bound generated from ridge regression. Default NULL.
##' @export
##' @importFrom Matrix norm
##' @importFrom rARPACK eigs_sym
SSEMiPALM = function(X, Y, B, F, Sk, sigma2, lambda, rho, Wb = 1 / abs(B), Wf = 1 / abs(F), p,
                     maxit = 100, inert = inert_opt("linear"), threshold = 1e-6,
                     verbose = TRUE, sparse = TRUE, trans = FALSE, B2norm = NULL,
                     strict = FALSE) {
  Xr = X
  Yr = Y
  centered = proc.centerSSEM(X, Y)
  X  = centered[["X"]]
  Y  = centered[["Y"]]
  mX = centered[["mX"]]
  mY = centered[["mY"]]
  q  = nrow(X)
  n  = ncol(Y)
  if (trans) {
    F  = t(F)
    Wf = t(Wf)
  }
  if(is.null(B2norm)) {
    B2norm = colSums(B**2)
    B = matrix(0, nrow = p, ncol = p)
  }
  Fx = F %*% X
  Y2norm = vector("numeric", p)
  X2eig  = vector("numeric", p)
  for (i in 1:p) {
    yi = Y[i, , drop = FALSE]
    xi = X[Sk[[i]], , drop = FALSE]
    Y2norm[i] = tcrossprod(yi)
    X2eig[i]  = max(eigen(tcrossprod(xi))$values)
  }

  ## block coordinate descent of SSEM (Col-major on B and Row-major on F)
  niter  = 1
  ImB    = diag(p) - B
  Det    = det(ImB)
  invImB = solve(ImB)
  sigma2 = sigma2SSEM(X, Y, B, F, n, p)
  LogLik = loglikSSEM(B, F, Wb, Wf, lambda, rho, sigma2, Det, n, p)
  ## for wolfe update
  gamma0  = 1
  eta0    = 1
  precB   = 0
  precF   = 0
  ninc = TRUE
  Bhist = list(B, B)
  Fhist = list(F, F)
  while (niter <= maxit) {
    tau = inert(niter)
    Bm  = Bhist[[2]] + tau * (Bhist[[2]] - Bhist[[1]])
    Lhist  = LogLik
    if (!ninc) {
      if (gamma0 <= 1e6) {
        gamma0 = min(gamma0 * 1.01, 1e6)
      }
    }
    for (i in 1:p) {
      ci = invImB[i, , drop = FALSE]
      bi = Bm[, i, drop = FALSE]
      Ri = Y - B[,-i] %*% Y[-i,] - Fx
      gradient = cwiseGradient4SSEM2B(n, ci, Y[i, ,drop = FALSE], Ri, Y2norm[i], sigma2)
      gi = gradient(bi)[-i, , drop = FALSE]
      z  = ci[,i] * Det                      ## cofactor of (i,i) in I-B
      c2 = sum((ci[,-i] * Det)**2)           ## cofactor of (-i, i) in I-B
      b2 = B2norm[i]
      li = cwiseLipschitz4SSEM2B(n, z, c2, b2, Y2norm[i], sigma2, ImB, i)
      if (is.infinite(li)) {
        next
      }
      ## Armoji scheme
      if (niter == 1) {
        li = (1 + 2 * tau) / (1 - 2 * tau) * li
      } else {
        li = (1 + 2 * tau) / 2 / (1 - tau) * li
      }
      gamma = gamma0
      while (TRUE) {
        ui = bi[-i,] - 1 / (gamma * li) * gi
        wi = Wb[-i, i]
        xi = proxLasso(ui, wi, lambda, gamma * li)
        Det_update = invImB[i, i] - tcrossprod(xi, invImB[i, -i])[1]
        if (Det_update != 0) {
          break
        }
        gamma = gamma * 1.01
      }
      ## Covariates update
      B[-i, i] = xi
      ImB = diag(p) - B
      Det = tcrossprod(ImB[,i], invImB[i, ,drop = FALSE])[1] * Det
      dB  = Bhist[[2]][, i, drop = FALSE] - B[, i, drop = FALSE]
      ## Woodbury matrix identity
      invImB = invImB - invImB %*% dB %*% invImB[i, , drop = FALSE] / (1 + invImB[i, , drop = FALSE] %*% dB)[1]
    }
    ## update block F
    eta = (1 + sqrt(1 + 4 * eta0**2)) / 2
    Fm  = Fhist[[2]] + (eta0 - 1) / eta * (Fhist[[2]] - Fhist[[1]])    # inertial term for F block
    eta0 = eta
    for (i in 1:p) {
      fi = Fm[i, Sk[[i]]]
      gradient = rwiseGradient4SSEM2F(ImB[i, , drop = FALSE], Y, X[Sk[[i]], , drop = FALSE], sigma2)
      gi = gradient(fi)
      li = X2eig[i] / sigma2
      ui = fi - 1 / li * gi
      wi = Wf[i, Sk[[i]]]
      xi = proxLasso(ui, wi, rho, li)
      F[i, Sk[[i]]] = xi
    }
    ## covariates update
    Fx = F %*% X
    sigma2 = sigma2SSEM(X, Y, B, F, n, p)
    ## convergenece
    Berr = norm(B - Bhist[[2]], "f") / (1 + norm(Bhist[[2]], "f"))
    Ferr = norm(F - Fhist[[2]], "f") / (1 + norm(Fhist[[2]], "f"))
    LogLik = loglikSSEM(B, F, Wb, Wf, lambda, rho, sigma2, Det, n, p)
    Lerr = abs(Lhist - LogLik) / (1 + abs(Lhist))
    ninc  = (Lhist - LogLik) < 0
    converged = if (strict) {
      Berr + Ferr <= threshold && Lerr <= threshold
    } else {
      Lerr <= threshold
    }
    if (verbose) {
      cat(sprintf("SSEM:\tniter = %d,\trelerr = %f,\tprevLik = %f, \tlogLik = %f,\tσ2 = %f\n", niter, Berr + Ferr, Lhist, LogLik, sigma2))
    }
    niter = niter + 1
    Bhist = list(Bhist[[2]], B)
    Fhist = list(Fhist[[2]], F)
    if (converged || niter > maxit || sigma2 < 1e-6 || (!ninc & gamma0 >= 1e6)) {
      precB = min(abs(B - Bhist[[1]]))
      precF = min(abs(F - Fhist[[1]]))
      break
    }
  } # blockwise coordinate descent
  mu = (diag(p) - B) %*% mY - sapply(1:p, function(i) {
    mX[Sk[[i]]] %*% F[i, Sk[[i]]]
  })
  sigma2 = sigma2SSEM(X, Y, B, F, n, p)
  if (sparse) {
    B = Matrix(B, sparse = T)
    F = Matrix(F, sparse = T)
  }
  list(B = B, F = F, mu = mu, sigma2 = sigma2, Det = Det, trans = TRUE)
}

##' Initialized Lambda coefficients (lambda, rho) for Sparse SEM
##' @title initLambdaSSEM
##' @param X  eQTL matrix
##' @param Y  Gene expression matrix
##' @param Sk eQTL index of genes
##' @param Wb weight matrix of B for adaptive lasso terms. Default as 1 / B(init) from eidge regression.
##' @param Wf weight matrix of F for adaptive lasso terms. Default as 1 / F(init) from ridge regression.
##' @param p number of genes
##' @param q number of eQTLs
initLambdaSSEM2F = function(X, Y, Sk, Wb, Wf, p, q) {
  centered = proc.centerSSEM(X, Y)
  X  = centered[["X"]]            # k x n
  Y  = centered[["Y"]]            # p x n
  mX = centered[["mX"]]
  mY = centered[["mY"]]
  n  = ncol(Y)
  B0 = matrix(0, nrow = p, ncol = p)
  F0 = matrix(0, nrow = p, ncol = q)
  sigma2 = sigma2SSEM(X, Y, B0, F0, n, p)
  # \nabla_B = \frac{1}{\sigma^2}|YY^T - 0XY^T| <= \lambda * Wb    B p x p
  # \nabla_F = \frac{1}{\sigma^2}|YX^T - 0XX^T| <= \rho * Wf       F p x k
  gF = abs(1 / sigma2 * tcrossprod(Y, X)) / Wf
  rho  = 0
  for(i in 1:p) {
    rho = max(rho, max(gF[i,Sk[[i]]]))
  }
  rho
}

##' Initialized Lambda coefficients (lambda, rho) for Sparse SEM
##' @title initLambdaSSEM
##' @param X  eQTL matrix
##' @param Y  Gene expression matrix
##' @param fit SEM fit Object as lambda = Inf
##' @param Wb weight matrix of B for adaptive lasso terms. Default as 1 / B(init) from eidge regression.
initLambdaSSEM2B = function(X, Y, fit, Wb) {
  centered = proc.centerSSEM(X, Y)
  X  = centered[["X"]]            # k x n
  Y  = centered[["Y"]]            # p x n
  max(abs(1 / fit$sigma2 * (tcrossprod(Y) - fit$F %*% tcrossprod(X, Y)) / Wb))
}

##' Optimized Sparse SEM by BIC
##' @title opt.SSEMiPALM
##' @param X  eQTL matrix
##' @param Y  Gene expression matrix
##' @param B  initialized GRN-matrix
##' @param F  initialized eQTL effect matrix
##' @param Sk eQTL index of genes
##' @param sigma2 initialized noise variance from ridge regression
##' @param Wb weight matrix of B for adaptive lasso terms. Default as 1 / B(init) from eidge regression.
##' @param Wf weight matrix of F for adaptive lasso terms. Default as 1 / F(init) from ridge regression.
##' @param nlambda  number of lambda
##' @param nrho number of rho
##' @param p  number of genes
##' @param wt adaptive lasso or not. Default as TRUE
##' @export
opt.SSEMiPALM = function(X, Y, B, F, Sk, sigma2, Wb = NULL, Wf = NULL, nlambda = 20, nrho = 20, p, wt = TRUE) {
  B2norm = colSums(B**2)
  q = nrow(X)                     # number of eQTLs
  c = length(unlist(Sk))          # number of candidate eQTL
  if (wt) {
    Wb = 1 / abs(B)
    Wf = 1 / abs(F)
  }
  rho_max = initLambdaSSEM2F(X, Y, Sk, Wb, t(Wf), p, q)
  rho_factors = 10 ** seq(0, -3, length.out = nrho) * rho_max
  params = vector("list", nlambda * nrho)
  fit = list()
  j   = 1
  for (irho in 1:nrho) {
    rho = rho_factors[irho]
    for (ilambda in 1:nlambda) {
      if (ilambda == 1) {
        fit[[j]] = SSEMiPALM(X = X, Y = Y, B = B, F = F, Sk = Sk, sigma2 = sigma2,
                         lambda = Inf, rho = rho, Wb = Wb, Wf = Wf, p = p, maxit = 200,
                         threshold = 1e-4, trans = TRUE, strict = FALSE, verbose = FALSE)
        lambda_max = initLambdaSSEM2B(X, Y, fit[[j]], Wb)
        lambda_factors = 10 ** seq(0, -3, length.out = nlambda) * lambda_max
        cat(sprintf("SSEM@lambda = %4f, rho = %4f\n", lambda_factors[ilambda], rho))

      } else {
        fit[[j]] = SSEMiPALM(X = X, Y = Y, B = B, F = F, Sk = Sk, sigma2 = sigma2,
                         lambda = lambda_factors[ilambda], rho = rho, Wb = Wb, Wf = Wf,
                         p = p, maxit = 200, threshold = 1e-4, trans = TRUE, strict = FALSE, verbose = FALSE)
        cat(sprintf("SSEM@lambda = %4f, rho = %4f\n", lambda_factors[ilambda], rho))
      }
      err = bayesianInfocriterion(X, Y, fit[[j]]$B, fit[[j]]$F, fit[[j]]$mu, fit[[j]]$Det, fit[[j]]$sigma2, p, c)
      params[[j]] = c(lambda_factors[ilambda], rho, err)
      j = j + 1
    }
  }
  BICmat = data.frame(do.call(rbind, params))
  colnames(BICmat) = c("lambda", "rho", "BIC")
  BICmin = which.min(BICmat$BIC)
  list(lambda = BICmat[BICmin, 1], rho = BICmat[BICmin, 2],
       fit = fit[[BICmin]], minBIC = min(BICmat$BIC), BICs = BICmat)
}


##' Elstic-net based Sparse Structural Equation Model (SEM) with iPALM for eQTL mapping and GRN inference
##' @title SSEMeNet
##' @description Implementing iPALM algorithm to solve the Elstic-net based sparse SEM for eQTL mapping and network inference
##' @param X  eQTL matrix
##' @param Y  Gene expression matrix
##' @param B  initialized GRN-matrix
##' @param F  initialized eQTL effect matrix
##' @param Sk  eQTL index of genes
##' @param sigma2 initialized noise variance from ridge regression
##' @param lambda Hyperparameter of lasso term in SSEM for B
##' @param alpha  Hyperparameter for adaptive elastic net term in SSEM for F. Default (1 - alpha) x rho x ||F||_2^2 + alpha x rho x Wf x ||F||_1^1
##' @param rho    Hyperparameter of adaptive elastic net lasso term in SSEM for F
##' @param Wb weight matrix of B for adaptive lasso terms. Default as 1 / B(init) from eidge regression.
##' @param Wf weight matrix of F for adaptive lasso terms. Default as 1 / F(init) from ridge regression.
##' @param p  number of genes
##' @param maxit maximum iteration number. Default 100
##' @param inert inertial function for iPALM. Default as k-1/k+2
##' @param threshold convergence threshold. Default 1e-6
##' @param verbose Default TRUE
##' @param sparse Sparse Matrix or not
##' @param strict Converge strictly or not. Default False
##' @param trans  Fs matrix is transposed to k x p or not. If Fs from ridge regression, trans = TRUE, else, trans = FALSE
##' @param B2norm B2norm matrices upper bound generated from ridge regression. Default NULL.
##' @export
SSEMeNet = function(X, Y, B, F, Sk, sigma2, lambda, alpha, rho, Wb = 1 / abs(B), Wf = 1 / abs(F), p,
                    maxit = 100, inert = inert_opt("linear"), threshold = 1e-6,
                    verbose = TRUE, sparse = TRUE, trans = FALSE, B2norm = NULL,
                    strict = FALSE) {
  Xr = X
  Yr = Y
  centered = proc.centerSSEM(X, Y)
  X  = centered[["X"]]
  Y  = centered[["Y"]]
  mX = centered[["mX"]]
  mY = centered[["mY"]]
  q  = nrow(X)
  n  = ncol(Y)
  if (trans) {
    F  = t(F)
    Wf = t(Wf)
  }
  if (strict) {
    maxit = min(2000, max(5 * p, maxit))
  }

  if(is.null(B2norm)) {
    B2norm = colSums(B**2)
    B = matrix(0, nrow = p, ncol = p)
  }
  Fx = F %*% X
  Y2norm = vector("numeric", p)
  X2eig  = vector("numeric", p)
  for (i in 1:p) {
    yi = Y[i, , drop = FALSE]
    xi = X[Sk[[i]], , drop = FALSE]
    Y2norm[i] = tcrossprod(yi)
    X2eig[i]  = max(eigen(tcrossprod(xi))$values)
  }

  ## block coordinate descent of SSEM (Col-major on B and Row-major on F)
  niter  = 1
  ImB    = diag(p) - B
  Det    = det(ImB)
  invImB = solve(ImB)
  sigma2 = sigma2SSEM(X, Y, B, F, n, p)
  LogLik = loglikENET(B, F, Wb, Wf, lambda, alpha, rho, sigma2, Det, n, p)
  ## for wolfe update
  gamma0  = 1
  eta0    = 1
  ninc = TRUE
  Bhist = list(B, B)
  Fhist = list(F, F)
  while (niter <= maxit) {
    tau = inert(niter)
    Bm  = Bhist[[2]] + tau * (Bhist[[2]] - Bhist[[1]])
    Lhist  = LogLik
    if (!ninc) {
      if (gamma0 < 1e6) {
        gamma0 = min(gamma0 * 1.01, 1e6)
      }
    }
    for (i in 1:p) {
      ci = invImB[i, , drop = FALSE]
      bi = Bm[, i, drop = FALSE]
      Ri = Y - B[,-i] %*% Y[-i,] - Fx
      gradient = cwiseGradient4SSEM2B(n, ci, Y[i, ,drop = FALSE], Ri, Y2norm[i], sigma2)
      gi = gradient(bi)[-i, , drop = FALSE]
      z  = ci[,i] * Det                      ## cofactor of (i,i) in I-B
      c2 = sum((ci[,-i] * Det)**2)           ## cofactor of (-i, i) in I-B
      b2 = B2norm[i]
      li = cwiseLipschitz4SSEM2B(n, z, c2, b2, Y2norm[i], sigma2, ImB, i)
      if (is.infinite(li)) {
        next
      }
      ## Armoji scheme
      if (niter == 1) {
        li = (1 + 2 * tau) / (1 - 2 * tau) * li
      } else {
        li = (1 + 2 * tau) / 2 / (1 - tau) * li
      }
      gamma = gamma0
      while (TRUE) {
        ui = bi[-i,] - 1 / (gamma * li) * gi
        wi = Wb[-i, i]
        xi = proxLasso(ui, wi, lambda, gamma * li)
        Det_update = invImB[i, i] - tcrossprod(xi, invImB[i, -i])[1]
        if (Det_update != 0) {
          break
        }
        gamma = gamma * 1.02
      }
      ## Covariates update
      B[-i, i] = xi
      ImB = diag(p) - B
      Det = tcrossprod(ImB[,i], invImB[i, ,drop = FALSE])[1] * Det
      dB  = Bhist[[2]][, i, drop = FALSE] - B[, i, drop = FALSE]
      ## Woodbury matrix identity
      invImB = invImB - invImB %*% dB %*% invImB[i, , drop = FALSE] / (1 + invImB[i, , drop = FALSE] %*% dB)[1]
    }
    ## update block F
    eta = (1 + sqrt(1 + 4 * eta0**2)) / 2
    Fm  = Fhist[[2]] + (eta0 - 1) / eta * (Fhist[[2]] - Fhist[[1]])    # inertial term for F block
    eta0 = eta
    for (i in 1:p) {
      fi = Fm[i, Sk[[i]]]
      gradient = rwiseGradient4SSEM2F(ImB[i, , drop = FALSE], Y, X[Sk[[i]], , drop = FALSE], sigma2)
      gi = gradient(fi)
      li = X2eig[i] / sigma2
      ui = fi - 1 / li * gi
      wi = Wf[i, Sk[[i]]]
      xi = proxElasticNet(ui, wi, alpha, rho, li)
      F[i, Sk[[i]]] = xi
    }
    ## covariates update
    Fx = F %*% X
    sigma2 = sigma2SSEM(X, Y, B, F, n, p)
    ## convergenece
    Berr = norm(B - Bhist[[2]], "f") / (1 + norm(Bhist[[2]], "f"))
    Ferr = norm(F - Fhist[[2]], "f") / (1 + norm(Fhist[[2]], "f"))
    LogLik = loglikENET(B, F, Wb, Wf, lambda, alpha, rho, sigma2, Det, n, p)
    Lerr = abs(Lhist - LogLik) / (1 + abs(Lhist))
    ninc  = Lhist > LogLik
    converged = if (strict) {
      Berr + Ferr <= threshold && Lerr <= threshold
    } else {
      Lerr <= threshold
    }
    if (verbose) {
      cat(sprintf("SSEM:\tniter = %d,\trelerr = %f,\tprevLik = %f, \tlogLik = %f,\tσ2 = %f\n", niter, Berr + Ferr, Lhist, LogLik, sigma2))
    }
    niter = niter + 1
    Bhist = list(Bhist[[2]], B)
    Fhist = list(Fhist[[2]], F)
    if (converged || niter > maxit || sigma2 < 1e-6 || (!ninc & gamma0 >= 1e6)) {
      break
    }
  } # blockwise coordinate descent
  mu = (diag(p) - B) %*% mY - sapply(1:p, function(i) {
    mX[Sk[[i]]] %*% F[i, Sk[[i]]]
  })
  sigma2 = sigma2SSEM(X, Y, B, F, n, p)
  if (sparse) {
    B = Matrix(B, sparse = T)
    F = Matrix(F, sparse = T)
  }
  list(B = B, F = F, mu = mu, sigma2 = sigma2, Det = Det, trans = TRUE)
}


##' Optimized Elastic Net based Sparse SEM by BIC
##' @title opt.SSEMeNet
##' @param X  eQTL matrix
##' @param Y  Gene expression matrix
##' @param B  initialized GRN-matrix
##' @param F  initialized eQTL effect matrix
##' @param Sk eQTL index of genes
##' @param sigma2 initialized noise variance from ridge regression
##' @param Wb weight matrix of B for adaptive lasso terms. Default as 1 / B(init) from eidge regression.
##' @param Wf weight matrix of F for adaptive lasso terms. Default as 1 / F(init) from ridge regression.
##' @param nlambda number of lambda
##' @param nalpha number of alpha. Default [0, 1]
##' @param nrho number of rho
##' @param p  number of genes
##' @param wt adaptive lasso or not. Default as TRUE
##' @export
opt.SSEMeNet = function(X, Y, B, F, Sk, sigma2, Wb = NULL, Wf = NULL, nlambda = 20, nalpha = 5, nrho = 20, p, wt = TRUE) {
  B2norm = colSums(B**2)
  q = nrow(X)                     # number of eQTLs
  c = length(unlist(Sk))          # number of candidate eQTL
  if (wt) {
    Wb = 1 / abs(B)
    Wf = 1 / abs(F)
  }
  rho_max = initLambdaSSEM2F(X, Y, Sk, Wb, t(Wf), p, q)
  # rho_max = alpha * rho_max
  rho_factors = 10 ** seq(0, -3, length.out = nrho) * rho_max
  alpha_factors = seq(0.1, 0.9, length.out = nalpha)
  alpha_factors = c(0.01, 0.05, alpha_factors, 0.95, 0.99, 1.0)
  nalpha = length(alpha_factors)
  params = vector("list", nlambda * nrho)
  fit = list()
  j   = 1
  for (ialpha in 1:nalpha) {
    alpha = alpha_factors[ialpha]
    for (irho in 1:nrho) {
      rho = rho_factors[irho] / alpha
      for (ilambda in 1:nlambda) {
        if (ilambda == 1) {
          fit[[j]] = SSEMeNet(X = X, Y = Y, B = B, F = F, Sk = Sk, sigma2 = sigma2,
                              lambda = Inf, alpha = alpha, rho = rho, Wb = Wb, Wf = Wf,
                              p = p, maxit = 200, threshold = 1e-4, trans = TRUE, verbose = FALSE)
          lambda_max = initLambdaSSEM2B(X, Y, fit[[j]], Wb)
          lambda_factors = 10 ** seq(0, -3, length.out = nlambda) * lambda_max
          cat(sprintf("SSEM@lambda = %4f, alpha = %4f, rho = %4f\n", lambda_factors[ilambda], alpha, rho))
        } else {
          fit[[j]] = SSEMeNet(X = X, Y = Y, B = B, F = F, Sk = Sk, sigma2 = sigma2,
                              lambda = lambda_factors[ilambda], alpha = alpha, rho = rho, Wb = Wb, Wf = Wf,
                              p = p, maxit = 200, threshold = 1e-4, trans = TRUE, verbose = FALSE)
          cat(sprintf("SSEM@lambda = %4f, alpha = %4f, rho = %4f\n", lambda_factors[ilambda], alpha, rho))
        }
        err = bayesianInfocriterion(X, Y, fit[[j]]$B, fit[[j]]$F, fit[[j]]$mu, fit[[j]]$Det, fit[[j]]$sigma2, p, c)
        params[[j]] = c(lambda_factors[ilambda], alpha, rho, err)
        j = j + 1
      }
    }
  }
  BICmat = data.frame(do.call(rbind, params))
  colnames(BICmat) = c("lambda", "alpha", "rho", "BIC")
  BICmin = which.min(BICmat$BIC)
  list(lambda = BICmat[BICmin, 1], alpha = BICmat[BICmin, 2], rho = BICmat[BICmin, 3],
       fit = fit[[BICmin]], minBIC = min(BICmat$BIC), BICs = BICmat)
}



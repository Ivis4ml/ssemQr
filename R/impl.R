##' @title implSSEM
##' @description implementor function of SSEM solver for eQTL mapping
##' @param data Data archive of experiment measurements, includeing
##' eQTL matrices, Gene expression matrices of different conditions,
##' marker of eQTLs and data generation SEM model
##' @return List of TPR and FDR
##' @export
implSSEM = function(data = NULL) {
  gamma = cv.ridgeRegression(
    data$Data$X,
    data$Data$Y,
    data$Data$Sk,
    ngamma = 50,
    nfold = 5,
    data$Vars$n,
    data$Vars$p,
    data$Vars$q
  )

  fit = ridgeRegression(
    data$Data$X,
    data$Data$Y,
    data$Data$Sk,
    gamma,
    data$Vars$n,
    data$Vars$p,
    data$Vars$q,
    trans = FALSE,
    sparse = FALSE
  )

  optBIC = opt.SSEMiPALM(
    X = data$Data$X,
    Y = data$Data$Y,
    B = fit$B,
    F = fit$F,
    Sk = data$Data$Sk,
    sigma2 = fit$sigma2,
    nlambda = 10,
    nrho = 10,
    p = data$Vars$p,
    wt = TRUE
  )
  fitc = SSEMiPALM(
    X = data$Data$X,
    Y = data$Data$X,
    B = fit$B,
    F = fit$F,
    Sk = data$Data$Sk,
    sigma2 = fit$sigma2,
    lambda = optBIC$lambda,
    rho = optBIC$rho,
    Wb = 1 / abs(fit$B),
    Wf = 1 / abs(fit$F),
    p = Ng,
    maxit = 200,
    trans = TRUE,
    strict = TRUE
  )

  res = data.frame(
    TPRB = TPR(fitc$B, data$Vars$B),
    FDRB = FDR(fitc$B, data$Vars$B),
    TPRF = TPR(fitc$F, data$Vars$F),
    FDRF = FDR(fitc$F, data$Vars$F)
  )
}

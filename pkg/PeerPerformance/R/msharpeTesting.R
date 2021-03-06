####################################################################################
## Set of R functions for the modified Sharpe ratio testing
####################################################################################

####################################################################################
## Test modified Sharpe ratio difference 
####################################################################################

.msharpeTesting = function(x, y, level = 0.90, na.neg = TRUE, control = list()){
  
  x = as.matrix(x)
  y = as.matrix(y)
  
  # process control parameters
  ctr = processControl(control)
  
  # check if enough data are available for testing
  dxy = x - y
  idx = (!is.nan(dxy)&!is.na(dxy))
  rets = cbind(x[idx], y[idx])
  T = sum(idx)
  if (T < ctr$minObs){
    stop("intersection of 'x' and 'y' is shorter than 'minObs'")
  }
  
  # msharpe testing
  if (ctr$type == 1) {
    # ==> asymptotic approach
    tmp = msharpeTestAsymptotic(rets, level, na.neg, ctr$hac, ctr$ttype)
  }
  else {
    # ==> bootstrap approach (iid and circular block bootstrap)
    if (ctr$bBoot == 0){
      ctr$bBoot = msharpeBlockSize(x, y, level, na.neg, ctr)
    }
    bsids = bootIndices(T, ctr$nBoot, ctr$bBoot)
    tmp = msharpeTestBootstrap(rets, level, na.neg, bsids, ctr$bBoot, ctr$ttype, ctr$pBoot)
  }
  
  # info on the funds  
  info = infoFund(rets, level = level, na.neg = na.neg)
  
  ## form output
  out = list(n = T, 
             msharpe  = info$msharpe,
             dmsharpe = as.vector(tmp$dmsharpe), 
             tstat    = as.vector(tmp$tstat), 
             pval     = as.vector(tmp$pval))  
  return(out)
  
}
msharpeTesting = compiler::cmpfun(.msharpeTesting)

####################################################################################
## Test modified Sharpe ratio difference using the asymptotic approach
####################################################################################

# difference of sharpe ratios
.msharpe.ratio.diff = function(X, Y = NULL, level, na.neg, ttype){
  if (is.null(Y)){
    Y = X[,2,drop=FALSE]
    X = X[,1,drop=FALSE]
  }
  m1X   = colMeans(X)
  m1Y   = colMeans(Y)
  X_    = sweep(x = X, MARGIN = 2, STATS = m1X, FUN = "-")
  Y_    = sweep(x = Y, MARGIN = 2, STATS = m1Y, FUN = "-")
  m2X   = colMeans(X_^2)
  m2Y   = colMeans(Y_^2)
  m3X   = colMeans(X_^3)
  m3Y   = colMeans(Y_^3)
  m4X   = colMeans(X_^4)
  m4Y   = colMeans(Y_^4)
  za    = qnorm(1 - level)
  skewX = m3X / m2X^(3/2)
  skewY = m3Y / m2Y^(3/2)
  kurtX = (m4X / m2X^2) - 3
  kurtY = (m4Y / m2Y^2) - 3
  mVaRX = -m1X + sqrt(m2X) * (-za - (1/6) * (za^2 - 1) * skewX - (1/24) * (za^3 - 3*za) * kurtX + (1/36) * (2*za^3 - 5*za) * skewX^2)
  mVaRY = -m1Y + sqrt(m2Y) * (-za - (1/6) * (za^2 - 1) * skewY - (1/24) * (za^3 - 3*za) * kurtY + (1/36) * (2*za^3 - 5*za) * skewY^2)
  if (na.neg){
    mVaRX[mVaRX < 0] = NA
    mVaRY[mVaRY < 0] = NA
  }
  if (ttype == 1){
    # test based on quotient
    mSR1 = m1X / mVaRX
    mSR2 = m1Y / mVaRY
  }
  else{
    # test based on product
    mSR1 = m1X * mVaRY
    mSR2 = m1Y * mVaRX
  }
  diff = mSR1 - mSR2
  return(diff)

}
msharpe.ratio.diff = compiler::cmpfun(.msharpe.ratio.diff)

.msharpeTestAsymptotic = function(rets, level, na.neg, hac, ttype){
  
  dmsharpe = msharpe.ratio.diff(rets, Y = NULL, level, na.neg, ttype)
  if (is.na(dmsharpe)){
    out = list(dmsharpe = NA, tstat = NA, se = NA, pval = NA)
    return(out)
  }
  se    = se.msharpe.asymptotic(rets, level, hac, ttype)
  tstat = dmsharpe / se
  pval  = 2 * pnorm(-abs(tstat)) # asymptotic normal p-value
  out   = list(dmsharpe = dmsharpe, tstat = tstat, se = se, pval = pval)
  return(out)
  
}
msharpeTestAsymptotic = compiler::cmpfun(.msharpeTestAsymptotic)

####################################################################################
# Asymptotic standard error
####################################################################################

.se.msharpe.asymptotic = function(X, level, hac, ttype){
  
  # estimation of (robust) Psi function; see Ledoit Wolf paper
  compute.Psi.hat = function(V.hat, hac){
    if (hac){
      T = length(V.hat[, 1])
      alpha.hat = compute.alpha.hat(V.hat)
      S.star = min(2.6614 * (alpha.hat * T)^0.2, T) # DA fix here to avoid >T values
      Psi.hat = compute.Gamma.hat(V.hat, 0)
      j = 1
      while (j < S.star){
        Gamma.hat = compute.Gamma.hat(V.hat, j)
        Psi.hat = Psi.hat + kernel.Parzen(j/S.star) * (Gamma.hat + t(Gamma.hat))
        j = j + 1
      }
      Psi.hat = (T/(T - 4)) * Psi.hat
    }
    else{
      Psi.hat = cov(V.hat)
    }
    return(Psi.hat)
    
  }
  
  # Parzen kernel
  kernel.Parzen = function(x){
    if (abs(x) <= 0.5) 
      result = 1 - 6 * x^2 + 6 * abs(x)^3
    else if (abs(x) <= 1) 
      result = 2 * (1 - abs(x))^3
    else result = 0
    return(result)
  }
  
  compute.alpha.hat = function(V.hat){
    p = ncol(V.hat)
    num = den = 0
    for (i in 1:p){
      fit = ar(V.hat[, i], 0, 1, method = "ols")
      rho.hat = as.numeric(fit[2])
      sig.hat = sqrt(as.numeric(fit[3]))
      num = num + 4 * rho.hat^2 * sig.hat^4/(1 - rho.hat)^8
      den = den + sig.hat^4/(1 - rho.hat)^4
    }
    return(num/den)
    
  }
  
  compute.Gamma.hat = function (V.hat, j){
    T = nrow(V.hat)
    p = ncol(V.hat)
    Gamma.hat = matrix(0, p, p)
    if (j >= T) 
      stop("j must be smaller than the row dimension!")
    for (i in ((j + 1):T)){
      Gamma.hat = Gamma.hat + tcrossprod(V.hat[i, ], V.hat[i - j, ]) 
    }
    Gamma.hat = Gamma.hat/T
    return(Gamma.hat)
    
  }
  
  T = nrow(X)
  m1   = colMeans(X)
  X_   = sweep(x = X, MARGIN = 2, STATS = m1, FUN = "-")
  m2   = colMeans(X_^2)
  m3   = colMeans(X_^3)
  m4   = colMeans(X_^4)
  g2   = m2 + m1^2
  g3   = m3 + 3 * m1 * g2 - 2 * m1^3
  g4   = m4 + 4 * m1 * g3 - 6 * m1^2 * g2 + 3 * m1^4
  skew = m3 / m2^(3/2)
  kurt = (m4 / m2^2) - 3
  
  # gradient for underlying moments
  dm1i = c(1, 0, 0, 0, 0, 0, 0, 0)
  dm1j = c(0, 0, 0, 0, 1, 0, 0, 0)
  dm2i = c(-2 * m1[1], 1, 0, 0, 0, 0, 0, 0)
  dm2j = c(0, 0, 0, 0, -2 * m1[2], 1, 0, 0)
  dm3i = c(-3 * g2[1] + 6 * m1[1]^2, -3 * m1[1], 1, 0, 0, 0, 0, 0)
  dm3j = c(0, 0, 0, 0, -3 * g2[2] + 6 * m1[2]^2, -3 * m1[2], 1, 0)
  dm4i = c(-4 * g3[1] + 12 * m1[1] * g2[1] - 12 * m1[1]^3, 6 * m1[1]^2, -4 * m1[1], 1, 0, 0, 0, 0)
  dm4j = c(0, 0, 0, 0, -4 * g3[2] + 12 * m1[2] * g2[2] - 12 * m1[2]^3, 6 * m1[2]^2, -4 * m1[2], 1)
  
  dsi = (m2[1]^(3/2) * dm3i - 1.5 * m3[1] * m2[1]^(1/2) * dm2i) / m2[1]^3
  dsj = (m2[2]^(3/2) * dm3j - 1.5 * m3[2] * m2[2]^(1/2) * dm2j) / m2[2]^3
  
  dki = (m2[1]^2 * dm4i - 2 * m4[1] * m2[1] * dm2i) / m2[1]^4
  dkj = (m2[2]^2 * dm4j - 2 * m4[2] * m2[2] * dm2j) / m2[2]^4
  
  za = qnorm(1 - level)
  
  # constants for speedup
  cst1i = 1 / (2 * sqrt(m2[1]))
  cst1j = 1 / (2 * sqrt(m2[2]))
  cst2i = sqrt(m2[1])
  cst2j = sqrt(m2[2])
  
  # dmVaRi
  dmVaRi = -dm1i - cst1i * dm2i * za 
  dmVaRi = dmVaRi - (1/6) * (za^2 - 1) * (cst1i * skew[1] * dm2i + cst2i * dsi)
  dmVaRi = dmVaRi + (1/36) * (2*za^3 - 5*za) * (cst1i * skew[1]^2 * dm2i + 2 * cst2i * skew[1] * dsi)
  dmVaRi = dmVaRi - (1/24) * (za^3 - 3*za) * (cst1i * kurt[1] * dm2i + cst2i * dki)
  
  # dmVaRj
  dmVaRj = -dm1j - cst1j * dm2j * za
  dmVaRj = dmVaRj - (1/6) * (za^2 - 1) * (cst1j * skew[2] * dm2j + cst2j * dsj)
  dmVaRj = dmVaRj + (1/36) * (2*za^3 - 5*za) * (cst1j * skew[2]^2 * dm2j + 2 * cst2j * skew[2] * dsj)
  dmVaRj = dmVaRj - (1/24) * (za^3 - 3*za) * (cst1j * kurt[2] * dm2j + cst2j * dkj)
  
  # mVaR
  mVaR = -m1 + sqrt(m2) * (-za - (1/6) * (za^2 - 1) * skew - (1/24) * (za^3 - 3*za) * kurt + (1/36) * (2*za^3 - 5*za) * skew^2)
  
  # gradient
  if (ttype == 1){
    # test based on quotient
    tmp1 = (mVaR[1] * dm1i - m1[1] * dmVaRi) / mVaR[1]^2
    tmp2 = (mVaR[2] * dm1j - m1[2] * dmVaRj) / mVaR[2]^2
  }
  else{
    # test based on product
    tmp1 = (mVaR[2] * dm1i + m1[1] * dmVaRj)
    tmp2 = (mVaR[1] * dm1j + m1[2] * dmVaRi)
  }
  
  gradient = tmp1 - tmp2
  V.hat = matrix(NA, T, 8)
  V.hat[,c(1,5)] = sweep(x = X, MARGIN = 2, STATS = m1, FUN = "-")
  V.hat[,c(2,6)] = sweep(x = X^2, MARGIN = 2, STATS = g2, FUN = "-")
  V.hat[,c(3,7)] = sweep(x = X^3, MARGIN = 2, STATS = g3, FUN = "-")
  V.hat[,c(4,8)] = sweep(x = X^4, MARGIN = 2, STATS = g4, FUN = "-")
  Psi.hat = compute.Psi.hat(V.hat, hac)
  se = as.numeric( sqrt(crossprod(gradient, Psi.hat %*% gradient) / T) )
  return(se)
  
}
se.msharpe.asymptotic = compiler::cmpfun(.se.msharpe.asymptotic)

####################################################################################
## Test modified Sharpe ratio difference using circular studentized boostrap of Ledoit and Wolf
####################################################################################

.msharpeTestBootstrap = function(rets, level, na.neg, bsids, b, ttype, pBoot, d = 0){
  
  T = nrow(rets)
  x = rets[,1,drop=FALSE]
  y = rets[,2,drop=FALSE]
  
  dmsharpe = as.numeric( msharpe.ratio.diff(x, y, level, na.neg, ttype) - d )
  if (is.na(dmsharpe)){
    out = list(dmsharpe = NA, tstat = NA, se = NA, bststat = NA, pval = NA)
    return(out)
  }
  se = se.msharpe.bootstrap(x, y, level, b, ttype)
  #se = se.msharpe.asymptotic(X = cbind(x, y), level = level, hac = TRUE, ttype = ttype)
  
  # bootstrap indices
  nBoot = ncol(bsids)
  bsidx = 1 + bsids %% T # ensure that the bootstrap indices match the length of the time series
  bsX   = matrix(x[bsidx], T, nBoot)
  bsY   = matrix(y[bsidx], T, nBoot)
  
  bsdmsharpe = msharpe.ratio.diff(bsX, bsY, level, na.neg = FALSE, ttype) # DA consider negative mVaR as well in the bootstrap
  bsse       = se.msharpe.bootstrap(bsX, bsY, level, b, ttype)
  tstat      = dmsharpe / se
  
  if (pBoot == 1){
    # first type p-value calculation
    bststat = abs(bsdmsharpe - dmsharpe) / bsse
    pval    = (sum(bststat >= abs(tstat)) + 1) / (nBoot + 1)
    #pval = sum(bststat >= abs(tstat)) / nBoot
  }
  else{
    # second type p-value calculation (as in Barras)
    bststat = (bsdmsharpe - dmsharpe) / bsse
    pval    = 2 * min(sum(bststat > tstat) + 1, sum(bststat < tstat) + 1) / (nBoot + 1)
    #pval = 2 * min(sum(bststat > tstat), sum(bststat < tstat)) / nBoot
  }
 
  out = list(dmsharpe = dmsharpe, tstat = tstat, se = se, bststat = bststat, pval = pval)
  return(out)
}
msharpeTestBootstrap = compiler::cmpfun(.msharpeTestBootstrap)

####################################################################################
# Bootstrap standard error
####################################################################################

.se.msharpe.bootstrap = function(X, Y, level, b, ttype) {
  
  ## Compute Psi with two approaches: 1) iid bootstrap, 2) circular block bootstrap
  compute.Psi.hat = function(V.hat, b) {
    T = length(V.hat[, 1])
    if (b == 1) {
      # ==> standard estimation
      Psi.hat = cov(V.hat)
    }
    else {
      # ==> block estimation
      l = floor(T/b)
      Psi.hat = matrix(0, 8, 8)
      for (j in (1:l)) {
        zeta = b^0.5 * colMeans(V.hat[((j - 1) * b + 1):(j * b),,drop=FALSE])
        Psi.hat = Psi.hat + tcrossprod(zeta) 
      }
      Psi.hat = Psi.hat / l
    }
    return(Psi.hat)
  }
  
  T = nrow(X)
  N = ncol(Y)
  m1X   = colMeans(X)
  m1Y   = colMeans(Y)
  X_    = sweep(x = X, MARGIN = 2, STATS = m1X, FUN = "-")
  Y_    = sweep(x = Y, MARGIN = 2, STATS = m1Y, FUN = "-")
  m2X   = colMeans(X_^2)
  m2Y   = colMeans(Y_^2)
  m3X   = colMeans(X_^3)
  m3Y   = colMeans(Y_^3)
  m4X   = colMeans(X_^4)
  m4Y   = colMeans(Y_^4)
  g2X   = m2X + m1X^2
  g2Y   = m2Y + m1Y^2
  g3X   = m3X + 3 * m1X * g2X - 2 * m1X^3
  g3Y   = m3Y + 3 * m1Y * g2Y - 2 * m1Y^3
  g4X   = m4X + 4 * m1X * g3X - 6 * m1X^2 * g2X + 3 * m1X^4 
  g4Y   = m4Y + 4 * m1Y * g3Y - 6 * m1Y^2 * g2Y + 3 * m1Y^4
  
  dm1X = matrix(rep(c(1, 0, 0, 0, 0, 0, 0, 0), N), 8, N, FALSE)
  dm1Y = matrix(rep(c(0, 0, 0, 0, 1, 0, 0, 0), N), 8, N, FALSE)
  dm2X = rbind(-2 * m1X, 1, 0, 0, 0, 0, 0, 0)
  dm2Y = rbind(0, 0, 0, 0, -2 * m1Y, 1, 0, 0)
  dm3X = rbind(-3 * g2X + 6 * m1X^2, -3 * m1X, 1, 0, 0, 0, 0 ,0)
  dm3Y = rbind(0, 0, 0, 0, -3 * g2Y + 6 * m1Y^2, -3 * m1Y, 1, 0)
  dm4X = rbind(-4 * g3X + 12 * m1X * g2X - 12 * m1X^3, 6 * m1X^2, -4 * m1X, 1, 0, 0, 0, 0)
  dm4Y = rbind(0, 0, 0, 0, -4 * g3Y + 12 * m1Y * g2Y - 12 * m1Y^3, 6 * m1Y^2, -4 * m1Y, 1)
  
  # matrix form
  m1X_  = matrix(m1X, nrow = 8, ncol = N, byrow = TRUE)
  m1Y_  = matrix(m1Y, nrow = 8, ncol = N, byrow = TRUE)
  m2X_  = matrix(m2X, nrow = 8, ncol = N, byrow = TRUE)
  m2Y_  = matrix(m2Y, nrow = 8, ncol = N, byrow = TRUE)
  m3X_  = matrix(m3X, nrow = 8, ncol = N, byrow = TRUE)
  m3Y_  = matrix(m3Y, nrow = 8, ncol = N, byrow = TRUE)
  m4X_  = matrix(m4X, nrow = 8, ncol = N, byrow = TRUE)
  m4Y_  = matrix(m4Y, nrow = 8, ncol = N, byrow = TRUE)
  
  dm1X_ = matrix(dm1X, nrow = 8, ncol = N, byrow = FALSE)
  dm1Y_ = matrix(dm1Y, nrow = 8, ncol = N, byrow = FALSE)
  dm2X_ = matrix(dm2X, nrow = 8, ncol = N, byrow = FALSE) 
  dm2Y_ = matrix(dm2Y, nrow = 8, ncol = N, byrow = FALSE) 
  dm3X_ = matrix(dm3X, nrow = 8, ncol = N, byrow = FALSE) 
  dm3Y_ = matrix(dm3Y, nrow = 8, ncol = N, byrow = FALSE)
  dm4X_ = matrix(dm4X, nrow = 8, ncol = N, byrow = FALSE) 
  dm4Y_ = matrix(dm4Y, nrow = 8, ncol = N, byrow = FALSE)
  
  dsX_ = (m2X_^(3/2) * dm3X_ - 1.5 * m3X_ * m2X_^(1/2) * dm2X_) / m2X_^3
  dsY_ = (m2Y_^(3/2) * dm3Y_ - 1.5 * m3Y_ * m2Y_^(1/2) * dm2Y_) / m2Y_^3
  dkX_ = (m2X_^2 * dm4X_ - 2 * m4X_ * m2X_ * dm2X_) / m2X_^4
  dkY_ = (m2Y_^2 * dm4Y_ - 2 * m4Y_ * m2Y_ * dm2Y_) / m2Y_^4
  
  skewX_ = m3X_ / m2X_^(3/2)
  skewY_ = m3Y_ / m2Y_^(3/2)
  kurtX_ = (m4X_ / m2X_^2) - 3
  kurtY_ = (m4Y_ / m2Y_^2) - 3
  
  za = qnorm(1 - level)
  
  # constants for speedup
  cst1X_ = 1 / (2 * sqrt(m2X_))
  cst1Y_ = 1 / (2 * sqrt(m2Y_))
  cst2X_ = sqrt(m2X_)
  cst2Y_ = sqrt(m2Y_)
  
  dmVaRX_ = -dm1X_ - za * cst1X_ * dm2X_ 
  dmVaRX_ = dmVaRX_ - (1/6) * (za^2 - 1) * (cst1X_ * skewX_ * dm2X_ + cst2X_ * dsX_)
  dmVaRX_ = dmVaRX_ + (1/36) * (2*za^3 - 5*za) * (cst1X_ * skewX_^2 * dm2X_ + 2 * cst2X_ * skewX_ * dsX_)
  dmVaRX_ = dmVaRX_ - (1/24) * (za^3 - 3*za) * (cst1X_ * kurtX_ * dm2X_ + cst2X_ * dkX_)
  
  # dmVaRY
  dmVaRY_ = -dm1Y_ - za * cst1Y_ * dm2Y_ 
  dmVaRY_ = dmVaRY_ - (1/6) * (za^2 - 1) * (cst1Y_ * skewY_ * dm2Y_ + cst2Y_ * dsY_)
  dmVaRY_ = dmVaRY_ + (1/36) * (2*za^3 - 5*za) * (cst1Y_ * skewY_^2 * dm2Y_ + 2 * cst2Y_ * skewY_ * dsY_)
  dmVaRY_ = dmVaRY_ - (1/24) * (za^3 - 3*za) * (cst1Y_ * kurtY_ * dm2Y_ + cst2Y_ * dkY_)
  
  # mVaR
  mVaRX_ = -m1X_ + sqrt(m2X_) * (-za - (1/6) * (za^2 - 1) * skewX_ - (1/24) * (za^3 - 3*za) * kurtX_ + (1/36) * (2*za^3 - 5*za) * skewX_^2)
  mVaRY_ = -m1Y_ + sqrt(m2Y_) * (-za - (1/6) * (za^2 - 1) * skewY_ - (1/24) * (za^3 - 3*za) * kurtY_ + (1/36) * (2*za^3 - 5*za) * skewY_^2)
  
  # gradient
  if (ttype == 1){
    # test based on quotient
    tmp1 = (mVaRX_ * dm1X_ - m1X_ * dmVaRX_) / mVaRX_^2
    tmp2 = (mVaRY_ * dm1Y_ - m1Y_ * dmVaRY_) / mVaRY_^2
  }
  else{
    # test based on product
    tmp1 = (mVaRY_ * dm1X_ + m1X_ * dmVaRY_)
    tmp2 = (mVaRX_ * dm1Y_ + m1Y_ * dmVaRX_)
  }
  
  # =======
  gradient = array(NA, c(8, 1, N))
  gradient[1:8,1,] = tmp1 - tmp2
  V.hat = array(NA, c(T, 8, N))
  V.hat[,1,] = sweep(x = X, MARGIN = 2, STATS = m1X, FUN = "-")
  V.hat[,5,] = sweep(x = Y, MARGIN = 2, STATS = m1Y, FUN = "-")
  V.hat[,2,] = sweep(x = X^2, MARGIN = 2, STATS = g2X, FUN = "-")
  V.hat[,6,] = sweep(x = Y^2, MARGIN = 2, STATS = g2Y, FUN = "-")
  V.hat[,3,] = sweep(x = X^3, MARGIN = 2, STATS = g3X, FUN = "-")
  V.hat[,7,] = sweep(x = Y^3, MARGIN = 2, STATS = g3Y, FUN = "-")
  V.hat[,4,] = sweep(x = X^4, MARGIN = 2, STATS = g4X, FUN = "-")
  V.hat[,8,] = sweep(x = Y^4, MARGIN = 2, STATS = g4Y, FUN = "-")
  Psi.hat = array(apply(X = V.hat, MARGIN = 3, FUN = compute.Psi.hat, b = b), c(8, 8, N))
  se = vector('double', N)
  for (i in 1 : N) {
    se[i] = sqrt(crossprod(gradient[,,i], Psi.hat[,,i] %*% gradient[,,i]) / T)
  }
  return(se)
}
se.msharpe.bootstrap = compiler::cmpfun(.se.msharpe.bootstrap)

####################################################################################
## Set of R functions used to compute pi and lambda
####################################################################################

####################################################################################
# Given a value of lambda, this function returns the pi0, pi+ and pi- of the funds
####################################################################################

# pval is N x N, where N is the number of funds
# if we want to estimate pi for a single fund, we need to provide pval as a 1 x N matrix

computePi = function (pval, dalpha, lambda, nBoot){
  
  # remove full NA
  ntot = nrow(pval)
  pos.na = which(apply(is.na(pval), 1, all))
  pos.ok = (1:ntot)
  npos.na = length(pos.na)
  if (npos.na > 0) {
    pos.ok = pos.ok[-pos.na]
    pval = pval[-pos.na, -pos.na]
    dalpha = dalpha[-pos.na, -pos.na]
    if (!is.null(lambda)) {
      lambda = lambda[-pos.na]
    }
  }
  if (is.null(lambda)) {
    lambda = PeerPerformance:::computeOptLambda(pval = pval, nBoot = nBoot)
  }
  pizero = PeerPerformance:::computePizero(pval, lambda)
  n1  = ncol(pval)
  n   = max(n1 - 1, 1)
  ni0 = pizero * n
  pipos = pineg = rep(0, nrow(pval))
  idx = (0 <= pizero) & (pizero < 1)
  pipos[idx] = ((1/n) * (rowSums(dalpha >= 0, na.rm = TRUE) - 0.5 * ni0))[idx]
  pineg[idx] = ((1/n) * (rowSums(dalpha < 0, na.rm = TRUE) - 0.5 * ni0))[idx]
  idx = pizero == 0
  pipos[idx & apply(dalpha >= 0, 1, all, na.rm = TRUE)] = 1
  pineg[idx & apply(dalpha < 0, 1, all, na.rm = TRUE)]  = 1
  idx = pipos < 0
  pipos[idx] = 0
  pineg[idx] = 1 - pizero[idx]
  idx = pineg < 0
  pineg[idx] = 0
  pipos[idx] = 1 - pizero[idx]
  idx = pipos > 1
  pipos[idx] = 1 - pizero[idx]
  pineg[idx] = 0
  idx = pineg > 1
  pineg[idx] = 1 - pizero[idx]
  pipos[idx] = 0
  # check if a delta is remaining
  delta = pizero + pipos + pineg - 1
  pipos = pipos - (delta/2)
  pineg = pineg - (delta/2)
  # internal assignation
  pizero_ = pipos_ = pineg_ = lambda_ = rep(NA, ntot)
  pizero_[pos.ok] = pizero
  pipos_[pos.ok] = pipos
  pineg_[pos.ok] = pineg
  lambda_[pos.ok] = lambda
  out = list(pizero = pizero_, pipos = pipos_, pineg = pineg_, 
             lambda = lambda_)
  return(out)
  
}

####################################################################################
# Given a value of lambda, this function returns the pizero of the funds
####################################################################################

computePizero = function(pval, lambda){
  
  n1 = nrow(pval)
  n2 = ncol(pval)
  mlambda = matrix(lambda, n1, n2, byrow = FALSE)
  
  pizero = rowMeans(pval >= mlambda, na.rm = TRUE)
  pizero = pizero * (1 / (1 - lambda))
  pizero[pizero > 1] = 1
  
  return(pizero)
}

####################################################################################
# Compute optimal lamba values 
####################################################################################

.computeOptLambda = function(pval, nBoot){
  
  n = nrow(pval) # number of funds
  
  vlambda = seq(0.3, 0.7, 0.02)
  nvlambda = length(vlambda)
  
  mpizero = matrix(data = NA, nrow = n, ncol = nvlambda)
  for (i in 1 : nvlambda){
    mpizero[,i] = computePizero(pval, vlambda[i])
  }
  # pi0hat
  vminpizero = apply(mpizero, 1, 'min')
  
  idx  = !is.nan(pval) & !is.na(pval)
  nObs = rowSums(idx)
  
  bsunif = matrix(runif(max(nObs) * nBoot), max(nObs), nBoot)    
  optlambda = rep(0.5, n)
  
  for (i in 1 : n) {
    nObsi = nObs[i]
    if (nObsi > 0) {
      bsidx = ceiling(bsunif[1:nObsi, 1:nBoot] * nObsi)
      pvalb = na.omit(pval[i,])
      pvalb = matrix(pvalb[bsidx], nrow = nBoot)
      
      mpizerob = matrix(data = NA, nrow = nBoot, ncol = nvlambda)
      for (j in 1 : nvlambda){
        mpizerob[,j] = computePizero(pvalb, vlambda[j])
      }
      
      vMSE = colSums((mpizerob - vminpizero[i])^2)
      
      optlambda[i] = vlambda[which.min(vMSE)]
    } 
  }
  
  return(optlambda)
}
computeOptLambda = cmpfun(.computeOptLambda)
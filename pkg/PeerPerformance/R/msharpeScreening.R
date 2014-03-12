####################################################################################
## Set of R functions for modified Sharpe screening
####################################################################################

msharpeScreening = function(X, rf = 0, level = 0.90, na.neg = TRUE, control = list()) {
  
  # process control
  ctr = processControl(control)
  
  # size of inputs and outputs
  T = nrow(X) 
  N = ncol(X)
  pval = dmsharpe = matrix(data = NA, N, N)
  if (length(rf) == 1){
    rf = rep(rf, N)
  }
  
  # determine which pairs can be compared (in a matrix way)
  Y = 1 * (!is.nan(X) & !is.na(X))
  YY = crossprod(Y) #YY = t(Y) %*% Y # row i indicates how many observations in common with column k
  YY[YY < ctr$minObs] = 0
  YY[YY > 0] = 1 
  liststocks = c(1:nrow(YY))[rowSums(YY) > ctr$minObsPi]  
  
  # determine bootstrap indices (do it before to speed up computations)
  bsids = bootIndices(T, ctr$nBoot, ctr$bBoot)
  
  if (length(liststocks) > 1){
    cl = makeCluster(c(rep("localhost", ctr$nCore)), type = "SOCK")
    
    liststocks = liststocks[1:(length(liststocks)-1)]
    
    z <- clusterApply(cl     = cl, 
                      x      = as.list(liststocks), 
                      fun    = msharpeScreeningi, 
                      rdata  = X,
                      rf     = rf,
                      level  = level,
                      T      = T, 
                      N      = N, 
                      na.neg = na.neg,
                      nBoot  = ctr$nBoot, 
                      bsids  = bsids, 
                      minObs = ctr$minObs,
                      type   = ctr$type,
                      hac    = ctr$hac,
                      b      = ctr$bBoot,
                      ttype  = ctr$ttype,
                      pBoot  = ctr$pBoot)
    stopCluster(cl)
    
    for (i in 1:length(liststocks)) {
      out = z[[i]]  
      id  = liststocks[i]
      pval[id,id:N]     = pval[id:N,id] = out[[2]][id:N]
      dmsharpe[id,id:N] = out[[1]][id:N]
      dmsharpe[id:N,id] = -out[[1]][id:N]
    }
  }
  
  # pi
  pi = computePi(pval = pval, dalpha = dmsharpe, lambda = ctr$lambda, nBoot = ctr$nBoot)
  
  # info on the funds  
  info = infoFund(X, rf = rf, level = level, na.neg = na.neg)
  
  # form output
  out = list(n        = info$nObs, 
             npeer    = colSums(!is.na(pval)),
             msharpe  = info$msharpe, 
             dmsharpe = dmsharpe, 
             pval     = pval,
             lambda   = pi$lambda,
             pizero   = pi$pizero,
             pipos    = pi$pipos,
             pineg    = pi$pineg)
  
  return(out)
}

## Sharpe ratio screening for fund i again its peers
.msharpeScreeningi = function(i, rdata, rf, level, T, N, nBoot, bsids, minObs, na.neg, type, hac, b, ttype, pBoot) {
  
  nPeer = N - i
  X = matrix(rdata[,i], nrow = T, ncol = nPeer)
  Y = matrix(rdata[, (i+1):N], nrow = T, ncol = nPeer)
  #rf.x = rf[i]
  #rf.y = rf[(i+1):N]
  
  dXY = X - Y
  idx = (!is.nan(dXY)&!is.na(dXY))
  X[!idx] = NA
  Y[!idx] = NA    
  nObs = colSums(idx)
  
  pvali = dmsharpei = rep(NA, N)
  
  k = 0
  for (j in (i + 1) : N) {
    k = k + 1
    if (nObs[k] < minObs) {
      next
    }
    rets = cbind(X[idx[,k],1], Y[idx[,k],k])
    
    if (type == 1) {
      #tmp = msharpeTestAsymptotic(rets, rf.x, rf.y[k], level, na.neg, hac, ttype)
      tmp = msharpeTestAsymptotic(rets, rf.x = 0, rf.y = 0, level, na.neg, hac, ttype)
    }
    else {
      #tmp = msharpeTestBootstrap(rets, rf.x, rf.y[k], level, na.neg, bsids, b, ttype, pBoot)
      tmp = msharpeTestBootstrap(rets, rf.x = 0, rf.y = 0, level, na.neg, bsids, b, ttype, pBoot)
    }
    
    dmsharpei[j] = tmp$dmsharpe
    pvali[j] = tmp$pval
  }
  
  out = list(dmsharpei = dmsharpei, pvali = pvali)
  return(out)   
}
msharpeScreeningi = cmpfun(.msharpeScreeningi)
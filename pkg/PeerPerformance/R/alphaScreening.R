####################################################################################
## Set of R function for alpha screening
####################################################################################

alphaScreening = function(X, factors = NULL, control = list()) {  
  
  # process control
  ctr = processControl(control)
  
  T = nrow(X) 
  N = ncol(X)
  pval = dalpha = matrix(data = NA, N, N)
  
  Y = 1 * (!is.nan(X) & !is.na(X))
  #YY = t(Y) %*% Y # row i indicates how many observations in common with column k
  YY = crossprod(Y)
  YY[YY < ctr$minObs] = 0
  YY[YY > 0] = 1
  
  liststocks = c(1:nrow(YY))[rowSums(YY) > ctr$minObsPi]  
  
  if (length(liststocks) > 1){
    cl = makeCluster(c(rep("localhost", ctr$nCore)), type = "SOCK")
    
    liststocks = liststocks[1:(length(liststocks)-1)]
    
    z <- clusterApply(cl      = cl, 
                      x       = as.list(liststocks), 
                      fun     = alphaScreeningi, 
                      rdata   = X, 
                      factors = factors, 
                      T       = T, 
                      N       = N)
    
    stopCluster(cl)
    
    for (i in 1:length(liststocks)){
      out = z[[i]]  
      id  = liststocks[i]
      pval[id,id:N]   = pval[id:N,id] = out[[2]][id:N]
      dalpha[id,id:N] = out[[1]][id:N]
      dalpha[id:N,id] = -out[[1]][id:N]
    }
  }
  
  # pi
  pi = computePi(pval = pval, dalpha = dalpha, lambda = ctr$lambda, nBoot = ctr$nBoot)
  
  # info on the funds  
  info = infoFund(X)
  
  # form output
  out = list(n       = info$nObs, 
             npeer   = colSums(!is.na(pval)),
             alpha   = info$alpha, 
             dalpha  = dalpha, 
             pval    = pval,
             lambda  = pi$lambda,
             pizero  = pi$pizero,
             pipos   = pi$pipos,
             pineg   = pi$pineg)
  
  return(out)
}

.alphaScreeningi = function(i, rdata, factors, T, N) {
  pvali = dalphai = rep(NA, N)
  
  nPeer = N - i
  X = matrix(rdata[,i], nrow = T, ncol = nPeer)
  Y = matrix(rdata[,(i+1):N], nrow = T, ncol = nPeer)
  
  dXY = X - Y
  idx = (!is.nan(dXY)&!is.na(dXY))
  X[!idx] = NA
  Y[!idx] = NA    
  #nObs = colSums(idx)
  
  if (is.null(factors)){
    fit = lm(dXY ~ 1)
  }
  else {
    fit = lm(dXY ~ 1 + factors)
  }
  sumfit = summary(fit)
  if (nPeer == 1) {
    pvali[N]   = sumfit$coef[1,4]
    dalphai[N] = sumfit$coef[1,1]
  } 
  else {
    k = 1
    for (j in (i + 1) : N) {
      pvali[j]   = sumfit[[k]]$coef[1,4]
      dalphai[j] = sumfit[[k]]$coef[1,1]
      k = k + 1
    }
  }
  
  out = list(dalphai = dalphai, pvali = pvali)
  return(out)
} 
alphaScreeningi = cmpfun(.alphaScreeningi)
library('HUM')
library('VUROCS')
library('mcca')


##################
### Function to calculate the asymptotic variance of HUM estimator
## Input:
# s_dat: a list, which contains the vectors of biomarker values for individual class labels (the order of classes in this list should be matched with 'cat_ord' and 'cat_freq')
# cat_freq: a vector of the category sample size
# humhat: the obtained HUM estimator value
## Output:
# s.square: the asymptotic variance of HUM estimator
asymvar.f <- function(s_dat, cat_freq, humhat) {
  
  n <- sum(cat_freq)
  M <- length(s_dat)
  
  w_dat <- NULL
  layer_in <- s_dat[[1]]
  layer_in <- layer_in[order(layer_in)]
  attri_in <- rep(1, cat_freq[1])
  w_dat[[1]] <- attri_in
  
  for (k in 2:M) {
    
    layer_out <- s_dat[[k]]
    layer_out <- layer_out[order(layer_out)]
    attri_out <- rep(0, cat_freq[k])
    
    j <- 1
    while (j <= cat_freq[k-1] && layer_in[j] < layer_out[1]) {
      attri_out[1] <- attri_out[1] + attri_in[j]
      j <- j + 1
    }
    
    for (i in 2:cat_freq[k]) {
      attri_out[i] <- attri_out[i-1]
      while (j <= cat_freq[k-1] && layer_in[j] < layer_out[i]) {
        attri_out[i] <- attri_out[i] + attri_in[j]
        j <- j + 1
      }
    }
    
    layer_in <- layer_out
    attri_in <- attri_out
    
    w_dat[[k]] <- attri_out
    
  }
  
  wt_dat <- NULL
  layer_in <- s_dat[[M]]
  layer_in <- layer_in[order(layer_in)]
  attri_in <- rep(1, cat_freq[M])
  wt_dat[[M]] <- attri_in
  
  for (k in (M-1):1) {
    
    layer_out <- s_dat[[k]]
    layer_out <- layer_out[order(layer_out)]
    attri_out <- rep(0, cat_freq[k])
    
    j <- cat_freq[k+1]
    while (j >= 1 && layer_in[j] > layer_out[cat_freq[k]]) {
      attri_out[cat_freq[k]] <- attri_out[cat_freq[k]] + attri_in[j]
      j <- j - 1
    }
    
    for (i in (cat_freq[k]-1):1) {
      attri_out[i] <- attri_out[i+1]
      while (j >= 1 && layer_in[j] > layer_out[i]) {
        attri_out[i] <- attri_out[i] + attri_in[j]
        j <- j - 1
      }
    }
    
    layer_in <- layer_out
    attri_in <- attri_out
    
    wt_dat[[k]] <- attri_out
    
  }
  
  L <- Reduce(c, w_dat) * Reduce(c, wt_dat)  # L: the number of M-path that passes through the v-th subject/vertex 
  
  s.square <- n * sum(L^2) / (prod(cat_freq))^2 - n * humhat^2 * sum(1/cat_freq)
  
  return(s.square)
  
}




##################
### Function to calculate the HUM estimate
## Input:
# dat: a dataframe of the dataset with two variables (class: the category index, X: the value of biomarker)
# cat_ord: a vector of the order of categories
# cat_freq: a vector of the category sample size
## Output:
# humC0: the HUM estimate
hum.f <- function(cat_ord, dat, cat_freq) {
  
  n <- sum(cat_freq)
  M <- length(cat_ord)
  
  layer_in <- dat$X[dat$class == cat_ord[1]]
  layer_in <- layer_in[order(layer_in)]
  attri_in <- rep(1, cat_freq[1])
  
  for (k in 2:M) {
    
    layer_out <- dat$X[dat$class == cat_ord[k]]
    layer_out <- layer_out[order(layer_out)]
    attri_out <- rep(0, cat_freq[k])
    
    j <- 1
    while (j <= cat_freq[k-1] && layer_in[j] < layer_out[1]) {
      attri_out[1] <- attri_out[1] + attri_in[j]
      j <- j + 1
    }
    
    for (i in 2:cat_freq[k]) {
      attri_out[i] <- attri_out[i-1]
      while (j <= cat_freq[k-1] && layer_in[j] < layer_out[i]) {
        attri_out[i] <- attri_out[i] + attri_in[j]
        j <- j + 1
      }
    }
    
    layer_in <- layer_out
    attri_in <- attri_out
    
  }
  
  humC0 <- sum(attri_out) / prod(cat_freq)
  
  return(humC0)
}




############################ generate data

###Example 1
dat.gen_ex1 <- function(M=4, nn) {
  
  xmean <- c(1, 2, 3, 5)
  xvar <- c(0.5, 1, 2, 1)
  xx <- lapply(1:M, function(i) rnorm(nn[i], xmean[i], sqrt(xvar[i])))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 2
dat.gen_ex2 <- function(M=4, nn) {
  
  xmean <- seq(1, M, 1)
  xx <- lapply(1:M, function(i) rnorm(nn[i], xmean[i], sqrt(0.5)))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 3
dat.gen_ex3 <- function(M=4, nn) {
  
  lambda <- c(8, 4, 2, 1/2, 1/4, 1/8)
  xx <- lapply(1:M, function(i) rexp(nn[i], lambda[i]))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 4
dat.gen_ex4 <- function(M=4, nn) {
  
  gshape <- c(1, 1, 2, 2, 3)
  gscale <- c(0.5, 1, 1, 2, 2)
  xx <- lapply(1:M, function(i) rgamma(n = nn[i], shape = gshape[i], scale = gscale[i]))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 5
dat.gen_ex5 <- function(M=4, nn) {
  
  xmean <- c(1, 2, 3, 5)
  xvar <- c(0.5, 1, 2, 1)
  xx <- lapply(1:M, function(i) rlnorm(nn[i], xmean[i], sqrt(xvar[i])))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 6
dat.gen_ex6 <- function(M=4, nn) {
  
  xshape <- c(0.5, 1, 2, 5)
  xscale <- c(1, 2, 3, 5)
  xx <- lapply(1:M, function(i) rweibull(nn[i], xshape[i], xscale[i]))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 7
dat.gen_ex7 <- function(M=4, nn) {
  
  xdf <- c(1, 2, 4, 8)
  xx <- lapply(1:M, function(i) rchisq(nn[i], xdf[i]))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}



sim.fun <- function(M = 4, csize = 20, nn=c(), dat.gen) {
  
  if (length(nn)==0) nn <- sapply(1:M, function(x) csize + rbinom(1, floor(csize/2), 0.5)) 
  
  ZZ <- dat.gen(M=M, nn=nn)
  
  dat <- ZZ$dat
  cat_ord <- ZZ$cat_ord
  
  t <- as.data.frame(table(dat$class))
  loc <- match(t[[1]], cat_ord)
  cat_freq <- t$Freq[loc]
  
  n <- sum(cat_freq)
  
  
  ## compute HUM
  
  # our method
  {
    ptm <- proc.time()
    humC0 <- hum.f(cat_ord = cat_ord, dat = dat, cat_freq = cat_freq)
    comptime0 <- unname(proc.time() - ptm)[1:3]
  }
  
  # HUM package 
  {  
    s_data <- NULL
    prodValue <- 1
    for(i in 1:length(cat_ord)){
      ind <- which(dat$class == cat_ord[i])
      vrem <- sort(dat$X[ind])
      s_data <- c(s_data, list(vrem))
      prodValue=prodValue*length(ind)
    }
    seqAll <- matrix(1:M, 1, M)
    ptm <- proc.time()
    humC_HUM <- CalcGene(s_data, seqAll, prodValue)$HUM
    comptime_HUM <- unname(proc.time() - ptm)[1:3]
  }
  
  # VUROCS package
  {
    ptm <- proc.time()
    humC_VUROCS <- VUS(y = dat$class, fx = dat$X)$val
    comptime_VUROCS <- unname(proc.time() - ptm)[1:3]
  }  
  
  # mcca package
  {
    ptm <- proc.time()
    humC_mcca <- hum(dat$class, dat$X, method = 'multinom')$measure
    comptime_mcca <- unname(proc.time() - ptm)[1:3]
  }
  
  
  ## compute asymptotic variance
  
  # Our method
  {
    s_dat <- NULL
    for(i in 1:length(cat_ord)){
      ind <- which(dat$class == cat_ord[i])
      vrem <- dat$X[ind]
      s_dat <- c(s_dat, list(vrem))
    }
    ptm <- proc.time()
    s.square <- asymvar.f(s_dat, cat_freq, humC0) # the asymptotic variance
    asymp.sd <- sqrt(s.square / sum(cat_freq) ) # the asymptotic sd
    comptimeSE <- unname(proc.time() - ptm)[1:3]
  }
  
  # VUROCS package
  {
    ptm <- proc.time()
    se_VUROCS <- sqrt(VUSvar(y = dat$class, fx = dat$X)$var)
    comptimeSE_VUROCS <- unname(proc.time() - ptm)[1:3]
  }
  
  
  
  return(list(humC0=humC0, comptime0=comptime0, 
              humC_HUM=humC_HUM, comptime_HUM=comptime_HUM,
              humC_VUROCS=humC_VUROCS, comptime_VUROCS=comptime_VUROCS, 
              humC_mcca=humC_mcca, comptime_mcca=comptime_mcca,
              asymp.sd = asymp.sd, comptimeSE = comptimeSE, 
              se_VUROCS=se_VUROCS, comptimeSE_VUROCS=comptimeSE_VUROCS))
  
}





res=NULL
for (i in 1:500) {
  sim.res <- sim.fun(M = 4, csize = 20, dat.gen = dat.gen_ex1)  # to simulate Example 1
#  sim.res <- sim.fun(M = 4, nn = c(100, 70, 30, 10), dat.gen = dat.gen_ex1)  # to simulate Example 7
  res[[i]]=sim.res
  print(i)
}


Comptime0 <- apply(t(sapply(res, function(x) {x$comptime0})), 2, sum)
Comptime_HUM <- apply(t(sapply(res, function(x) {x$comptime_HUM})), 2, sum)
Comptime_VUROCS <- apply(t(sapply(res, function(x) {x$comptime_VUROCS})), 2, sum)
Comptime_mcca <- apply(t(sapply(res, function(x) {x$comptime_mcca})), 2, sum)

ComptimeSE <- apply(t(sapply(res, function(x) {x$comptimeSE})), 2, sum)
ComptimeSE_VUROCS <- apply(t(sapply(res, function(x) {x$comptimeSE_VUROCS})), 2, sum)

HumC0 <- sapply(res, function(x) {x$humC0})
HumC_HUM <- sapply(res, function(x) {x$humC_HUM})
HumC_VUROCS <- sapply(res, function(x) {x$humC_VUROCS})
HumC_mcca <- sapply(res, function(x) {x$humC_mcca})

SE <- sapply(res, function(x) {x$asymp.sd})
SE_VUROCS <- sapply(res, function(x) {x$se_VUROCS})




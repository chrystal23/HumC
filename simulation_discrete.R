library('HUM')
library('VUROCS')
library('mcca')


##################
### Function to calculate the asymptotic variance of HUM estimator considering tied events
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
  attri_in <- matrix(1, nrow = cat_freq[1], ncol = 2)
  w_dat[[1]] <- attri_in
  
  for (k in 2:M) {
    
    layer_out <- s_dat[[k]]
    layer_out <- layer_out[order(layer_out)]
    attri_out <- matrix(0, nrow = cat_freq[k], ncol = k+1)
    
    j <- 1
    while (j <= cat_freq[k-1] && layer_in[j] < layer_out[1]) {
      attri_out[1,2] <- attri_out[1,2] + attri_in[j,1]
      j <- j + 1
    }
    while (j <= cat_freq[k-1] && layer_in[j] == layer_out[1]) {
      attri_out[1,-(1:2)] <- attri_out[1,-(1:2)] + attri_in[j,-1] / (2:k)
      j <- j + 1
    }
    
    for (i in 2:cat_freq[k]) {
      if(layer_out[i]==layer_out[i-1]) {
        attri_out[i,] <- attri_out[i-1,]
        next
      } 
      attri_out[i,2] <- attri_out[i-1,2] + sum(attri_out[i-1,-(1:2)] * (2:k))
      while (j <= cat_freq[k-1] && (layer_in[j] < layer_out[i])) {
        attri_out[i,2] <- attri_out[i,2] + attri_in[j,1]
        j <- j + 1
      }
      while (j <= cat_freq[k-1] && (layer_in[j] == layer_out[i])) {
        attri_out[i,-(1:2)] <- attri_out[i,-(1:2)] + attri_in[j,-1] / (2:k)
        j <- j + 1
      }
    }
    
    attri_out[,1] <- rowSums(attri_out[,-1])
    layer_in <- layer_out
    attri_in <- attri_out
    
    w_dat[[k]] <- attri_out
    
  }
  
  wt_dat <- NULL
  layer_in <- s_dat[[M]]
  layer_in <- layer_in[order(layer_in)]
  attri_in <- matrix(1, nrow = cat_freq[M], ncol = 2)
  wt_dat[[M]] <- attri_in
  
  for (k in (M-1):1) {
    
    layer_out <- s_dat[[k]]
    layer_out <- layer_out[order(layer_out)]
    attri_out <- matrix(0, nrow = cat_freq[k], ncol = M-k+2)
    
    j <- cat_freq[k+1]
    while (j >= 1 && layer_in[j] > layer_out[cat_freq[k]]) {
      attri_out[cat_freq[k],2] <- attri_out[cat_freq[k],2] + attri_in[j,1]
      j <- j - 1
    }
    while (j >= 1 && layer_in[j] == layer_out[cat_freq[k]]) {
      attri_out[cat_freq[k],-(1:2)] <- attri_out[cat_freq[k],-(1:2)] + attri_in[j,-1] / (2:(M-k+1))
      j <- j - 1
    }
    
    for (i in (cat_freq[k]-1):1) {
      if(layer_out[i]==layer_out[i+1]) {
        attri_out[i,] <- attri_out[i+1,]
        next
      } 
      attri_out[i,2] <- attri_out[i+1,2] + sum(attri_out[i+1,-(1:2)] * (2:(M-k+1)))
      while (j >= 1 && layer_in[j] > layer_out[i]) {
        attri_out[i,2] <- attri_out[i,2] + attri_in[j,1]
        j <- j - 1
      }
      while (j >= 1 && layer_in[j] == layer_out[i]) {
        attri_out[i,-(1:2)] <- attri_out[i,-(1:2)] + attri_in[j,-1] / (2:(M-k+1))
        j <- j - 1
      }
    }
    
    attri_out[,1] <- rowSums(attri_out[,-1])
    layer_in <- layer_out
    attri_in <- attri_out
    
    wt_dat[[k]] <- attri_out
    
  }
  
  L <- sapply(1:M, function(d) {
    if (d==M) A <- matrix(1, nrow = M, ncol = 1)
    else A <- t(sapply(1:d, function(a) sapply(1:(M-d+1), function(b) factorial(a)*factorial(b)/factorial(a+b-1))))
    diag( w_dat[[d]][,-1,drop=F] %*% A %*% t(wt_dat[[d]][,-1,drop=F]) )
  })
  L <- Reduce(c, L) # L: the weighted number of M-path that passes through the v-th subject/vertex 
  
  s.square <- n * sum(L^2) / (prod(cat_freq))^2 - n * humhat^2 * sum(1/cat_freq)
  
  return(s.square)
  
}




##################
### Function to calculate the HUM estimate considering tied events
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
  attri_in <- matrix(1, nrow = cat_freq[1], ncol = 2)
  
  for (k in 2:M) {
    
    layer_out <- dat$X[dat$class == cat_ord[k]]
    layer_out <- layer_out[order(layer_out)]
    attri_out <- matrix(0, nrow = cat_freq[k], ncol = k+1)
    
    j <- 1
    while (j <= cat_freq[k-1] && layer_in[j] < layer_out[1]) {
      attri_out[1,2] <- attri_out[1,2] + attri_in[j,1]
      j <- j + 1
    }
    while (j <= cat_freq[k-1] && layer_in[j] == layer_out[1]) {
      attri_out[1,-(1:2)] <- attri_out[1,-(1:2)] + attri_in[j,-1] / (2:k)
      j <- j + 1
    }
    
    for (i in 2:cat_freq[k]) {
      if(layer_out[i]==layer_out[i-1]) {
        attri_out[i,] <- attri_out[i-1,]
        next
      } 
      attri_out[i,2] <- attri_out[i-1,2] + sum(attri_out[i-1,-(1:2)] * (2:k))
      while (j <= cat_freq[k-1] && (layer_in[j] < layer_out[i])) {
        attri_out[i,2] <- attri_out[i,2] + attri_in[j,1]
        j <- j + 1
      }
      while (j <= cat_freq[k-1] && (layer_in[j] == layer_out[i])) {
        attri_out[i,-(1:2)] <- attri_out[i,-(1:2)] + attri_in[j,-1] / (2:k)
        j <- j + 1
      }
    }
    
    attri_out[,1] <- rowSums(attri_out[,-1])
    layer_in <- layer_out
    attri_in <- attri_out
    
  }
  
  humC0 <- sum(attri_out[,1]) / prod(cat_freq)
  
  return(humC0)
}




############################ generate data

###Example 8
dat.gen_ex8 <- function(M=4, nn) {
  
  p <- c(1/16, 1/8, 1/4, 1/2)
  xx <- lapply(1:M, function(i) rbinom(nn[i], 40, p[i]))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 9
dat.gen_ex9 <- function(M=4, nn) {
  
  p <- 1:8/16
  xx <- lapply(1:M, function(i) rbinom(nn[i], 80, p[i]))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 10
dat.gen_ex10 <- function(M=4, nn) {
  
  p <- seq(0.1,0.9,0.2)
  xx <- lapply(1:M, function(i) rbinom(nn[i], 1, p[i]))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 11
dat.gen_ex11 <- function(M=4, nn) {
  
  lambda <- c(1, 2, 4, 8)
  xx <- lapply(1:M, function(i) rpois(nn[i], lambda[i]))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 12
dat.gen_ex12 <- function(M=4, nn) {
  
  p <- c(1/2, 1/4, 1/8, 1/16, 1/32)
  xx <- lapply(1:M, function(i) rgeom(nn[i], p[i]))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 13
dat.gen_ex13 <- function(M=4, nn) {
  
  m <- c(5, 10, 15, 20)
  xx <- lapply(1:M, function(i) rhyper(nn[i], m[i], 40-m[i], 10))
  
  dat <- c()
  for (i in 1:M) dat <- rbind(dat, cbind(i, xx[[i]]))
  dat <- as.data.frame(dat)
  colnames(dat) <- c('class', 'X')
  
  cat_ord <- 1:M
  
  return(list(dat=dat, cat_ord=cat_ord))
}

###Example 14
dat.gen_ex14 <- function(M=4, nn) {
  
  p <- c(1/2, 1/4, 1/8, 1/16)
  xx <- lapply(1:M, function(i) rnbinom(nn[i], 5, p[i]))
  
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
  sim.res <- sim.fun(M = 4, csize = 20, dat.gen = dat.gen_ex8)  # to simulate Example 8
  #  sim.res <- sim.fun(M = 4, nn = c(100, 70, 30, 10), dat.gen = dat.gen_ex14)  # to simulate Example 14
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




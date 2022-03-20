

library('gtools')

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







######################## BMI outcome

load('exposome.RData')

fac_id <- rep(0, dim(exposome)[2])
for (i in 1:dim(exposome)[2]) { if (is.factor(exposome[[i]])) fac_id[i]=1}
fac_id <- which(fac_id == 1) # id's of all the factor variables

## determine the category order and calculate its corresponding HUM for each exposure variable 
# (results are recorded in the list humC)

dat <- data.frame(class = phenotype$hs_bmi_c_cat)
M <- 4

cat_ord <- 1:M
t <- as.data.frame(table(dat$class))
loc <- match(cat_ord, t[[1]])
cat_freq <- t$Freq[loc]

humC <- NULL

for (i in 2:dim(exposome)[2]) {
  
  if (i %in% fac_id) next
  
  dat$X <- exposome[[i]]
  
  perm <- permutations(M, M, 1:M)
  humperm <- apply(perm, 1, function(x) hum.f(x, dat, cat_freq[x]))
  
  id <- which.max(humperm)
  humC[[i]] <- list(hum=humperm[id], cat_ord=perm[id,])
  
  print(i)
  
}


## pick out the exposures with high HUM (>0.1) (results recorded in the list humC_picked)

humC0 <- sapply(humC[-1], function(x) if (is.null(x$hum)) 0 else x$hum )
hid <- which(humC0>0.1)

humC_picked <- humC[hid+1]
expo_picked <- colnames(exposome)[hid+1]
for (i in 1:length(humC_picked)) {
  humC_picked[[i]]$exposure <- expo_picked[i]
}


#save(humC, humC_picked, file = 'humC_bmi.RData')




## record the computation time

library('HUM')
library('VUROCS')
library('mcca')

datt <- data.frame(class = phenotype$hs_bmi_c_cat)
M <- 4

cat_ord <- 1:M
t <- as.data.frame(table(datt$class))
loc <- match(cat_ord, t[[1]])
cat_freq <- t$Freq[loc]

comptime0 <- rep(0, 3)
comptime_HUM <- rep(0, 3)
comptime_VUROCS <- rep(0, 3)
comptimeSE <- rep(0, 3)
comptimeSE_VUROCS <- rep(0, 3)


humC0 <- c()
humC_HUM <- c()
humC_VUROCS <- c()
asymp.sd <- c()
se_VUROCS <- c()


for (j in 2:dim(exposome)[2]) {
  
  if (j %in% fac_id) {
    humC0 <- c(humC0, NA)
    humC_HUM <- c(humC_HUM, NA)
    humC_VUROCS <- c(humC_VUROCS, NA)
    asymp.sd <- c(asymp.sd, NA)
    se_VUROCS <- c(se_VUROCS, NA)
    next
  }
  
  datt$X <- exposome[[j]]
  cat_ord <- humC[[j]]$cat_ord
  t <- as.data.frame(table(datt$class))
  loc <- match(cat_ord, t[[1]])
  cat_freq <- t$Freq[loc]
  
  ### by HUM package
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
    humC_HUM <- c(humC_HUM, CalcGene(s_data, seqAll, prodValue)$HUM)  # the HUM estimate
    comptime_HUM <- comptime_HUM + unname(proc.time() - ptm)[1:3]  # the computation time
  }
  
  ### by VUROCS package
  {
    datvus <- datt
    vusid <- NULL
    for (m in 1:M) {
      vusid[[m]] <- which(datvus$class==cat_ord[m])
    }
    for (m in 1:M) {
      datvus$class[vusid[[m]]] <- m
    }
    ptm <- proc.time()
    humC_VUROCS <- c(humC_VUROCS, VUS(y = datvus$class, fx = datvus$X)$val)
    comptime_VUROCS <- comptime_VUROCS + unname(proc.time() - ptm)[1:3]
  }
  {
    ptm <- proc.time()
    se_VUROCS <- c(se_VUROCS, sqrt(VUSvar(y = as.numeric(datvus$class), fx = datvus$X)$var))
    comptimeSE_VUROCS <- comptimeSE_VUROCS + unname(proc.time() - ptm)[1:3]
  }
  
  dat <- datvus
  cat_ord <- 1:M
  
  ### By our method
  {
    ptm <- proc.time()
    humC0 <- c(humC0, hum.f(cat_ord = cat_ord, dat = dat, cat_freq = cat_freq))
    comptime0 <- comptime0 + unname(proc.time() - ptm)[1:3]
    
  }
  
  ## compute asymptotic variance
  {
    s_dat <- NULL
    for(i in 1:length(cat_ord)){
      ind <- which(dat$class == cat_ord[i])
      vrem <- dat$X[ind]
      s_dat <- c(s_dat, list(vrem))
    }
    ptm <- proc.time()
    s.square <- asymvar.f(s_dat, cat_freq, humC0[length(humC0)]) # the asymptotic variance
    asymp.sd <- c(asymp.sd, sqrt(s.square / sum(cat_freq) ) ) # the asymptotic sd
    comptimeSE <- comptimeSE + unname(proc.time() - ptm)[1:3]
  }
  
  
  print(j)
  
}


Comptime_bmi <- cbind(comptime0, comptime_HUM, comptime_VUROCS, comptimeSE, comptimeSE_VUROCS)
HumC_bmi <- cbind(humC0, humC_HUM, humC_VUROCS, asymp.sd, se_VUROCS)

#save(Comptime_bmi, HumC_bmi, file = 'time_bmi.RData')



## draw boxplots for the significant exposure biomarkers


library('ggplot2')

boxdata <- cbind(exposome['hs_hcb_cadj_Log2'], phenotype['hs_bmi_c_cat'])
boxdata <- boxdata[boxdata[[1]]>0,]

p <- ggplot(data=boxdata, aes(x=hs_bmi_c_cat, y=hs_hcb_cadj_Log2, fill=hs_bmi_c_cat))
p1 <- p + geom_boxplot() +theme(legend.position = 'none') + ylab('hcb_cadj') + xlab('bmi_cat') #+ coord_cartesian(ylim = c(0, 8))
p1



legd <- cowplot::get_legend(p + geom_boxplot() + ylab('hcb_cadj') + xlab('bmi_cat') + 
                              theme(legend.position = 'top') +
                              scale_fill_discrete(name = 'bmi_cat', labels = c('1 (13)', '2 (904)', '3 (253)', '4 (131)')))



boxdata <- cbind(exposome['hs_sumPCBs5_cadj_Log2'], phenotype['hs_bmi_c_cat'])

p <- ggplot(data=boxdata, aes(x=hs_bmi_c_cat, y=hs_sumPCBs5_cadj_Log2, fill=hs_bmi_c_cat))
p2 <- p + geom_boxplot() +theme(legend.position = 'none') + ylab('sumPCBs5_cadj') + xlab('bmi_cat') #+ coord_cartesian(ylim = c(0, 8))
p2


boxdata <- cbind(exposome['hs_pcb180_cadj_Log2'], phenotype['hs_bmi_c_cat'])

p <- ggplot(data=boxdata, aes(x=hs_bmi_c_cat, y=hs_pcb180_cadj_Log2, fill=hs_bmi_c_cat))
p3 <- p + geom_boxplot() +theme(legend.position = 'none') + ylab('pcb180_cadj') + xlab('bmi_cat')  #+ coord_cartesian(ylim = c(0, 8))
p3

cowplot::plot_grid(legd, cowplot::plot_grid(p1, p2, p3, ncol=3, align = 'h'), nrow = 2, rel_heights = c(1,10) )
ggsave('bmi.eps', width = 18, height = 7, units = 'cm', dpi=800, device = cairo_ps)








rm(list = setdiff(ls(), c('hum.f', 'asymvar.f')) )
gc()

######################## IQ outcome

load('exposome.RData')

fac_id <- rep(0, dim(exposome)[2])
for (i in 1:dim(exposome)[2]) { if (is.factor(exposome[[i]])) fac_id[i]=1}
fac_id <- which(fac_id == 1) # id's of all the factor variables

phenotype$IQ <- phenotype$hs_correct_raven 

for (i in 1:dim(phenotype)[1]) {
  if (phenotype$hs_correct_raven[i] <= 15) phenotype$IQ[i] <- 1
  else if (phenotype$hs_correct_raven[i] <= 20) phenotype$IQ[i] <- 2
  else if (phenotype$hs_correct_raven[i] <= 25) phenotype$IQ[i] <- 3
  else if (phenotype$hs_correct_raven[i] <= 30) phenotype$IQ[i] <- 4
  else phenotype$IQ[i] <- 5
}

phenotype$IQ <- as.factor(phenotype$IQ)


## determine the category order and calculate its corresponding HUM for each exposure variable 
# (results are recorded in the list humC_IQ)

dat_IQ <- data.frame(class = phenotype$IQ)
M <- 5

cat_ord <- 1:M
t <- as.data.frame(table(dat_IQ$class))
loc <- match(cat_ord, t[[1]])
cat_freq_IQ <- t$Freq[loc]

humC_IQ <- NULL
varid <- c()

for (i in 2:dim(exposome)[2]) {
  
  if (i %in% fac_id) next
  
  dat_IQ$X <- exposome[[i]]
  
  perm <- permutations(M, M, 1:M)
  
  humperm <- apply(perm, 1, function(x) hum.f(x, dat_IQ, cat_freq_IQ[x])) 
    
  id <- which.max(humperm)
  humC_IQ[[i]] <- list(hum=humperm[id], cat_ord=perm[id,])
  
  print(i)

}


## pick out the exposures with high HUM (>1/24) (results recorded in the list humC_picked_IQ)

humC0_IQ <- sapply(humC_IQ[-1], function(x) if (is.null(x$hum)) 0 else x$hum )
hid_IQ <- which(humC0_IQ>1/24)

humC_picked_IQ <- humC_IQ[hid_IQ+1]
expo_picked_IQ <- colnames(exposome)[hid_IQ+1]
for (i in 1:length(humC_picked_IQ)) {
  humC_picked_IQ[[i]]$exposure <- expo_picked_IQ[i]
}


#save(humC_IQ, humC_picked_IQ, file = 'humC_IQ.RData')




## record the computation time

library('HUM')
library('VUROCS')
library('mcca')

datt <- data.frame(class = phenotype$IQ)
M <- 5

cat_ord <- 1:M
t <- as.data.frame(table(datt$class))
loc <- match(cat_ord, t[[1]])
cat_freq <- t$Freq[loc]

comptime0 <- rep(0, 3)
comptime_VUROCS <- rep(0, 3)
comptimeSE <- rep(0, 3)
comptimeSE_VUROCS <- rep(0, 3)


humC0 <- c()
humC_VUROCS <- c()
asymp.sd <- c()
se_VUROCS <- c()


for (j in 2:dim(exposome)[2]) {
  
  if (j %in% fac_id) {
    humC0 <- c(humC0, NA)
    humC_VUROCS <- c(humC_VUROCS, NA)
    asymp.sd <- c(asymp.sd, NA)
    se_VUROCS <- c(se_VUROCS, NA)
    next
  }
  
  datt$X <- exposome[[j]]
  cat_ord <- humC_IQ[[j]]$cat_ord
  t <- as.data.frame(table(datt$class))
  loc <- match(cat_ord, t[[1]])
  cat_freq <- t$Freq[loc]
  
  
  ### by VUROCS package
  {
    datvus <- datt
    vusid <- NULL
    for (m in 1:M) {
      vusid[[m]] <- which(datvus$class==cat_ord[m])
    }
    for (m in 1:M) {
      datvus$class[vusid[[m]]] <- m
    }
    ptm <- proc.time()
    humC_VUROCS <- c(humC_VUROCS, VUS(y = datvus$class, fx = datvus$X)$val)
    comptime_VUROCS <- comptime_VUROCS + unname(proc.time() - ptm)[1:3]
  }
  {
    ptm <- proc.time()
    se_VUROCS <- c(se_VUROCS, sqrt(VUSvar(y = as.numeric(datvus$class), fx = datvus$X)$var))
    comptimeSE_VUROCS <- comptimeSE_VUROCS + unname(proc.time() - ptm)[1:3]
  }
  
  dat <- datvus
  cat_ord <- 1:M
  
  ### By our method
  {
    ptm <- proc.time()
    humC0 <- c(humC0, hum.f(cat_ord = cat_ord, dat = dat, cat_freq = cat_freq))
    comptime0 <- comptime0 + unname(proc.time() - ptm)[1:3]
    
  }
  
  ## compute asymptotic variance
  {
    s_dat <- NULL
    for(i in 1:length(cat_ord)){
      ind <- which(dat$class == cat_ord[i])
      vrem <- dat$X[ind]
      s_dat <- c(s_dat, list(vrem))
    }
    ptm <- proc.time()
    s.square <- asymvar.f(s_dat, cat_freq, humC0[length(humC0)]) # the asymptotic variance
    asymp.sd <- c(asymp.sd, sqrt(s.square / sum(cat_freq) ) ) # the asymptotic sd
    comptimeSE <- comptimeSE + unname(proc.time() - ptm)[1:3]
  }
  
  
  print(j)
  
}


Comptime_IQ <- cbind(comptime0, comptime_VUROCS, comptimeSE, comptimeSE_VUROCS)
HumC_IQ <- cbind(humC0, humC_VUROCS, asymp.sd, se_VUROCS)

#save(Comptime_IQ, HumC_IQ, file = 'time_IQ.RData')



#################### draw


library('ggplot2')

boxdata <- cbind(exposome['hs_pcb138_madj_Log2'], phenotype['IQ'])
boxdata <- boxdata[boxdata[[1]]>0,]

p <- ggplot(data=boxdata, aes(x=IQ, y=hs_pcb138_madj_Log2, fill=IQ))
p1 <- p + geom_boxplot() +theme(legend.position = 'none') + ylab('pcb138_madj') + xlab('IQ_cat') #+ coord_cartesian(ylim = c(0, 8))
p1


legd <- cowplot::get_legend(p + geom_boxplot() + ylab('pcb138_madj') + xlab('IQ_cat') + 
                              theme(legend.position = 'top') +
                              scale_fill_discrete(name = 'IQ_cat', labels = c('1 (99)', '2 (174)', '3 (269)', '4 (310)', '5 (449)')))


boxdata <- cbind(exposome['hs_pfhxs_c_Log2'], phenotype['IQ'])

p <- ggplot(data=boxdata, aes(x=IQ, y=hs_pfhxs_c_Log2, fill=IQ))
p2 <- p + geom_boxplot() +theme(legend.position = 'none') + ylab('pfhxs_c') + xlab('IQ_cat') #+ coord_cartesian(ylim = c(0, 8))
p2


boxdata <- cbind(exposome['hs_pm25_yr_hs_h_None'], phenotype['IQ'])

p <- ggplot(data=boxdata, aes(x=IQ, y=hs_pm25_yr_hs_h_None, fill=IQ))
p3 <- p + geom_boxplot() +theme(legend.position = 'none') + ylab('pm25_yr_hs') + xlab('IQ_cat')  #+ coord_cartesian(ylim = c(0, 8))
p3

cowplot::plot_grid(legd, cowplot::plot_grid(p1, p2, p3, ncol=3, align = 'h'), nrow = 2, rel_heights = c(1,10) ) 
ggsave('iq.eps', width = 18, height = 7, units = 'cm', dpi=800, device = cairo_ps)











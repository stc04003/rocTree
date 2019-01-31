library(rocTree)
library(survival)

## check hazard estimation
truehaz1.1 <- function(dat) {
    haz <- with(dat, 2 * Y * exp(0.5 * (z2 + z4 + z6 + z8 + z10) - 0.5 * (z1 + z3 + z5 + z7 + z9)))
    approxfun(dat$Y, haz, "constant", yleft = 0)
    ## stepfun(dat$Y[ord], c(0, haz[ord]))
}

truehaz1.2 <- function(dat) {
    haz <- with(dat, 2 * Y * exp(2 * (z1 + z2)))
    approxfun(dat$Y, haz, "constant", yleft = 0)
    ## stepfun(dat$Y[ord], c(0, haz[ord]))
}

truehaz2.1 <- function(dat) {
    dat <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    haz <- with(dat, 2 * Y  * exp(2 * z1 + 2 * z2))
    approxfun(dat$Y, haz, "constant", yleft = 0)
    ## stepfun(dat$Y[ord], c(0, haz[ord]))
}

truehaz2.2 <- function(dat) {
    dat <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    haz <- with(dat, exp(2 * z1 + 2 * z2))
    approxfun(dat$Y, haz, "constant", yleft = 0)
    ## stepfun(dat$Y[ord], c(0, haz[ord]))
}

do.haz <- function(n, cen, sce = 1.1) {
    dat <- simu(n, cen, sce)
    dat.test <- do.call(rbind, lapply(1:50, function(e) simuTest(dat)))
    dat.test$id <- rep(1:50, each = n)
    ctrl <- list(CV = TRUE, tau = quantile(unique(dat$Y), .95), minsp = 20, minsp2 = 5,
                 parallel = TRUE)
    if (sce %in% c(1.1, 1.2))
        fm <- Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10
    else fm <- Surv(Y, death) ~ z1 + z2
    fit <- rocForest(fm, data = dat, id = id, control = ctrl)
    pred <- predict(fit, dat.test, type = "hazard")    
    pred0 <- predict(fit, dat.test, type = "hazard0")    
    list(dat = dat, fit = fit, pred = pred, pred0 = pred0)
}

set.seed(1)
foo <- do.haz(100, 0, 1.2)


##################################################################################
## Benchmark

do.tree <- function(n, cen, sce = 1.1) {
    dat <- simu(n, cen, sce)
    dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    ctrl <- list(CV = TRUE, tau = quantile(dat0$Y, .95), minsp = 20, minsp2 = 5)
    tt <- seq(0, ctrl$tau, length = 100)
    ## Fitting
    if (sce %in% c(1.1, 1.2))
        fm <- Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10
    else fm <- Surv(Y, death) ~ z1 + z2
    fit <- rocTree(fm, data = dat, id = id, control = ctrl)
    fit.dcon <- rocTree(fm, data = dat, id = id, splitBy = "dCON", control = ctrl)
    print(fit)
    print(fit.dcon)
}

debug(rocTree)
traceback()

set.seed(1)
system.time(do.tree(100, 0, 1.1))

## Root                      
##  ¦--2) z1 <= 0.36000      
##  ¦   ¦--4) z10 <= 0.30000*
##  ¦   °--5) z10 > 0.30000* 
##  °--3) z1 > 0.36000       
##      ¦--6) z2 <= 0.59000* 
##      °--7) z2 > 0.59000*  
## Root                     
##  ¦--2) z1 <= 0.36000*    
##  °--3) z1 > 0.36000      
##      ¦--6) z2 <= 0.59000*
##      °--7) z2 > 0.59000* 
##  user  system elapsed 
## 7.884   0.160   8.081

set.seed(1)
system.time(do.tree(100, 0.25, 1.3))

## Root                     
##  ¦--2) z1 <= 0.33000     
##  ¦   ¦--4) z1 <= 0.10000*
##  ¦   °--5) z1 > 0.10000* 
##  °--3) z1 > 0.33000      
##      ¦--6) z2 <= 0.89000*
##      °--7) z2 > 0.89000* 
## Root                     
##  ¦--2) z1 <= 0.33000     
##  ¦   ¦--4) z1 <= 0.10000*
##  ¦   °--5) z1 > 0.10000* 
##  °--3) z1 > 0.33000*     
##    user  system elapsed 
##   1.395   0.102   1.509     

set.seed(1)
system.time(do.tree(100, 0.25, 2.3))

## Root                     
##  ¦--2) z2 <= 0.77000     
##  ¦   ¦--4) z1 <= 0.51000*
##  ¦   °--5) z1 > 0.51000* 
##  °--3) z2 > 0.77000*     
## Root                 
##  ¦--2) z2 <= 0.77000*
##  °--3) z2 > 0.77000* 
##    user  system elapsed 
##   1.636   0.112   1.765 

set.seed(1)
system.time(do.tree(100, 0.50, 2.5))

## Root                     
##  ¦--2) z1 <= 0.51000*    
##  °--3) z1 > 0.51000      
##      ¦--6) z1 <= 0.72000*
##      °--7) z1 > 0.72000* 
## Root                 
##  ¦--2) z1 <= 0.51000*
##  °--3) z1 > 0.51000* 
##    user  system elapsed 
##   1.325   0.110   1.455 

set.seed(1)
system.time(do.tree(100, 0.50, 1.5))

## Root                 
##  ¦--2) z1 <= 0.80000*
##  °--3) z1 > 0.80000* 
## Root                 
##  ¦--2) z1 <= 0.80000*
##  °--3) z1 > 0.80000* 
##    user  system elapsed 
##   1.151   0.098   1.264 

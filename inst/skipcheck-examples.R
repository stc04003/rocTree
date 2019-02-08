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
## Benchmark -- tree
##################################################################################

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

##################################################################################
## Benchmark -- forest
##################################################################################
 

do.Forest <- function(n, cen, sce = 1.1) {
    dat <- simu(n, cen, sce)
    dat0 <- dat[cumsum(with(dat, unlist(lapply(split(id, id), length), use.names = FALSE))),]
    ctrl <- list(CV = TRUE, tau = quantile(dat0$Y, .95), minsp = 3, minsp2 = 1, B = 100)
    if (sce %in% c(1.1, 1.2))
        fm <- Surv(Y, death) ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10
    else fm <- Surv(Y, death) ~ z1 + z2
    fit <- rocForest(fm, data = dat, id = id, control = ctrl)
    fit.dcon <- rocForest(fm, data = dat, id = id, splitBy = "dCON", control = ctrl)
    print(fit, 1)
    print(fit.dcon, 1)    
}

set.seed(1)
system.time(do.Forest(100, 0, 1.1))
  
## Root                                       
##  ¦--2) z10 <= 0.29897                      
##  ¦   ¦--4) z8 <= 0.22680*                  
##  ¦   °--5) z8 > 0.22680                    
##  ¦       ¦--10) z1 <= 0.10309              
##  ¦       ¦   ¦--20) z7 <= 0.40206*         
##  ¦       ¦   °--21) z7 > 0.40206*          
##  ¦       °--11) z1 > 0.10309*              
##  °--3) z10 > 0.29897                       
##      ¦--6) z4 <= 0.61856                   
##      ¦   ¦--12) z10 <= 0.76289             
##      ¦   ¦   ¦--24) z4 <= 0.60825          
##      ¦   ¦   ¦   ¦--48) z4 <= 0.17526*     
##      ¦   ¦   ¦   °--49) z4 > 0.17526*      
##      ¦   ¦   °--25) z4 > 0.60825*          
##      ¦   °--13) z10 > 0.76289              
##      ¦       ¦--26) z7 <= 0.51546          
##      ¦       ¦   ¦--52) z8 <= 0.72165*     
##      ¦       ¦   °--53) z8 > 0.72165*      
##      ¦       °--27) z7 > 0.51546*          
##      °--7) z4 > 0.61856                    
##          ¦--14) z3 <= 0.40206              
##          °--15) z3 > 0.40206               
##              ¦--30) z3 <= 0.88660          
##              ¦   ¦--60) z4 <= 0.83505      
##              ¦   ¦   ¦--120) z4 <= 0.82474*
##              ¦   ¦   °--121) z4 > 0.82474* 
##              ¦   °--61) z4 > 0.83505*      
##              °--31) z3 > 0.88660*          
  
## Root                                  
##  ¦--2) z1 <= 0.54000                  
##  ¦   ¦--4) z4 <= 0.68000              
##  ¦   ¦   ¦--8) z1 <= 0.27000          
##  ¦   ¦   ¦   ¦--16) z2 <= 0.30000*    
##  ¦   ¦   ¦   °--17) z2 > 0.30000*     
##  ¦   ¦   °--9) z1 > 0.27000           
##  ¦   ¦       ¦--18) z2 <= 0.47000*    
##  ¦   ¦       °--19) z2 > 0.47000*     
##  ¦   °--5) z4 > 0.68000               
##  ¦       ¦--10) z7 <= 0.65000*        
##  ¦       °--11) z7 > 0.65000          
##  ¦           ¦--22) z3 <= 0.65000*    
##  ¦           °--23) z3 > 0.65000*     
##  °--3) z1 > 0.54000                   
##      ¦--6) z10 <= 0.68000             
##      ¦   ¦--12) z10 <= 0.26000        
##      ¦   ¦   ¦--24) z2 <= 0.67000*    
##      ¦   ¦   °--25) z2 > 0.67000*     
##      ¦   °--13) z10 > 0.26000*        
##      °--7) z10 > 0.68000              
##          ¦--14) z7 <= 0.69000         
##          ¦   ¦--28) z5 <= 0.64000     
##          ¦   ¦   ¦--56) z3 <= 0.11000*
##          ¦   ¦   °--57) z3 > 0.11000* 
##          ¦   °--29) z5 > 0.64000*     
##          °--15) z7 > 0.69000          
##    user  system elapsed 
##  12.820   0.080  12.906 

set.seed(1)
system.time(do.Forest(100, 0, 1.3))

## Root                                   
##  ¦--2) z1 <= 0.27473                   
##  ¦   ¦--4) z1 <= 0.14286               
##  ¦   ¦   ¦--8) z2 <= 0.51648           
##  ¦   ¦   °--9) z2 > 0.51648            
##  ¦   °--5) z1 > 0.14286                
##  ¦       ¦--10) z1 <= 0.17582          
##  ¦       °--11) z1 > 0.17582           
##  ¦           ¦--22) z1 <= 0.24176*     
##  ¦           °--23) z1 > 0.24176*      
##  °--3) z1 > 0.27473                    
##      ¦--6) z1 <= 0.82418               
##      ¦   ¦--12) z2 <= 0.64835          
##      ¦   ¦   ¦--24) z1 <= 0.56044      
##      ¦   ¦   ¦   ¦--48) z2 <= 0.054945*
##      ¦   ¦   ¦   °--49) z2 > 0.054945  
##      ¦   ¦   °--25) z1 > 0.56044*      
##      ¦   °--13) z2 > 0.64835           
##      ¦       ¦--26) z2 <= 0.83516*     
##      ¦       °--27) z2 > 0.83516*      
##      °--7) z1 > 0.82418                
##          ¦--14) z2 <= 0.31868*         
##          °--15) z2 > 0.31868           
##              ¦--30) z2 <= 0.92308*     
##              °--31) z2 > 0.92308*        
## Root                                   
##  ¦--2) z1 <= 0.42000                   
##  ¦   ¦--4) z2 <= 0.60000               
##  ¦   ¦   ¦--8) z2 <= 0.38000           
##  ¦   ¦   ¦   ¦--16) z2 <= 0.27000*     
##  ¦   ¦   ¦   °--17) z2 > 0.27000*      
##  ¦   ¦   °--9) z2 > 0.38000            
##  ¦   ¦       ¦--18) z2 <= 0.50000*     
##  ¦   ¦       °--19) z2 > 0.50000       
##  ¦   ¦           ¦--38) z1 <= 0.040000*
##  ¦   ¦           °--39) z1 > 0.040000* 
##  ¦   °--5) z2 > 0.60000                
##  ¦       ¦--10) z2 <= 0.81000          
##  ¦       ¦   ¦--20) z1 <= 0.18000      
##  ¦       ¦   °--21) z1 > 0.18000*      
##  ¦       °--11) z2 > 0.81000           
##  ¦           ¦--22) z2 <= 0.85000*     
##  ¦           °--23) z2 > 0.85000*      
##  °--3) z1 > 0.42000                    
##      ¦--6) z2 <= 0.66000               
##      ¦   ¦--12) z2 <= 0.25000*         
##      ¦   °--13) z2 > 0.25000           
##      °--7) z2 > 0.66000                
##          ¦--14) z2 <= 0.72000*         
##          °--15) z2 > 0.72000           

##    user  system elapsed 
##   6.244   0.072   6.327 

set.seed(1)
system.time(do.Forest(100, 0.25, 2.3))

## Root                              
##  ¦--2) z2 <= 0.40206              
##  ¦   ¦--4) z1 <= 0.43299          
##  ¦   °--5) z1 > 0.43299           
##  ¦       ¦--10) z2 <= 0.24742*    
##  ¦       °--11) z2 > 0.24742*     
##  °--3) z2 > 0.40206               
##      ¦--6) z2 <= 0.80412          
##      ¦   ¦--12) z2 <= 0.59794*    
##      ¦   °--13) z2 > 0.59794*     
##      °--7) z2 > 0.80412           
##          ¦--14) z1 <= 0.17526*    
##          °--15) z1 > 0.17526      
##              ¦--30) z2 <= 0.81443*
##              °--31) z2 > 0.81443  
  
## Root                                  
##  ¦--2) z2 <= 0.40816                  
##  ¦   ¦--4) z1 <= 0.27551              
##  ¦   °--5) z1 > 0.27551*              
##  °--3) z2 > 0.40816                   
##      ¦--6) z2 <= 0.88776              
##      ¦   ¦--12) z2 <= 0.54082         
##      ¦   ¦   ¦--24) z1 <= 0.19388*    
##      ¦   ¦   °--25) z1 > 0.19388*     
##      ¦   °--13) z2 > 0.54082          
##      ¦       ¦--26) z2 <= 0.80612     
##      ¦       ¦   ¦--52) z2 <= 0.73469*
##      ¦       ¦   °--53) z2 > 0.73469* 
##      ¦       °--27) z2 > 0.80612*     
##      °--7) z2 > 0.88776               
##          ¦--14) z2 <= 0.93878         
##          °--15) z2 > 0.93878*         

##    user  system elapsed 
##   4.820   0.096   4.924 

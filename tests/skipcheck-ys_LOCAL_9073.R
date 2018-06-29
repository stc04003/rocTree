#######################################################################
## Load package
#######################################################################

library(rocTree)
library(survival)

#######################################################################
## Codes to generate data; scenario III on page 18
#######################################################################

datGen <- function(n) {
    ## 3 covariates, X1 and x2 are used in the paper, X3 is a noise
    lambda <- .1
    P <- 2
    U <- runif(n)
    X <- matrix(runif(n * P, 0, 1), ncol = P, nrow = n)
    a <- X[, 1]
    b <- X[, 2]
    betat <- 2
    beta <- 1
    ks <- runif(n, 1, 2)
    it <- runif(n, 1, 2)
    T <- log(1 + (betat * ks * (-log(U))) / lambda / exp(beta * a + betat * it)) / betat / ks
    C <- runif(n, 0, .8)
    Y <- pmin(T, C)
    E <- T <= C ## censoring indicator
    dat <- NULL ## ID, Time, Status, X1, X2, X3
    for (i in 1:n) {
        dat <- rbind(dat, cbind(rep(i, rank(Y)[i]), Y[which(Y <= Y[i])], 0, a[i],
                                Y[which(Y <= Y[i])] * ks[i] + it[i], runif(1)))
    }
    dat <- data.frame(dat)
    colnames(dat) <- c("ID", "Time", "Status", "X1", "X2", "X3")
    dat <- dat[order(dat$ID, dat$Time),]
    rownames(dat) <- NULL
    dat$Status[cumsum(aggregate(Time ~ ID, dat, length)[,2])] <- E
    dat
}

n <- 100
set.seed(1234)
dat <- datGen(n)

system.time(foo <- rocTree(Surv(Time, Status) ~ X1 + X2, id = ID, data = dat))
##   user  system elapsed 
##  0.293   0.032   0.345 

foo
##  ROC-guided survival tree
##  node), split
##    * denotes terminal node
## Root                          
##  ¦--2) X2 <= 0.70000          
##  ¦   ¦--4) X1 <= 0.50000*     
##  ¦   °--5) X1 > 0.50000*      
##  °--3) X2 > 0.70000           
##      ¦--6) X1 <= 0.38000*     
##      °--7) X1 > 0.38000       
##          ¦--14) X1 <= 0.71000*
##          °--15) X1 > 0.71000* 

system.time(foo1 <- rocTree(Surv(Time, Status) ~ X1 + X2, id = ID, data = dat, control = list(split.method = "dCON")))
##   user  system elapsed 
##  0.273   0.036   0.314

foo1
##  ROC-guided survival tree
##  node), split
##    * denotes terminal node
## Root                          
##  ¦--2) X2 <= 0.70000          
##  ¦   ¦--4) X1 <= 0.50000*     
##  ¦   °--5) X1 > 0.50000*      
##  °--3) X2 > 0.70000           
##      ¦--6) X1 <= 0.71000      
##      ¦   ¦--12) X2 <= 0.91000*
##      ¦   °--13) X2 > 0.91000* 
##      °--7) X1 > 0.71000*

## ------------------------------------------------------------------------------------------------
## Timing comparison
## ------------------------------------------------------------------------------------------------
do <- function(n, smethod) {
    dat <- datGen(n)
    time1 <- system.time(rocTree(Surv(Time, Status) ~ X1 + X2, id = ID, data = dat))[3]
    time2 <- system.time(rocTree(Surv(Time, Status) ~ X1 + X2, id = ID, data = dat, control = list(split.method = "dCON")))[3]
    c(time1, time2)
}

set.seed(1234)
summary(t(replicate(100, do(100))))

## Not much gain in timing
## elapsed          elapsed      
## Min.   :0.2590   Min.   :0.2640  
## 1st Qu.:0.3272   1st Qu.:0.3147  
## Median :0.3510   Median :0.3440  
## Mean   :0.3561   Mean   :0.3542  
## 3rd Qu.:0.3862   3rd Qu.:0.3820  
## Max.   :0.5610   Max.   :0.5720

set.seed(1234)
summary(t(replicate(100, do(200))))
## elapsed         elapsed     
## Min.   :1.958   Min.   :2.159  
## 1st Qu.:2.955   1st Qu.:2.835  
## Median :3.307   Median :3.129  
## Mean   :3.297   Mean   :3.171  
## 3rd Qu.:3.618   3rd Qu.:3.461  
## Max.   :4.722   Max.   :4.840  

## ------------------------------------------------------------------------------------------------
## With CV
## ------------------------------------------------------------------------------------------------

set.seed(1)
system.time(foo <- rocTree(Surv(Time, Status) ~ X1 + X2, id = ID, data = dat, control = list(CV = TRUE, nflds = 10)))

## Parallel computing is only faster for large n or large number of nflds
set.seed(1)
system.time(foo2 <- rocTree(Surv(Time, Status) ~ X1 + X2, id = ID, data = dat, control = list(CV = TRUE, nflds = 10, parallel = TRUE)))

identical(foo, foo2) ## True

## Print tree: 3 ways to print
foo2 ## default
print(foo2, classic = TRUE)
print(foo2, classic = TRUE, top2bottom = TRUE)

## Plot tree: few different options
plot(foo2)
plot(foo2, control = list(rankdir = "LR"))
plot(foo2, control = list(nodeOnly = TRUE))
plot(foo2, control = list(nodeOnly = TRUE, rankdir = "LR"))
plot(foo2, control = list(nodeOnly = TRUE, shape = "box"))
plot(foo2, control = list(nodeOnly = TRUE, shape = "box", rankdir = "LR"))

plot(foo2, output = "visNetwork")
plot(foo2, output = "visNetwork", control = list(nodeOnly = TRUE, shape = "box"))

## Plot hazard functions
plotTreeHaz(foo2)



###############################################################################################
## fitting rpart and party; using covariates at baseline.
###############################################################################################

## -------------------------------------------------------------------------------------------
## Prepare data for rpart::rtree and partykit::ctree
## -------------------------------------------------------------------------------------------

n <- 100
n3 <- 1000
tt <- seq(0, .5, length = 100)
set.seed(1234)
dat <- datGen(n)[,-6]
dat0 <- dat[cumsum(with(dat, unlist(lapply(split(ID, ID), length), use.names = FALSE))), ]
set.seed(1234)
dat.test <- datGen(n3)[,-6]
dat.test0 <- dat.test[cumsum(with(dat.test,unlist(lapply(split(ID, ID), length), use.names = FALSE))),]
dat.test0[,2:3] <- NA
rownames(dat0) <- rownames(dat.test) <- rownames(dat.test0) <- NULL

## -------------------------------------------------------------------------------------------
## fit rocTree
## -------------------------------------------------------------------------------------------
system.time(foo <- rocTree(Surv(Time, Status) ~ X1 + X2, id = ID, data = dat))

foo.pred <- predict(foo, dat)
ggplot(foo.pred, aes(x = Time, y = Surv)) + geom_smooth()

e

## -------------------------------------------------------------------------------------------
## Fitting rpart::rtree
## -------------------------------------------------------------------------------------------

library(rpart)
fit.rtree <- rpart(Surv(Time, Status) ~ X1 + X2, data = dat0)
ft <- fit.rtree$cptable
cp <- ft[ft[, 4] == min(ft[, 4]), 1][1]
fit.rtree <- prune(fit.rtree, cp)
r3 <- predict(fit.rtree, newdata = dat.test)
fit.rtree$frame$nd <- 1:dim(fit.rtree$frame)[1]
df.rtree <- merge(data.frame(id = 1:n3, r3 = r3), fit.rtree$frame,
                by.x = "r3", by.y = "yval", sort = FALSE)
df.rtree <- df.rtree[order(df.rtree$id),]
nd3 <- df.rtree$nd


## -------------------------------------------------------------------------------------------
## Fitting partykit::ctree
## -------------------------------------------------------------------------------------------

library(party)
library(partykit)

fit.ctree <- ctree(Surv(Time, Status) ~ X1 + X2, data = dat0)

rn <- 0
if (length(fit.ctree$node) > 1) c3 <- predict(fit.ctree, dat.test, type = "prob") else rm <- 1

SY3c <- SY3r <- SY3tr <- SY3 <- matrix(ncol = length(tt), nrow = n3)
for (i in 1:n3) {
  if (rn == 0) {
    fc3 <- stepfun(c3[[i]]$time, c(1, c3[[i]]$surv))
  } else {
    c3.km <- survfit(Surv(Time, Status) ~ 1, data = dat0)
    fc3 <- stepfun(c3.km$time, c(1, c3.km$surv))
  }
  km <- survfit(Surv(Time, Status) ~ 1, dat0[fit.rtree$where == nd3[i],])
  fr3 <- stepfun(km$time, c(1, km$surv))
  SY3r[i,] <- fr3(tt)
  SY3c[i,] <- fc3(tt)
  SY3[i,] <- sapply(tt, function(u)  exp(-c(sum(dl3[Y <= u, i]))))
}




## 
set.seed(12)
dat <- datGen(10)[,-6]
head(dat)

system.time(foo <- rocTree(Surv(Time, Status) ~ X1 + X2, id = ID, data = dat))

predict(foo, dat)

sim3.1 <- function(n, cen = 0) {
    e <- rbinom(n, 1, .5)
    u <- rexp(n, 5)
    z2 <- runif(n)
    Time <- rep(NA, n)    
    for (i in 1:n) {
        sol <- rexp(1)
        if (e[i] == 1)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u[i]) * .1 * exp(z2[i]) * x -
                (x >= u[i]) * .1 * exp(z2[i]) * (u[i] + exp(1) * (x - u[i])),
                interval = c(0, 50))$root
        if (e[i] == 0)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u[i]) * .1 * exp(z2[i] + 1) * x -
                (x > u[i]) * .1 * exp(z2[i]) * (exp(1) * u[i] + x + u[i]),
                interval = c(0, 50))$root
    }
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 1.73)
    if (cen == .50) cens <- runif(n, 0, 0.83)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = Y[Y <= Y[x]], death = 1 * (Time[x] <= cens[x]), z2 = z2[x], e = e[x], u = u[x])))
    dat$z1 <- with(dat, e * (Y < u) + (1 - e) * (Y >= u))
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "e", "u")])
}

sim3.2 <- function(n, cen = 0) {
    e <- rbinom(n, 1, .5)
    u <- matrix(rexp(3 * n, 5), n)
    u1 <- u[,1]
    u2 <- u[,2]
    u3 <- u[,3]
    z2 <- runif(n)
    Time <- rep(NA, n)
    for (i in 1:n) {
        sol <- rexp(1)
        if (e[i] == 1)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u1[i]) * .1 * exp(z2[i]) * x - 
                (x >= u1[i]) * (x < u2[i]) * .1 * exp(z2[i]) * (u1[i] + exp(1) * (x - u1[i])) -
                (x >= u2[i]) * (x < u3[i]) * .1 * exp(z2[i]) * (u1[i] + exp(1) * (u2[i] - u1[i]) + x - u2[i]) -
                (x >= u3[i]) * .1 * exp(z2[i]) * (u1[i] + u3[i] - u2[i] + exp(1) * (u2[i] - u1[i] + x - u3[i])),
                interval = c(0, 50))$root
        if (e[i] == 0)
            Time[i] <- uniroot(f = function(x)
                sol - (x < u1[i]) * .1 * exp(z2[i]) * exp(1) * x -
                (x >= u1[i]) * (x < u2[i]) * .1 * exp(z2[i]) * (exp(1) * u1[i] + x - u1[i]) - 
                (x >= u2[i]) * (x < u3[i]) * .1 * exp(z2[i]) * (exp(1) * (u1[i] + x - u2[i]) + u2[i] - u1[i]) -
                (x >= u3[i]) * .1 * exp(z2[i]) * (exp(1) * (u1[i] + u3[i] - u2[i]) + u2[i] - u1[i] - u3[i] + x),
                interval = c(0, 50))$root
    }
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 2.11)
    if (cen == .50) cens <- runif(n, 0, 1.02)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = Y[Y <= Y[x]], death = 1 * (Time[x] <= cens[x]), z2 = z2[x], e = e[x], u1 = u1[x], u2 = u2[x], u3 = u3[x])))
    dat$z1 <- with(dat, e * (u1 <= Y) * (Y < u2) + e * (u3 <= Y) + (1 - e) * (Y < u1) + (1 - e) * (u2 <= Y) * (Y < u3))
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "e", "u1", "u2", "u3")])
}

sim3.3 <- function(n, cen = 0) {
    k <- runif(n, 1, 2)
    b <- runif(n, 1, 2)
    z2 <- runif(n)
    Time <- log(10 * rexp(n) * exp(-z2 - b) * k + 1) / k
    if (cen == 0) cens <- rep(Inf, n)
    if (cen == .25) cens <- runif(n, 0, 1.11)
    if (cen == .50) cens <- runif(n, 0, 0.56)
    Y <- pmin(Time, cens)
    dat <- do.call(rbind, lapply(1:n, function(x)
        data.frame(id = x, Y = Y[Y <= Y[x]], death = 1 * (Time[x] <= cens[x]), z2 = z2[x], k = k[x], b = b[x])))
    dat$z1 <- with(dat, k * Time + b)
    return(dat[order(dat$id, dat$Y), c("id", "Y", "death", "z1", "z2", "k", "b")])
}

sim3.1(5, 0)
sim3.1(5, 0.25)
sim3.1(5, 0.5)

sim3.2(5, 0)
sim3.2(5, 0.25)
sim3.2(5, 0.5)

sim3.3(5, 0)
sim3.3(5, 0.25)
sim3.3(5, 0.5)

## Compare with YS's code

#######################################################################
## Load package
#######################################################################

library(rocTree)
library(survival)

#######################################################################
## Mimicking Yifei's draft
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

set.seed(1)
dat <- datGen(200)

## data at baselinef to compare to Yifei's code
dat0 <- dat[cumsum(with(dat, unlist(lapply(split(ID, ID), length), use.names = FALSE))), ]
summary(dat0)
## Checked

foo <- rocTree(Surv(Time, Status) ~ X1 + X2, data = dat, id = ID)

## checking xlist
str(foo$xlist)
summary(c(foo$xlist[[1]]))
summary(c(foo$xlist[[2]]))

str(foo$xlist0)
summary(c(foo$xlist0[[1]]))
summary(c(foo$xlist0[[2]]))
## Checked

## checking new data (testing set)
## constructing new data
set.seed(1)
newdata <- data.frame(ID = rep(1:20, each = length(unique(dat$Time))),
                      Time = sort(unique(dat$Time)),
                      X1 = rep(runif(20), each = length(unique(dat$Time))))
## newdata$X2 <- newdata$Time * seq(1, 2, length = 20) + seq(1, 2, length = 20)
newdata$X2 <- newdata$Time * rep(seq(1, 2, length = 20), each = 200) +
    rep(seq(1, 2, length = 20), each = 200)

unique(newdata$X1)
summary(newdata$X2)
##  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 1.001   1.413   1.717   1.736   2.017   3.267 
## Checked

## checking xlist with X32.list
str(predict(foo, newdata))
tmp <- predict(foo, newdata)
str(tmp$xlist)
summary(c(tmp$xlist[[1]]))
summary(c(tmp$xlist[[2]]))
## checked

set.seed(1)
system.time(foo2 <- rocForest(Surv(Time, Status) ~ X1 + X2, data = dat, id = ID))

set.seed(1)
system.time(foo22 <- rocForest(Surv(Time, Status) ~ X1 + X2, data = dat, id = ID,
                               control = list(parallel = TRUE)))


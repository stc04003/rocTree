library(rocTree)
library(MASS)
source("newSim.R")

do <- function(n, cen, sce) {
    txt1 <- paste("dat <- dCov", sce, "(", n, ",", cen, ")", sep = "")
    eval(parse(text = txt1))
    if (sce == 1) .t0 <- seq(0, 1.53, length = 1000)
    if (sce == 2) .t0 <- seq(0, 30.23, length = 1000)
    if (sce == 3) .t0 <- seq(0, 30.22, length = 1000)
    .Y0 <- dat$Time[cumsum(1:n)]
    fit1 <- rocTree(Surv(Time, death) ~ Z.1 + Z.2 + Z.3 + Z.4 + Z.5 + Z.6 + Z.7 + Z.8 + Z.9 + Z.10,
                    id = ID, data = dat, ensemble = FALSE)
    fit2 <- rocTree(Surv(Time, death) ~ Z.1 + Z.2 + Z.3 + Z.4 + Z.5 + Z.6 + Z.7 + Z.8 + Z.9 + Z.10,
                    id = ID, data = dat)
    err1 <- err2 <- rep(NA, 200)
    for (i in 1:200) {
        ## newX at each .Y0
        if (sce == 1) {
            newX <- dCov1(1, 0)
            trueSv <- 1 - with(newX, ecdf(log(10 * rexp(1e3) * exp(-Z.1 - b) * k + 1)/k)(.t0))
            newX1 <- data.frame(
                Time = .Y0,
                Z.1 = newX$Z.1,
                Z = t(sapply(.Y0, function(yy)
                    mvrnorm(1, rep(yy * newX$k + newX$b, 10),
                            0.5 * 0.9^outer(1:10, 1:10, function(a, b) abs(a - b))))))
        }
        if (sce == 2) {
            newX <- dCov2(1, 0)
            invF <- function(x, k, b, z2, u) {
                10 * log(u) + x - cos(k * x + b + z2)/k + cos(b + z2)/k
            }
            u <- runif(1e3)
            tt <- sapply(u, function(uu)
                with(newX,uniroot(f = invF, interval = c(0,500), k = k, b = b, z2 = Z.1, u = uu)$root))
            trueSv <- 1 - ecdf(tt)(.t0)
            newX1 <- data.frame(
                Time = .Y0,
                Z.1 = newX$Z.1,
                Z = t(sapply(.Y0, function(yy)
                    mvrnorm(1, rep(yy * newX$k + newX$b, 10),
                            0.5 * 0.9^outer(1:10, 1:10, function(a, b) abs(a - b))))))
        }
        if (sce == 3) {
            newX <- dCov3(1, 0)
            invF <- function(x, k, b, z2, u) {
                10 * log(u) + x - cos(k * x * (2 * (x > 5) - 1) + b + z2)/k + cos(b + z2)/k
            }
            u <- runif(1e5)
            tt <- sapply(u, function(uu)
                with(newX,uniroot(f = invF, interval = c(0,500), k = k, b = b, z2 = Z.1, u = uu)$root))
            trueSv <- 1 - ecdf(tt)(.t0)
            newX1 <- data.frame(
                Time = .Y0,
                Z.1 = newX$Z.1, 
                Z = t(sapply(.Y0, function(yy)
                    mvrnorm(1, rep(yy * newX$k * (2 * (yy > 5) - 1) + newX$b, 10),
                            0.5 * 0.9^outer(1:10, 1:10, function(a, b) abs(a - b))))))
        }
        ## newX <- newX[, grep("Z.", names(newX))]
        ref <- stepfun(.t0, c(1, trueSv))(.t0)       
        pred1 <- predict(fit1, newX1)$survFun
        pred2 <- predict(fit2, newX1)$survFun
        err1[i] <- mean(abs(pred1(.t0) - ref))
        err2[i] <- mean(abs(pred2(.t0) - ref))
        ## print(i)
    }
    c(mean(err1), mean(err2))
}

system.time(print(do(200, .25, 1)))
system.time(print(do(200, .25, 2)))
set.seed(0); system.time(print(do(200, .25, 3)))

library(parallel)

cl <- makePSOCKcluster(16)
setDefaultCluster(cl)
invisible(clusterExport(NULL, "do"))
invisible(clusterExport(NULL, "dCov1"))
invisible(clusterExport(NULL, "dCov2"))
invisible(clusterExport(NULL, "dCov3"))
invisible(clusterExport(NULL, "dCov4"))
invisible(clusterExport(NULL, "dCov5"))
invisible(clusterExport(NULL, "ord"))
invisible(clusterExport(NULL, "ordDat"))
invisible(clusterExport(NULL, "trueSurv"))
invisible(clusterExport(NULL, "ti2td"))
invisible(clusterEvalQ(NULL, library(survival)))
invisible(clusterEvalQ(NULL, library(rocTree)))
invisible(clusterEvalQ(NULL, library(MASS)))

sim1.200.025 <- parSapply(NULL, 1:100, function(z) do(200, 0.25, 1))
sim2.200.025 <- parSapply(NULL, 1:100, function(z) do(200, 0.25, 2))
sim3.200.025 <- parSapply(NULL, 1:100, function(z) do(200, 0.25, 3))

stopCluster(cl)


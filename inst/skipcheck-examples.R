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

foo <- do.haz(200, 0, 1.2)


with(foo$pred$pred[[1]], 
     plot(Time, haz, 's', ylim = c(0, max(truehaz1.2(foo$dat)(foo$pred$pred[[1]]$Time)))))
invisible(sapply(1:50, function(e) with(foo$pred$pred[[e]],
                                        lines(Time, haz, 's', col = "lightgray"))))
lines(foo$pred$pred[[1]]$Time, truehaz1.2(foo$dat)(foo$pred$pred[[1]]$Time), 's', col = 2)

with(foo$pred0$pred[[1]],
     plot(Time, haz, 's', ylim = c(0, max(truehaz1.2(foo$dat)(foo$pred$pred[[1]]$Time)))))
invisible(sapply(1:50, function(e) with(foo$pred0$pred[[e]],
                                        lines(Time, haz, 's', col = "lightgray"))))
lines(foo$pred$pred[[1]]$Time, truehaz1.2(foo$dat)(foo$pred$pred[[1]]$Time), 's', col = 2)
invisible(sapply(1:50, function(e) with(foo$pred$pred[[e]],
                                        lines(Time, haz, 's', col = "gray55"))))
legend("topleft", c("True", "Original", "Modified"),
       col = c(2, "lightgray", "gray55"), lty = 1, lwd = 2, bty = "n")

foo$pred$W[[1]][,1]
foo$pred$W[[1]][,50]
diag(foo$pred0$W[[1]])



debug(rocTree:::predict.rocForest)
undebug(rocTree:::predict.rocForest)
foo <- do.haz(20, 0, 1.2)

str(foo$pred$W[[1]])
str(foo$pred0$W[[1]])
summary(c(foo$pred$W[[1]]))
summary(c(foo$pred0$W[[1]]))

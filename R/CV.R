con.cv <- function(y, id, ind, d, x0, beta.seq, parm) {
    yi <- unlist(lapply(split(y, id), max), use.names = FALSE)
    di <- unlist(lapply(split(d, id), max), use.names = FALSE)
    n <- length(unique(id))
    p <- length(x0)
    y1 <- yi[-ind]
    y2 <- yi[ind]
    d1 <- di[-ind]
    d2 <- di[ind]
    x1 <- x2 <- x12 <- list()
    n1 <- length(y1)
    n2 <- length(y2)
    for (i in 1:p) {
        x1[[i]] <- x0[[i]][-ind, -ind]
        x2[[i]] <- x0[[i]][ind, ind]
        x12[[i]] <- x0[[i]][ind, -ind]
        if (!parm@disc[i]) {
            x1tmp <- x1[[i]]
            x2tmp <- x12[[i]]
            x2tmp[rowSums(!is.na(x2tmp)) == 0,] <- 1
            x1[[i]] <- t(apply(x1tmp, 1, rank, ties.method = 'min', na.last = "keep")) / rowSums(!is.na(x1tmp))
            x2[[i]] <- t(sapply(1:n2, function(x) findInterval(x2[[i]][x,], sort(x2tmp[x,])) / sum(!is.na(x2tmp[x,]))))
            x12[[i]] <-t(sapply(1:n2, function(x) findInterval(x12[[i]][x,], sort(x2tmp[x,])) / sum(!is.na(x2tmp[x,]))))
        }
    }
    if (parm@splitBy == "CON") 
        return(CVSeq(y1, d1, x1, y2, d2, x2, x12, beta.seq, parm))
    if (parm@splitBy == "dCON")
        return(CVRP(y1, d1, x1, y2, d2, x2, x12, beta.seq, parm))
}

#' Cross-validation, evaluated when CV = TRUE
#'
#' The resulting tree should be smaller than the full grew tree in rocTree.
#' This is the same procedure with splitBy = "dCON" or "CON".
#'
#' @keywords internal
#' @noRd
CVSeq <- function(Y1, E1, X1.list, Y2, E2, X2.list, X12.list, beta.seq, parm) {
    ## hN1 <- hN2 <- parm@hN
    hN1 <- hN2 <- parm@h * 1.2
    M <- parm@M
    tau <- parm@tau
    minsp <- parm@minsp
    minsp2 <- parm@minsp2
    N1 <- length(X1.list[[1]][1,])
    N2 <- length(X2.list[[1]][1,])
    fitc <- survfit(Surv(c(Y1, Y2), c(1 - E1, 1 - E2)) ~ 1)
    fit <- survfit(Surv(Y1, E1) ~ 1)
    fit$surv <- c(1, fit$surv)
    fitc$surv <- c(1, fitc$surv)
    sc <- fitc$surv[findInterval(Y1, fitc$time)]
    ss <- fit$surv[findInterval(Y1, fit$time)]
    Smat <- matk2 <- outer(ifelse(E1, Y1, NA), Y1, "<=") / matrix(sc, N1, N1, TRUE)
    fmat <- matk <- t(E1 * sapply(ifelse(E1, Y1, NA), K2, vec = Y1, h = hN1) / sc / hN1)
    if (min(rowMeans(fmat)[E1 == 1]) < 0) 
        fmat[E1 == 1 & rowMeans(fmat) < 0, ] <- fmat[which(rowMeans(fmat) > 0)[1], ]
    EE1 <- Y1 <= tau * E1
    EE2 <- Y2 <= tau * E2
    const <- 1 / rowMeans(fmat)[EE1] / sc[EE1] / N1
    sc2 <- fitc$surv[findInterval(Y2, fitc$time)]
    Smat2 <- matk2 <- outer(ifelse(E2, Y2, NA), Y2, "<=") / matrix(sc2, N2, N2, TRUE)
    fmat2 <- matk <- t(E2 * sapply(ifelse(E2, Y2, NA), K2, vec = Y2, h = hN2) / sc2 / hN2)
    if (min(rowMeans(fmat2)[E2 == 1]) < 0) 
        fmat[E2 == 1 & rowMeans(fmat2) < 0, ] <- fmat2[which(rowMeans(fmat2) > 0)[1], ]
    const2 <- 1 / rowMeans(fmat2)[EE2] / sc2[EE2] / N2
    Smat3 <- outer(ifelse(E2, Y2, NA), Y1, "<=") 
    fmat3 <- t(E1 * sapply(ifelse(E2, Y2, NA), K2, vec = Y1, h = hN1) / hN1)
    fTree <- STree <- matrix(NA, M, sum(EE1))
    fTree2 <- STree2 <- fTree3 <- STree3 <- matrix(NA, M, sum(EE2))
    fTree[1, ] <- rowMeans(fmat)[EE1]
    STree[1, ] <- rowMeans(Smat)[EE1]
    fTree2[1, ] <- rowMeans(fmat2)[EE2]
    STree2[1, ] <- rowMeans(Smat2)[EE2]
    fTree3[1, ] <- rowMeans(fmat3)[EE2]
    STree3[1, ] <- rowMeans(Smat3)[EE2]
    treeMat <- data.frame(nd = 1:M, terminal = 0, u = NA, u2 = NA, p = NA, cut = NA)
    treeMat[1, ] <- c(1, 1, 1, 1, NA, NA)
    ndInd1 <- matrix(1, N1, N1)
    ndInd2 <- matrix(1, N2, N2)
    ndInd12 <- matrix(1, N2, N1)
    ndInd1[lower.tri(ndInd1)] <- 0
    ndInd2[lower.tri(ndInd2)] <- 0
    ndInd12[is.na(X12.list[[1]])] <- 0
    tlst <- which(Y1 >= tau)[1] - 1
    while (sum(treeMat$terminal == 1) > 0) {
        sp <- splitSeq(X1.list, Y1, E1, fmat, Smat, treeMat, ndInd1, const, fTree, STree, parm)
        if (sp[1] * 2 < M & !is.na(sp[2])) {
            ndInd1[ndInd1 == sp[1] & X1.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd1[ndInd1 == sp[1] & X1.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            ndInd12[ndInd12 == sp[1] & X12.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd12[ndInd12 == sp[1] & X12.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            fTree[sp[1] * 2, ] <- rowSums(fmat[EE1, diag(ndInd1) == 2 * sp[1], drop = FALSE]) / N1
            fTree[sp[1] * 2 + 1, ] <- rowSums(fmat[EE1, diag(ndInd1) == 2 * sp[1] + 1, drop = FALSE]) / N1
            STree[sp[1] * 2, ] <- rowSums(as.matrix(Smat * (ndInd1 == 2 * sp[1])))[EE1] / N1
            STree[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat * (ndInd1 == 2 * sp[1] + 1)))[EE1] / N1
            fTree2[sp[1] * 2, ] <- rowSums(fmat2[EE2, diag(ndInd2) == 2 * sp[1], drop = FALSE]) / N2
            fTree2[sp[1] * 2 + 1, ] <- rowSums(fmat2[EE2, diag(ndInd2) == 2 * sp[1] + 1, drop = FALSE]) / N2
            STree2[sp[1] * 2, ] <- rowSums(as.matrix(Smat2 * (ndInd2 == 2 * sp[1])))[EE2] / N2
            STree2[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat2 * (ndInd2 == 2 * sp[1] + 1)))[EE2] / N2
            ## rTree2[sp[1]*2,] <- fTree2[sp[1]*2,]/STree2[sp[1]*2,]
            ## rTree2[sp[1]*2+1,] <- fTree2[sp[1]*2+1,]/STree2[sp[1]*2+1,]
            fTree3[sp[1] * 2, ] <- rowSums(fmat3[EE2, diag(ndInd1) == 2 * sp[1], drop = FALSE]) / N1
            fTree3[sp[1] * 2 + 1, ] <- rowSums(fmat3[EE2, diag(ndInd1) == 2 * sp[1] + 1, drop = FALSE]) / N1
            STree3[sp[1] * 2, ] <- rowSums(as.matrix(Smat3 * (ndInd12 == 2 * sp[1])))[EE2] / N1
            STree3[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat3 * (ndInd12 == 2 * sp[1] + 1)))[EE2] / N1
            ## rTree3[sp[1]*2,] <- fTree3[sp[1]*2,]/STree3[sp[1]*2,]
            ## rTree3[sp[1]*2+1,] <- fTree3[sp[1]*2+1,]/STree3[sp[1]*2+1,]
            treeMat[sp[1], 2] <- 0
            treeMat[sp[1], 5:6] <- sp[2:3]
            treeMat[sp[1] * 2, ] <- c(sp[1] * 2, 1,
                                      mean(diag(ndInd1) == sp[1] * 2 & EE1),
                                      min(rowMeans(ndInd1[Y1 <= tau, ] == sp[1] * 2)), NA, NA)
            treeMat[sp[1] * 2 + 1, ] <- c(sp[1] * 2 + 1, 1,
                                          mean(diag(ndInd1) == (sp[1] * 2 + 1) & EE1),
                                          min(rowMeans(ndInd1[Y1 <= tau, ] == sp[1] * 2 + 1)), NA, NA)
            treeMat$terminal[which(treeMat$u <= minsp / N1 & treeMat$u2 <= minsp2 / N1)] <- 2
            ## conTree <- conTree + sp[4]
        } else {
            if (sum(!is.na(treeMat$p)) > 1) break
            else treeMat$terminal[treeMat$nd == sp[1]] <- 2
            break
        }
    }
    if (all(is.na(treeMat$p))) {
        cat("\nNo splits.\n")
        return(NULL)
    }
    ## prune
    treeMatTerm <- treeMat$nd[treeMat$terminal >= 1 & is.na(treeMat$p)]
    ## sort the node
    left <- 2 ^ round(max(log(treeMatTerm, 2)))
    right <- 2 * left - 1
    treeMatTerm2 <- treeMatTerm
    for (i in 1:length(treeMatTerm)) {
        nd <- treeMatTerm[i]
        treeMatTerm2[i] <- nd * 2 ^ ceiling(log(left / nd) / log(2))
    }
    treeMatTerm <- treeMatTerm[order(treeMatTerm2)]
    ## exhaustive searching of all subtrees
    ndTerm <- treeMatTerm
    sizeTree <- length(ndTerm)
    termList <- list()
    indList <- list()
    conList <- list()
    conList2 <- list()
    termList[[1]] <- matrix(ndTerm, nrow = 1)
    indList[[1]] <- matrix(1, nrow = 1, ncol = sizeTree)
    fuTerm <- fTree[ndTerm[order(ndTerm)], ]
    SuTerm <- STree[ndTerm[order(ndTerm)], ]
    fuTerm2 <- fTree2[ndTerm[order(ndTerm)], ]
    SuTerm2 <- STree2[ndTerm[order(ndTerm)], ]
    fuTerm3 <- fTree3[ndTerm[order(ndTerm)], ]
    SuTerm3 <- STree3[ndTerm[order(ndTerm)], ]
    conList[[1]] <- .C("con", 
                       as.integer(length(fuTerm) / sum(EE1)),
                       as.integer(sum(EE1)),
                       as.double(t(fuTerm)),
                       as.double(t(SuTerm)),
                       as.double(const), 
                       double(1))[[6]]
    conList2[[1]] <- .C("con3", 
                        as.integer(length(fuTerm2) / sum(EE2)),
                        as.integer(sum(EE2)),
                        as.double(t(fuTerm2)),
                        as.double(t(SuTerm2)),
                        as.double(t(fuTerm3)),
                        as.double(t(SuTerm3)),
                        as.double(const2),
                        double(1))[[8]]
    ## con3(fuTerm3, SuTerm3, fuTerm2, SuTerm2, const2)
    for (k in 1:(sizeTree - 1)) {
        sizeNew <- sizeTree - k
        mat1 <- termList[[k]]
        mat2 <- indList[[k]]
        num <- 0
        termMat <- matrix(ncol = sizeNew, nrow = 1)
        indMat <- matrix(ncol = sizeNew, nrow = 1)
        consub <- c()
        con2sub <- c()
        for (l in 1:dim(mat1)[1]) {
            ndTerml <- mat1[l, ]
            indl <- mat2[l, ]
            ## change here
            locl <- which(diff(ndTerml) == 1 & ndTerml[-length(ndTerml)] %% 2 == 0 &
                          indl[-length(ndTerml)] == 1)
            for (q in locl) {
                ndNew <- ndTerml
                ndNew[q] <- ndNew[q] / 2
                ndNew <- ndNew[-(q + 1)]
                if (l == 1 & q == locl[1]) indNew <- rep(1, sizeNew)
                else indNew <- rep(0:1, c(q - 1, sizeNew - q + 1))
                termMat <- rbind(termMat, ndNew)
                indMat <- rbind(indMat, indNew)
                num <- num + 1
                fuTerm <- fTree[treeMat$nd %in% ndNew, ]
                SuTerm <- STree[treeMat$nd %in% ndNew, ]
                fuTerm2 <- fTree2[treeMat$nd %in% ndNew, ]
                SuTerm2 <- STree2[treeMat$nd %in% ndNew, ]
                fuTerm3 <- fTree3[treeMat$nd %in% ndNew, ]
                SuTerm3 <- STree3[treeMat$nd %in% ndNew, ]
                consub <- c(consub, .C("con", 
                                       as.integer(length(fuTerm) / sum(EE1)),
                                       as.integer(sum(EE1)),
                                       as.double(t(fuTerm)),
                                       as.double(t(SuTerm)),
                                       as.double(const), 
                                       double(1))[[6]])
                con2sub <- c(con2sub, .C("con3", 
                                         as.integer(length(fuTerm2) / sum(EE2)),
                                         as.integer(sum(EE2)),
                                         as.double(t(fuTerm2)),
                                         as.double(t(SuTerm2)),
                                         as.double(t(fuTerm3)),
                                         as.double(t(SuTerm3)),
                                         as.double(const2),
                                         double(1))[[8]])
            }
        }
        termMat <- termMat[-1, ]
        indMat <- indMat[-1, ]
        if (is.vector(termMat)) {
            termList[[k + 1]] <- matrix(termMat, nrow = 1)
            indList[[k + 1]] <- matrix(indMat, nrow = 1)
        } else {
            termList[[k + 1]] <- termMat
            indList[[k + 1]] <- indMat
        }
        conList[[k + 1]] <- consub
        conList2[[k + 1]] <- con2sub
    }
    con1 <- unlist(lapply(conList, max))
    sz1 <- sizeTree:1
    j <- 1
    res <- matrix(c(1, 0), nrow = 1)
    while (j < sizeTree) {
        alpha <- (con1[j] - con1[-(1:j)]) / (1:(sizeTree - j))
        j <- j + tail(which(alpha == min(alpha)), 1)
        res <- rbind(res, c(j, min(alpha)))
    }
    con2opt <- (unlist(lapply(conList2, max)))[res[, 1]]
    con1opt <- con1[res[, 1]]
    res <- cbind(res, con1opt, con2opt)
    res <- rbind(res, c(100, Inf, 0, 0))
    cvconbb <- function(bb) res[max(1, which(res[, 2] >= bb)[1] - 1), 4]
    cvcon2 <- sapply(beta.seq[-length(beta.seq)], cvconbb)
    c(unlist(cvcon2), tail(res[, 4], 2)[1])
}


CVRP <- function(Y1, E1, X1.list, Y2, E2, X2.list, X12.list, beta.seq, parm) {
    ## hN1 <- hN2 <- parm@hN
    hN1 <- hN2 <- parm@h * 1.2
    M <- parm@M
    tau <- parm@tau
    minsp <- parm@minsp
    minsp2 <- parm@minsp2
    N1 <- length(X1.list[[1]][1,])
    N2 <- length(X2.list[[1]][1,])
    fitc <- survfit(Surv(c(Y1, Y2), c(1 - E1, 1 - E2)) ~ 1)
    fit <- survfit(Surv(Y1, E1) ~ 1)
    fit$surv <- c(1, fit$surv)
    fitc$surv <- c(1, fitc$surv)
    sc <- fitc$surv[findInterval(Y1, fitc$time)]
    ss <- fit$surv[findInterval(Y1, fit$time)]
    Smat <- matk2 <- outer(ifelse(E1, Y1, NA), Y1, "<=") / matrix(sc, N1, N1, TRUE)
    fmat <- matk <- t(E1 * sapply(ifelse(E1, Y1, NA), K2, vec = Y1, h = hN1) / sc / hN1)
    if (min(rowMeans(fmat)[E1 == 1]) < 0) 
        fmat[E1 == 1 & rowMeans(fmat) < 0, ] <- fmat[which(rowMeans(fmat) > 0)[1], ]
    EE1 <- Y1 <= tau * E1
    EE2 <- Y2 <= tau * E2
    const <- 1 / rowMeans(fmat)[EE1] / sc[EE1] / N1
    sc2 <- fitc$surv[findInterval(Y2, fitc$time)]
    Smat2 <- matk2 <- outer(ifelse(E2, Y2, NA), Y2, "<=") / matrix(sc2, N2, N2, TRUE)
    fmat2 <- matk <- t(E2 * sapply(ifelse(E2, Y2, NA), K2, vec = Y2, h = hN2) / sc2 / hN2)
    if (min(rowMeans(fmat2)[E2 == 1]) < 0) 
        fmat[E2 == 1 & rowMeans(fmat2) < 0, ] <- fmat2[which(rowMeans(fmat2) > 0)[1], ]
    const2 <- 1 / rowMeans(fmat2)[EE2] / sc2[EE2] / N2
    Smat3 <- outer(ifelse(E2, Y2, NA), Y1, "<=") 
    fmat3 <- t(E1 * sapply(ifelse(E2, Y2, NA), K2, vec = Y1, h = hN1) / hN1)
    fTree <- STree <- matrix(NA, M, sum(EE1))
    fTree2 <- STree2 <- fTree3 <- STree3 <- matrix(NA, M, sum(EE2))
    fTree[1, ] <- rowMeans(fmat)[EE1]
    STree[1, ] <- rowMeans(Smat)[EE1]
    fTree2[1, ] <- rowMeans(fmat2)[EE2]
    STree2[1, ] <- rowMeans(Smat2)[EE2]
    fTree3[1, ] <- rowMeans(fmat3)[EE2]
    STree3[1, ] <- rowMeans(Smat3)[EE2]
    treeMat <- data.frame(nd = 1:M, terminal = 0, u = NA, u2 = NA, p = NA, cut = NA)
    treeMat[1, ] <- c(1, 1, 1, 1, NA, NA)
    ndInd1 <- matrix(1, N1, N1)
    ndInd2 <- matrix(1, N2, N2)
    ndInd12 <- matrix(1, N2, N1)
    ndInd1[lower.tri(ndInd1)] <- 0
    ndInd2[lower.tri(ndInd2)] <- 0
    ndInd12[is.na(X12.list[[1]])] <- 0
    tlst <- which(Y1 >= tau)[1] - 1
    while (sum(treeMat$terminal == 1) > 0) {
        sp <- splitRP(X1.list, Y1, E1, fmat, Smat, treeMat, ndInd1, const, fTree, STree, parm)
        if (sp[1] * 2 < M & !is.na(sp[2])) {
            ndInd1[ndInd1 == sp[1] & X1.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd1[ndInd1 == sp[1] & X1.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            ndInd12[ndInd12 == sp[1] & X12.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
            ndInd12[ndInd12 == sp[1] & X12.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
            fTree[sp[1] * 2, ] <- rowSums(fmat[EE1, diag(ndInd1) == 2 * sp[1], drop = FALSE]) / N1
            fTree[sp[1] * 2 + 1, ] <- rowSums(fmat[EE1, diag(ndInd1) == 2 * sp[1] + 1, drop = FALSE]) / N1
            STree[sp[1] * 2, ] <- rowSums(as.matrix(Smat * (ndInd1 == 2 * sp[1])))[EE1] / N1
            STree[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat * (ndInd1 == 2 * sp[1] + 1)))[EE1] / N1
            fTree2[sp[1] * 2, ] <- rowSums(fmat2[EE2, diag(ndInd2) == 2 * sp[1], drop = FALSE]) / N2
            fTree2[sp[1] * 2 + 1, ] <- rowSums(fmat2[EE2, diag(ndInd2) == 2 * sp[1] + 1, drop = FALSE]) / N2
            STree2[sp[1] * 2, ] <- rowSums(as.matrix(Smat2 * (ndInd2 == 2 * sp[1])))[EE2] / N2
            STree2[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat2 * (ndInd2 == 2 * sp[1] + 1)))[EE2] / N2
            ## rTree2[sp[1]*2,] <- fTree2[sp[1]*2,]/STree2[sp[1]*2,]
            ## rTree2[sp[1]*2+1,] <- fTree2[sp[1]*2+1,]/STree2[sp[1]*2+1,]
            fTree3[sp[1] * 2, ] <- rowSums(fmat3[EE2, diag(ndInd1) == 2 * sp[1], drop = FALSE]) / N1
            fTree3[sp[1] * 2 + 1, ] <- rowSums(fmat3[EE2, diag(ndInd1) == 2 * sp[1] + 1, drop = FALSE]) / N1
            STree3[sp[1] * 2, ] <- rowSums(as.matrix(Smat3 * (ndInd12 == 2 * sp[1])))[EE2] / N1
            STree3[sp[1] * 2 + 1, ] <- rowSums(as.matrix(Smat3 * (ndInd12 == 2 * sp[1] + 1)))[EE2] / N1
            ## rTree3[sp[1]*2,] <- fTree3[sp[1]*2,]/STree3[sp[1]*2,]
            ## rTree3[sp[1]*2+1,] <- fTree3[sp[1]*2+1,]/STree3[sp[1]*2+1,]
            treeMat[sp[1], 2] <- 0
            treeMat[sp[1], 5:6] <- sp[2:3]
            treeMat[sp[1] * 2, ] <- c(sp[1] * 2, 1,
                                      mean(diag(ndInd1) == sp[1] * 2 & EE1),
                                      min(rowMeans(ndInd1[Y1 <= tau, ] == sp[1] * 2)), NA, NA)
            treeMat[sp[1] * 2 + 1, ] <- c(sp[1] * 2 + 1, 1,
                                          mean(diag(ndInd1) == (sp[1] * 2 + 1) & EE1),
                                          min(rowMeans(ndInd1[Y1 <= tau, ] == sp[1] * 2 + 1)), NA, NA)
            treeMat$terminal[which(treeMat$u <= minsp / N1 & treeMat$u2 <= minsp2 / N1)] <- 2
            ## conTree <- conTree + sp[4]
        } else {
            if (sum(!is.na(treeMat$p)) > 1) break
            else treeMat$terminal[treeMat$nd == sp[1]] <- 2
            break
        }
    }
    if (all(is.na(treeMat$p))) {
        cat("\nNo splits.\n")
        return(NULL)
    }
    ## prune
    treeMatTerm <- treeMat$nd[treeMat$terminal >= 1 & is.na(treeMat$p)]
    ## sort the node
    left <- 2 ^ round(max(log(treeMatTerm, 2)))
    right <- 2 * left - 1
    treeMatTerm2 <- treeMatTerm
    for (i in 1:length(treeMatTerm)) {
        nd <- treeMatTerm[i]
        treeMatTerm2[i] <- nd * 2 ^ ceiling(log(left / nd) / log(2))
    }
    treeMatTerm <- treeMatTerm[order(treeMatTerm2)]
    ## exhaustive searching of all subtrees
    ndTerm <- treeMatTerm
    sizeTree <- length(ndTerm)
    termList <- list()
    indList <- list()
    conList <- list()
    conList2 <- list()
    termList[[1]] <- matrix(ndTerm, nrow = 1)
    indList[[1]] <- matrix(1, nrow = 1, ncol = sizeTree)
    fuTerm <- fTree[ndTerm[order(ndTerm)], ]
    SuTerm <- STree[ndTerm[order(ndTerm)], ]
    fuTerm2 <- fTree2[ndTerm[order(ndTerm)], ]
    SuTerm2 <- STree2[ndTerm[order(ndTerm)], ]
    fuTerm3 <- fTree3[ndTerm[order(ndTerm)], ]
    SuTerm3 <- STree3[ndTerm[order(ndTerm)], ]
    conList[[1]] <- .C("con", 
                       as.integer(length(fuTerm) / sum(EE1)),
                       as.integer(sum(EE1)),
                       as.double(t(fuTerm)),
                       as.double(t(SuTerm)),
                       as.double(const), 
                       double(1))[[6]]
    conList2[[1]] <- .C("con3", 
                        as.integer(length(fuTerm2) / sum(EE2)),
                        as.integer(sum(EE2)),
                        as.double(t(fuTerm2)),
                        as.double(t(SuTerm2)),
                        as.double(t(fuTerm3)),
                        as.double(t(SuTerm3)),
                        as.double(const2),
                        double(1))[[8]]
    ## con3(fuTerm3, SuTerm3, fuTerm2, SuTerm2, const2)
    for (k in 1:(sizeTree - 1)) {
        sizeNew <- sizeTree - k
        mat1 <- termList[[k]]
        mat2 <- indList[[k]]
        num <- 0
        termMat <- matrix(ncol = sizeNew, nrow = 1)
        indMat <- matrix(ncol = sizeNew, nrow = 1)
        consub <- c()
        con2sub <- c()
        for (l in 1:dim(mat1)[1]) {
            ndTerml <- mat1[l, ]
            indl <- mat2[l, ]
            ## change here
            locl <- which(diff(ndTerml) == 1 & ndTerml[-length(ndTerml)] %% 2 == 0 &
                          indl[-length(ndTerml)] == 1)
            for (q in locl) {
                ndNew <- ndTerml
                ndNew[q] <- ndNew[q] / 2
                ndNew <- ndNew[-(q + 1)]
                if (l == 1 & q == locl[1]) indNew <- rep(1, sizeNew)
                else indNew <- rep(0:1, c(q - 1, sizeNew - q + 1))
                termMat <- rbind(termMat, ndNew)
                indMat <- rbind(indMat, indNew)
                num <- num + 1
                fuTerm <- fTree[treeMat$nd %in% ndNew, ]
                SuTerm <- STree[treeMat$nd %in% ndNew, ]
                fuTerm2 <- fTree2[treeMat$nd %in% ndNew, ]
                SuTerm2 <- STree2[treeMat$nd %in% ndNew, ]
                fuTerm3 <- fTree3[treeMat$nd %in% ndNew, ]
                SuTerm3 <- STree3[treeMat$nd %in% ndNew, ]
                consub <- c(consub, .C("con", 
                                       as.integer(length(fuTerm) / sum(EE1)),
                                       as.integer(sum(EE1)),
                                       as.double(t(fuTerm)),
                                       as.double(t(SuTerm)),
                                       as.double(const), 
                                       double(1))[[6]])
                con2sub <- c(con2sub, .C("con3", 
                                         as.integer(length(fuTerm2) / sum(EE2)),
                                         as.integer(sum(EE2)),
                                         as.double(t(fuTerm2)),
                                         as.double(t(SuTerm2)),
                                         as.double(t(fuTerm3)),
                                         as.double(t(SuTerm3)),
                                         as.double(const2),
                                         double(1))[[8]])
            }
        }
        termMat <- termMat[-1, ]
        indMat <- indMat[-1, ]
        if (is.vector(termMat)) {
            termList[[k + 1]] <- matrix(termMat, nrow = 1)
            indList[[k + 1]] <- matrix(indMat, nrow = 1)
        } else {
            termList[[k + 1]] <- termMat
            indList[[k + 1]] <- indMat
        }
        conList[[k + 1]] <- consub
        conList2[[k + 1]] <- con2sub
    }
    con1 <- unlist(lapply(conList, max))
    sz1 <- sizeTree:1
    j <- 1
    res <- matrix(c(1, 0), nrow = 1)
    while (j < sizeTree) {
        alpha <- (con1[j] - con1[-(1:j)]) / (1:(sizeTree - j))
        j <- j + tail(which(alpha == min(alpha)), 1)
        res <- rbind(res, c(j, min(alpha)))
    }
    con2opt <- (unlist(lapply(conList2, max)))[res[, 1]]
    con1opt <- con1[res[, 1]]
    res <- cbind(res, con1opt, con2opt)
    res <- rbind(res, c(100, Inf, 0, 0))
    cvconbb <- function(bb) res[max(1, which(res[, 2] >= bb)[1] - 1), 4]
    cvcon2 <- sapply(beta.seq[-length(beta.seq)], cvconbb)
    c(unlist(cvcon2), tail(res[, 4], 2)[1])
}

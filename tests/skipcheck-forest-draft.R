set.seed(1)
N <- 200

P <- 2
disc <- c(0,0)
U <- runif(N)
X <- matrix(runif(N*(P),0,1), ncol = P, nrow = N)
a <- X[,1]
b <- X[,2]
betat <- 2
beta <- 1
lambda <- .1
ks <- runif(N,1,2)
it <- runif(N,1,2)
T <- log(1+(betat*ks*(-log(U)))/lambda/exp(beta*a+betat*it))/betat/ks  
C <- runif(N, 0, .8)
tau <- .4
Y <- pmin(T,C)
E <- T<=C
  
X <- X[order(Y),]
a <- a[order(Y)]
b <- b[order(Y)]
E <- E[order(Y)]
ks <- ks[order(Y)]
it <- it[order(Y)]
Y <- Y[order(Y)]

X.list = list()
cdf0.list = list()
X0.list = list()
P1 <- 1
for(p in 1:P1) {
    X.list[[p]] <- matrix(rep(X[,p],each = N), ncol = N, nrow = N)
    X.list[[p]][lower.tri(X.list[[p]], diag = FALSE)] <- NA
    X0.list[[p]] <- X.list[[p]]
    if(disc[p] == 0){
        for(j in 1:N)
        {
            x <- X.list[[p]][j,j:N]
            cdf0.list[[(p-1)*N+j]] <- ecdf(x)
            X.list[[p]][j,j:N] <- (ecdf(x))(x)
        }
    }
}
P2 <- 2
for(p in (P1+1):P2)
{
    X.list[[p]] <- matrix(NA, ncol = N, nrow = N)
    X.list[[p]] <- matrix(rep(Y,N), ncol = N, nrow = N)*matrix(rep(ks,each = N), ncol = N, nrow = N)+matrix(rep(it,each = N), ncol = N, nrow = N)
    X.list[[p]][lower.tri(X.list[[p]], diag = FALSE)] <- NA
    X0.list[[p]] <- X.list[[p]]
    if(disc[p] == 0)
    {
        for(j in 1:N)
        {
            x <- X.list[[p]][j,j:N]
            cdf0.list[[(p-1)*N+j]] <- ecdf(x)
            X.list[[p]][j,j:N] <- (ecdf(x))(x)
        }
    }
}


## New data (testing set)

set.seed(1)
N3 <- 20
X3 <- matrix(runif(N3*(P),0,1), ncol = P, nrow = N3)
a3 <- X3[,1]
b3 <- X3[,2]
k3 <- seq(from = 1, to = 2, length = N3)#runif(N3,1,2)
it3 <- seq(from = 1, to = 2, length = N3)#runif(N3,1,2)
X32.list = list()

for(p in 1:P1) {
    X32.list[[p]] <- matrix(rep(X3[,p],each = N), ncol = N3, nrow = N)
    if(disc[p] == 0) {
        for(j in 1:N) {
            X32.list[[p]][j,] <- (cdf0.list[[(p-1)*N+j]])(X32.list[[p]][j,])
        }
    }
}
for(p in (P1+1):P2)
{
    X32.list[[p]] <- matrix(NA, ncol = N3, nrow = N)
    X32.list[[p]] <- matrix(rep(Y,N3), ncol = N3, nrow = N)*matrix(rep(k3,each = N), ncol = N3, nrow = N)+matrix(rep(it3,each = N), ncol = N3, nrow = N)
    if(disc[p] == 0)
    {
        for(j in 1:N)
        {
            X32.list[[p]][j,] <- (cdf0.list[[(p-1)*N+j]])(X32.list[[p]][j,])
        }
    }
}


control <- list(minsp = 6, minsp2 = 6,
                M = 500, tau = tau, disc = c(0,0),
                Trace = FALSE, hN = tau/15)
#################################################################
#                       some functions                          #
#################################################################
# This is a new function
# Y1,E1,X1.list is from L1
# X2.list is from L2
# X32.list is the covariates of test observations
# Y is the Y of the training sample 
# idb2 is the ids of L2
grow4 <- function(Y1, E1, X1.list, X2.list, X32.list, Y, idb2, control) {
  tau <- control$tau
  M <- control$M
  hN <- control$hN
  minsp <- control$minsp
  minsp2 <- control$minsp2
  
  N1 <- length(X1.list[[1]][1, ])
  fit <- survfit(Surv(Y1, E1) ~ 1)
  fit$surv <- c(1, fit$surv)
  ss <- fit$surv[findInterval(Y1, fit$time)]
  Smat <- matk2 <- outer(ifelse(E1, Y1, NA), Y1, "<=") ## opposite of Y >= Y[i]
  fmat <- matk <- t(E1 * sapply(ifelse(E1, Y1, NA), K2, vec = Y1, h = hN) / hN)
  if (min(rowMeans(fmat)[E1 == 1]) < 0) 
    fmat[E1 == 1 & rowMeans(fmat) < 0, ] <- fmat[which(rowMeans(fmat) > 0)[1], ]
  Stau <- fit$surv[findInterval(tau, fit$time) + 1]
  ss <- ss[Y1 <= tau & E1 == 1]
  fall <- rowMeans(fmat)
  ## const is actually dF/f
  const <- -diff(c(ss, Stau)) / (rowMeans(fmat)[Y1 <= tau * E1])
  ## initialization
  fTree <- STree <- matrix(NA, M, sum(Y1 <= E1 * tau))
  fTree[1, ] <- rowMeans(fmat)[Y1 <= tau * E1]
  STree[1, ] <- rowMeans(Smat)[Y1 <= tau * E1]
  ## Define a tree object
  treeMat <- matrix(nrow = M, ncol = 6)
  colnames(treeMat) <- c("nd", "terminal", "u", "u2", 
                         "p", "cut")
  ## u3 is u2 in L2
  ## terminal = 0 - internal
  ## terminal = 1 - terminal can be split
  ## terminal = 2 - terminal cannot be split
  treeMat[, 1] <- 1:M
  treeMat[, 2] <- 0
  treeMat[1, ] <- c(1, 1, 1, 1, NA, NA)
  ## node number of each observation
  ## each col is one subject
  ndInd <- matrix(1, N1, N1)
  ndInd[lower.tri(ndInd)] <- 0
  
  N2 <- dim(X2.list[[1]])[2]
  ndInd2 <- matrix(1, N, N2)
  ndInd2[outer(Y,Y2,FUN=">")] <- 0
  
  conTree <- sum(0.5 * const * ss * rowMeans(fmat)[Y1 <= tau * E1])
  
  while (sum(treeMat[, 2] == 1) > 0) {
### STEVEN: split function need to be modified so that m feature is selected at each split
    sp <- split3(X1.list, Y1, E1, fmat, Smat, treeMat, ndInd, const, fTree, STree, control)
    if (sp[1] * 2 < M & !is.na(sp[2])) {
      ndInd[ndInd == sp[1] & X1.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
      ndInd[ndInd == sp[1] & X1.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
      
      ## update L2 simutaneously
      ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] <= sp[3]] <- sp[1] * 2
      ndInd2[ndInd2 == sp[1] & X2.list[[sp[2]]] > sp[3]] <- sp[1] * 2 + 1
      
      treeMat[sp[1], 2] <- 0
      treeMat[sp[1], 5:6] <- sp[2:3]
      treeMat[sp[1] * 2, ] <- c(sp[1] * 2, 1,
                                mean(diag(ndInd) == sp[1] * 2 & Y1 <= tau * E1),
                                min(rowMeans(ndInd[Y1 <= tau, ] == sp[1] * 2)), 
                                NA, NA)
      treeMat[sp[1] * 2 + 1, ] <- c(sp[1] * 2 + 1, 1,
                                    mean(diag(ndInd) == sp[1] * 2 + 1 & Y1 <= tau * E1),
                                    min(rowMeans(ndInd[Y1 <= tau, ] == sp[1] * 2 + 1)), 
                                    NA, NA)
      ## mean(ndInd[tlst,] == sp[1]*2)
      treeMat[treeMat[, 3] < minsp / N1 & treeMat[, 4] < minsp2 / N1, 2] <- 2
      
      
      if(sum(diag(ndInd) == 2 * sp[1])>1)
      {
        fTree[sp[1] * 2, ] <- rowSums(fmat[Y1 <= tau * E1, diag(ndInd) == 2 * sp[1]]) / N1
        
      }else{
        fTree[sp[1] * 2, ] <- sum(fmat[Y1 <= tau * E1, diag(ndInd) == 2 * sp[1]]) / N1
        
      }
      if(sum(diag(ndInd) == 2 * sp[1]+1)>1)
      {
        fTree[sp[1] * 2 + 1, ] <- rowSums(fmat[Y1 <= tau * E1, diag(ndInd) == 2 * sp[1] + 1]) / N1
      }else{
        fTree[sp[1] * 2 + 1, ] <- sum(fmat[Y1 <= tau * E1, diag(ndInd) == 2 * sp[1] + 1]) / N1
      }
      STree[sp[1] * 2, ] <- rowSums(Smat * (ndInd == 2 * sp[1]))[Y1 <= tau * E1] / N1
      STree[sp[1] * 2 + 1, ] <- rowSums(Smat * (ndInd == 2 * sp[1] + 1))[Y1 <= tau * E1] / N1
      
      conTree <- conTree + sp[4]
    } else {
      treeMat[treeMat[, 1] == sp[1], 2] <- 2
      break
    }
    if (control$Trace) print(sp)
  }
  # Terminal nodes 
  ndTerm <- treeMat[treeMat[,2]==2,1]
  
  # SzL2 is a matrix for L2, counting the size of terminal nodes at each time point
  # each col is for one terminal node
  # each row is for one time point
  SzL2 <- matrix(ncol = length(ndTerm), nrow = N)
  for(k in 1:length(ndTerm))
  {
    SzL2[,k] <- apply(ndInd2, 1, function(x){sum(x==ndTerm[k])})
  }

  tmp <- sapply(ndTerm, function(x) rowSums(ndInd2 == x))
  
  # final tree
  tree2 <- treeMat[!is.na(treeMat[,3]),]
  
  # find the node label of the new observations in the tree
  ndInd32 <- matrix(1, N, N3)
  for(i in 1:dim(tree2)[1])
  {
    sp <- tree2[i,c(1,5:6,2)]
    if(sp[4] == 0)
    {
      ndInd32[ndInd32 == sp[1] & X32.list[[sp[2]]]<=sp[3]] <- sp[1]*2
      ndInd32[ndInd32 == sp[1] & X32.list[[sp[2]]]>sp[3]] <- sp[1]*2+1
    }else{
      # If it is a terminal node, do nothing
      # Nodes who are descendents of the terminal nodes will not appear
    }
  }
  
  # W2 is an 3d array, W2[i,,] is the weight matrix for the ith new obs.
  # each row is one subject
  W2 <- array(0, dim = c(N3,N,N))
  for(i in 1:N3)
  {
    ndi <- ndInd32[,i]
    # the nodes over time of the ith subjects (new obs.)
    for(j in 1:N)
    {
      # at jth time point, subject i enters node ndi[j]
      # update the weights for observations in L2 that falls into ndi[j]
      # If the observation falls into a node that no data from L2 is in that node, this 
      # means that there is no neighbor in L2, so the weight is 0 (later we do normalization so that the sum is 1)
      W2[i,j,idb2[ndInd2[j,] == ndi[j]]] <- 1/SzL2[j,ndTerm == ndi[j]]
    }
  }
  list(treeMat = tree2, W2 = W2)
}


# I think I did not change this function
split3 <- function(X, Y, E, fmat, Smat, treeMat, ndInd, const, fTree, STree, control)  {
  N <- dim(X[[1]])[1]
  P <- length(X)
  disc <- control$disc
  minsp <- control$minsp
  minsp2 <- control$minsp2
  M <- control$M
  tau <- control$tau
  ## all the terminal nodes that can be split; not all the terminal nodes
  nd.terminal <- treeMat[treeMat[, 2] == 1, 1]
  sopt <- matrix(NA, M, 2)
  dconopt <- rep(0, M)
  ## lnd <- length(nd.terminal)-1
  lt <- sum(Y <= tau * E)
  for (m in nd.terminal) {
    ## need to change with discrete data
    dconList <- list()
    cutList <- list()
    f <- fTree[(treeMat[, 2] >= 1) & (treeMat[, 1] != m), ]
    S <- STree[(treeMat[, 2] >= 1) & (treeMat[, 1] != m), ]
    ## size of nodes
    fm <- fTree[m, ]
    Sm <- STree[m, ]
    rm <- fm / Sm
    r <- f / S
    r[is.na(r)] <- Inf
    rm[is.na(r)] <- Inf
    for (p in 1:P) {
      ## if discrete
      if (disc[p] == 0) {
        cutAll <- sort(unique(X[[p]][1, ndInd[1, ] == m]))
      } else {
        cutAll <- sort(unique(X[[p]][ndInd == m]))
      }
      cutAll <- cutAll[-length(cutAll)]
      ## if there is no potential cut off, skip
      if (length(cutAll) == 0) {
        dconList[[p]] <- -1
        next
      }
      cutList[[p]] <- cutAll
      dconList[[p]] <- .C("cutSearch",
                          as.integer(N), ## n
                          as.integer(length(cutAll)), ## cL
                          as.integer(m), ## m
                          as.integer(which(Y <= tau * E) - 1), ## y
                          as.integer(sum(Y <= tau * E)), ## Ny
                          as.double(minsp), ## minsp
                          as.double(minsp2), ## minsp2
                          as.double(ifelse(is.na(X[[p]]), 0, X[[p]])), ## X
                          as.double(ndInd), ## ndInd
                          as.double(cutAll), ## cut
                          as.double(ifelse(is.na(fmat), 0, fmat)), ## fmat
                          as.double(ifelse(is.na(Smat), 0, Smat)), ## Smat
                          as.double(const), ## const cL by 1
                          as.double(t(f)), ## length(nd.terminal) by Ny
                          as.double(t(S)), ## length(nd.terminal) by Ny
                          as.double(t(ifelse(r == Inf, 99999, r))), ## length(nd.terminal) by Ny
                          as.double(ifelse(rm == Inf | is.na(rm), 99999, rm)), ## Ny by 1
                          as.integer(sum(treeMat[,2] >= 1 & treeMat[,1] != m)), ## sum(treeMat[,2] >= 1 & treeMat[,1] != m)
                          out = as.double(double(length(cutAll))), PACKAGE = "rocTree")$out
    } ## end P
    dconmaxP <- unlist(lapply(dconList, max))
    if (max(dconmaxP) < 0) {
      sopt[m, ] <- c(P + 1, 999)
      dconopt[m] <- -1
    } else {
      dconopt[m] <- max(dconmaxP)
      indm <- which.max(dconmaxP)
      cutindm <- cutList[[indm[1]]]
      dconindm <- dconList[[indm[1]]]
      sopt[m, ] <- c(indm[1], cutindm[which.max(dconindm)])
    }
  }
  nd.split <- which.max(dconopt)
  sopt2 <- sopt[nd.split, ]
  dconopt2 <- max(dconopt)
  c(nd.split, sopt2, dconopt2)
  ## which node to split, which variable, which cut off, which dcon
}



K1 <- function(u) {
  0.75 * (1 - u ^ 2) * (abs(u) < 1)
}

K2 <- function(s, vec, h) {
  if (is.na(s)) return(rep(NA, length(vec)))
  if (s < h) return(Kq((s - vec) / h, s / h))
  else return(K1((s - vec) / h))
}

Kq <- function(x, q) {
  sigk1 <- sqrt(0.2)
  2 / (q + 1) * K1(2 / (q + 1) * (x - (q - 1) / 2)) *
    (1 + ((q - 1) / (q + 1) / sigk1) ^ 2 + 2 / sigk1 ^ 2 * (1 - q) / (1 + q) ^ 2 * x)
}


#################################################################
#                           forests                             #
#################################################################

  matk <- matrix(ncol = N, nrow = N)
  matk2 <- matrix(ncol = N, nrow = N)
  matk3 <- matrix(ncol = N, nrow = N)
  # each row is one time point
  for(i in 1:N)
  {
    matk2[i,] <- (Y>=Y[i])
    matk3[i,] <- E*(Y == Y[i])
  }
  
  # weight matrix for prediction of new N3 observations
  # We have a weight matrix for each new observation, each col is one subject in the training sample, each row is a time point
  WW <- array(0, dim = c(N3,N,N))
  # resampling
B <- 500

  for(b in 1:B) {
    # subsampling
    S <- round(N*0.5)
    if(S%%2 != 0) {S <- S+1}
    idb <- sample(1:N, size = S)
    # divide the data into two halves
    S1 <- round(S/2)
    id1 <- sample(1:S, size = S1)
    
    # idb1 is the id in 1st half
    idb1 <- sort(idb[id1])
    idb2 <- sort(idb[-id1])
    
    Y1 <- Y[idb1]
    E1 <- E[idb1]
    Y2 <- Y[idb2]
    E2 <- E[idb2]
    # X1.list is the X of L1 on Y1
    # X2.list is the X of L2 on Y
    # each row is a time point, each col is a subject
    X1.list <- list()
    X2.list <- list()
    for(p in 1:P)
    {
      Xp <- X.list[[p]]
      X1p <- Xp[idb1,idb1]
      X2p <- Xp[,idb2]
      X1.list[[p]] <- X1p
      X2.list[[p]] <- X2p
    }
    fit4 <- grow4(Y1, E1, X1.list, X2.list, X32.list, Y, idb2, control)
      ## fit4 <- try(grow4(Y1, E1, X1.list, X2.list, X32.list, Y, idb2, control), silent = TRUE)
    ### STEVEN: Sometimes I got error messages here, which is from the splitting function
    ### It might be that the node size is too small for valid calculation. Could you check what went wrong? Thanks!!!
    if(is.list(fit4))
    {WW <- WW + fit4$W2}
    else{print(paste("error in the ", b, "th subsample, not used for weights", sep = ""))}
  }
  
## check the results
  
  
  Ltrue3 <- function(tt,k0,it0,a0)
  {
    lambda*exp(beta*a0+betat*it0)*(exp(betat*k0*tt)-1)/(betat*k0)
  }
  
  par(mfrow = c(4,5))
  tt <- Y
  # prediction for the 20 new observations
  for(i in 1:N3)
  {
    # normalize the weight so the row sum is 1
    Wi <- WW[i,,]/rowSums(WW[i,,])
    # Spred is our predicted survival function on each time point in Y
    Spred <- exp(-cumsum(rowSums(matk3*Wi)/rowSums(matk2*Wi)))
    # rowSums(matk2*Wi) is actually 1
    Strue <- exp(-Ltrue3(tt,k3[i],it3[i],a3[i]))
    plot(tt[tt<tau],Strue[tt<tau], ylim = c(0,1), type = "l")
    points(tt,Spred, col = 2, type = "l")
  }
  
  # check the weight matrix for i th new observation
  i <- 18
  Wi <- WW[i,,]/rowSums(WW[i,,])
  my_palette <- colorRampPalette(c("yellow","red","purple","blue","green"))(n = 1000)
  library(devtools)
  source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
  heatmap.3(Wi,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none', 
            col = my_palette, ColSideColors = cbind(color.scale(X.list[[1]][1,]), color.scale(X.list[[2]][1,])))



###########################################
# Orthonormalized kernel PLS fit algorithms
#
###########################################


okernelpls.fit <- function(X, Y, ncomp, stripped=FALSE, scale=FALSE, center=TRUE, kernelfun='linkern', ...) {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    if (!stripped) {
        dnX <- dimnames(X)
        dnY <- dimnames(Y)
    }
    dimnames(X) <- dimnames(Y) <- NULL
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]
    Xmeans <- colMeans(X)
    X <- scale(X, scale=scale, center=center)
    Ymeans <- colMeans(Y)
    Y <- scale(Y, scale=scale, center=center)

    # define a linear kernel for response vectors
    K_y <- tcrossprod(Y)
    if (kernelfun == "linkern") {
        K_x <- tcrossprod(X)
    } else {
        library(kernlab)
        K_x <- kernelMatrix(match.fun(kernelfun), X)
    }
    
    Amat <- K_x %*% K_y %*% K_x
    Amat0 <- Amat
    Bmat <- K_x %*% K_x
    TT <- U <- matrix(0, ncol = ncomp, nrow = nobj)
    B <- array(0, c(npred, nresp, ncomp))
    In <- diag(nobj)
    nits <- numeric(ncomp)
    A <- matrix(0, nobj, ncomp)
    if (!stripped) {
        fitted <- array(0, dim = c(nobj, nresp, ncomp))
        Xresvar <- numeric(ncomp)
        Xtotvar <- sum(diag(K_x))
    }
    for (a in 1:ncomp) {
        # generalized eigenvalue decomposition of feature and response kernels
        out  <- .geigen(Amat, Bmat)
        vals <- out$values
        vecs <- out$vectors
        val.max <- which.max(vals)
        alpha   <- vecs[,val.max,drop=FALSE]
        alpha <- alpha/ norm(matrix(alpha), 'F')
        beta  <- K_y %*% alpha
        utmp  <- beta/c(crossprod(alpha, beta))
        wpw   <- sqrt(c(crossprod(utmp, K_x) %*% utmp))
        TT[, a] <- alpha * wpw
        U[, a]  <- utmp * wpw
        A[,a]   <- alpha
        Amat <- Amat - (vals[val.max]  * (Bmat %*% alpha %*% t(alpha) %*% Bmat))
    }
    n <- length(X)
    m <- length(X)
    ret <- t(t(K_x - rowSums(K_x)/m) - rowSums(K_x)/m) + sum(K_x)/(m * n)
    TT  <- scale(tcrossprod(ret) %*% TT, scale=TRUE, center=TRUE)[,,drop=FALSE]
    TTtTinv <- TT %*% diag(1/colSums(TT * TT), ncol = ncol(TT))
    if (kernelfun != 'linkern') {  
        XW <- tcrossprod(ret) %*% U
        W <- scale(crossprod(X, XW), center=FALSE)
        W <- W/rep(sqrt(colSums(W * W)), each = npred)
        P <- crossprod(tcrossprod(XW, W), TTtTinv)
    } else {
        W <- crossprod(X, U)
        W <- W/rep(sqrt(colSums(W * W)), each = npred)
        P <- crossprod(X, TTtTinv)
    }

    Q <- crossprod(Y, TTtTinv)
    if (ncomp == 1) {
        R <- W
    }
    else {
        PW <- crossprod(P, W)
        if (nresp == 1) {
            PWinv <- diag(ncomp)
            bidiag <- -PW[row(PW) == col(PW) - 1]
            for (a in 1:(ncomp - 1)) PWinv[a, (a + 1):ncomp] <- cumprod(bidiag[a:(ncomp -  1)])
        }
        else {
            PWinv <- backsolve(PW, diag(ncomp))
        }
        R <- W %*% PWinv
    }
    for (a in 1:ncomp) {
        B[, , a] <- tcrossprod(R[, 1:a, drop = FALSE], Q[, 1:a, drop = FALSE])
    }
    if (stripped) {
        list(coefficients = B, Xmeans = Xmeans, Ymeans = Ymeans)
    } else {
        for (a in 1:ncomp) fitted[, , a] <- tcrossprod(TT[, 1:a, drop = FALSE], Q[, 1:a, drop = FALSE])
        residuals <- -fitted + c(Y)
        fitted <- fitted + rep(Ymeans, each = nobj)
        Xvar <- diff(-c(Xtotvar, Xresvar))
        objnames <- dnX[[1]]
        if (is.null(objnames)) 
            objnames <- dnY[[1]]
        prednames <- dnX[[2]]
        respnames <- dnY[[2]]
        compnames <- paste("Comp", 1:ncomp)
        nCompnames <- paste(1:ncomp, "comps")
        dimnames(TT) <- dimnames(U) <- list(objnames, compnames)
        dimnames(R) <- dimnames(W) <- dimnames(P) <- list(prednames, 
            compnames)
        dimnames(Q) <- list(respnames, compnames)
        dimnames(B) <- list(prednames, respnames, nCompnames)
        dimnames(fitted) <- dimnames(residuals) <- list(objnames, 
            respnames, nCompnames)
        names(Xvar) <- compnames
        class(TT) <- class(U) <- "scores"
        class(P) <- class(W) <- class(Q) <- "loadings"
        library(MASS)
        modbeta <- MASS::ginv(TT) %*% Y
        list(coefficients = B, scores = TT, loadings = P, loading.weights = W, 
            Yscores = U, Yloadings = Q, projection = R, Xmeans = Xmeans, 
            Ymeans = Ymeans, fitted.values = fitted, residuals = residuals, modbeta=modbeta,
            Xvar = Xvar, Xtotvar = Xtotvar, nits = nits, A=A, Kx=K_x, Ky=K_y, ret=ret)
    }


#    if (kernelfun != 'linkern') {  
#        n <- length(X)
#        m <- length(X)
#        ret <- t(t(K_x - rowSums(K_x)/m) - rowSums(K_x)/m) + sum(K_x)/(m * n)
#        XU <- ret %*% t(ret) %*% A
#        # mimic feature loadings using kernel factor analysis on the scores
#        # and computing inner prod with X and scores
##        fa  <- kfa(XU, features=npred)
##        U <- fa@xmatrix
#        U <- scale(t(X) %*% XU, center=FALSE)
#    } else {
#        ret   <- X
#        U     <- t(ret) %*% A
#        XU <- ret %*% U
#    }

#        rownames(U) <- colnames(X)
#        XUXU_inv <- tryCatch(solve(t(XU) %*% XU), error=function(e) {
#           require(MASS) ; ginv(t(XU) %*% XU)    })

#        B <-  XUXU_inv %*% t(XU) %*% Y
#        residuals <- Y - (XU %*% B)


#    list(coefficients = B, scores = XU, loadings = U, Xmeans = Xmeans, fitted.values=XU%*%B,
#         Ymeans = Ymeans, residuals = residuals,  Bmat=Bmat, Amat=Amat0,
#         Xvar = colSums(U * U), Xtotvar = sum(X * X), ret=ret)
}



.geigen <- function(A, B, rank=ncol(A)) {
    svd1 <- svd(B)
    F <- svd1$v
    D <- diag(svd1$d[1:rank]^-.5)
    T1 <- F[,1:rank] %*% D
    svd2 <- svd(t(T1) %*% A %*% T1)
    T2 <- svd2$v
    T3 <- T1 %*% T2
    list(vectors=T3, values=diag(t(T3) %*% A %*% T3))
}




#kpls <- function (X, Y, ncomp, stripped = FALSE, tol = .Machine$double.eps^0.5, maxit = 100, ...) {
#    Y <- as.matrix(Y)
#    if (!stripped) {
#        dnX <- dimnames(X)
#        dnY <- dimnames(Y)
#    }
#    dimnames(X) <- dimnames(Y) <- NULL
#    nobj <- dim(X)[1]
#    npred <- dim(X)[2]
#    nresp <- dim(Y)[2]
#    TT <- U <- matrix(0, ncol = ncomp, nrow = nobj)
#    B <- array(0, c(npred, nresp, ncomp))
#    In <- diag(nobj)
#    nits <- numeric(ncomp)
#    if (!stripped) {
#        fitted <- array(0, dim = c(nobj, nresp, ncomp))
#        Xresvar <- numeric(ncomp)
#    }
#    Xmeans <- colMeans(X)
#    X <- X - rep(Xmeans, each = nobj)
#    Ymeans <- colMeans(Y)
#    Y <- Y - rep(Ymeans, each = nobj)
#    XXt <- tcrossprod(X)
#    YYt <- tcrossprod(Y)
#    if (!stripped) 
#        Xtotvar <- sum(diag(XXt))
#    for (a in 1:ncomp) {
#        XXtYYt <- XXt %*% YYt
#        XXtYYt <- XXtYYt %*% XXtYYt
#        t.a.old <- Y[, 1]
#        nit <- 0
#        repeat {
#            nit <- nit + 1
#            t.a <- XXtYYt %*% t.a.old
#            t.a <- t.a/sqrt(c(crossprod(t.a)))
#            if (sum(abs((t.a - t.a.old)/t.a), na.rm = TRUE) < 
#                tol) 
#                break
#            else t.a.old <- t.a
#            if (nit >= maxit) {
#                warning("No convergence in", maxit, "iterations\n")
#                break
#            }
#        }
#        nits[a] <- nit
#        u.a <- YYt %*% t.a
#        utmp <- u.a/c(crossprod(t.a, u.a))
#        wpw <- sqrt(c(crossprod(utmp, XXt) %*% utmp))
#        TT[, a] <- t.a * wpw
#        U[, a] <- utmp * wpw
#        G <- In - tcrossprod(t.a)
#        XXt <- G %*% XXt %*% G
#        YYt <- G %*% YYt %*% G
#        if (!stripped) 
#            Xresvar[a] <- sum(diag(XXt))
#    }
#    W <- crossprod(X, U)
#    W <- W/rep(sqrt(colSums(W * W)), each = npred)
#    TTtTinv <- TT %*% diag(1/colSums(TT * TT), ncol = ncol(TT))
#    P <- crossprod(X, TTtTinv)
#    Q <- crossprod(Y, TTtTinv)
#    if (ncomp == 1) {
#        R <- W
#    }
#    else {
#        PW <- crossprod(P, W)
#        if (nresp == 1) {
#            PWinv <- diag(ncomp)
#            bidiag <- -PW[row(PW) == col(PW) - 1]
#            for (a in 1:(ncomp - 1)) PWinv[a, (a + 1):ncomp] <- cumprod(bidiag[a:(ncomp - 
#                1)])
#        }
#        else {
#            PWinv <- backsolve(PW, diag(ncomp))
#        }
#        R <- W %*% PWinv
#    }
#    for (a in 1:ncomp) {
#        B[, , a] <- tcrossprod(R[, 1:a, drop = FALSE], Q[, 1:a, 
#            drop = FALSE])
#    }
#    if (stripped) {
#        list(coefficients = B, Xmeans = Xmeans, Ymeans = Ymeans)
#    }
#    else {
#        for (a in 1:ncomp) fitted[, , a] <- tcrossprod(TT[, 1:a, 
#            drop = FALSE], Q[, 1:a, drop = FALSE])
#        residuals <- -fitted + c(Y)
#        fitted <- fitted + rep(Ymeans, each = nobj)
#        Xvar <- diff(-c(Xtotvar, Xresvar))
#        objnames <- dnX[[1]]
#        if (is.null(objnames)) 
#            objnames <- dnY[[1]]
#        prednames <- dnX[[2]]
#        respnames <- dnY[[2]]
#        compnames <- paste("Comp", 1:ncomp)
#        nCompnames <- paste(1:ncomp, "comps")
#        dimnames(TT) <- dimnames(U) <- list(objnames, compnames)
#        dimnames(R) <- dimnames(W) <- dimnames(P) <- list(prednames, 
#            compnames)
#        dimnames(Q) <- list(respnames, compnames)
#        dimnames(B) <- list(prednames, respnames, nCompnames)
#        dimnames(fitted) <- dimnames(residuals) <- list(objnames, 
#            respnames, nCompnames)
#        names(Xvar) <- compnames
#        class(TT) <- class(U) <- "scores"
#        class(P) <- class(W) <- class(Q) <- "loadings"
#        list(coefficients = B, scores = TT, loadings = P, loading.weights = W, 
#            Yscores = U, Yloadings = Q, projection = R, Xmeans = Xmeans, 
#            Ymeans = Ymeans, fitted.values = fitted, residuals = residuals, 
#            Xvar = Xvar, Xtotvar = Xtotvar, nits = nits)
#    }
#}




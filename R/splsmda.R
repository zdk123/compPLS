#######################
# multilabel sPLS-DA
#
#########################


splsmda <- function (x, y, probMethod = "softmax", prior = NULL, ...) {
    caret:::requireNamespaceQuietStop("spls")
    funcCall <- match.call(expand.dots = TRUE)
    if (!is.matrix(x)) 
        x <- as.matrix(x)
    if (probMethod == "softmax") {
        if (!is.null(prior)) 
            warning("Priors are ignored unless probMethod = \"Bayes\"")
    }
    if (is.factor(y)) {
        obsLevels <- list(levels(y))
        oldY <- list(y)
        y <- .class2ind(y)
    } else if (is.data.frame(y)) {
      ## coerce to columns to factors
        for (i in 1:ncol(y)) y[,i] <- droplevels(as.factor(y[,i]))
        obsLevels <- lapply(y, levels)
        oldY <- y
        y <- do.call('cbind', lapply(y, .class2ind))
    }
    else {
        if (is.matrix(y)) {
#            test <- apply(y, 1, sum)
#            if (any(test != 1)) 
#                stop("the rows of y must be 0/1 and sum to 1")
            obsLevels <- list(colnames(y))
            if (is.null(obsLevels)) 
                stop("the y matrix must have column names")
            oldY <- obsLevels[apply(y, 1, which.max)]
        }
        else stop("y must be a matrix, data.frame or a factor")
    }
    tmpData <- data.frame(n = paste("row", 1:nrow(y), sep = ""))
    tmpData$y <- y
    tmpData$x <- x
    out <- .spls(x, y, ...)
    out$obsLevels <- obsLevels
    out$probMethod <- probMethod
    if (probMethod == "Bayes") {
        caret:::requireNamespaceQuietStop("klaR")
        makeModels <- function(x, y, pri) {
            probModel <- klaR::NaiveBayes(x, y, prior = pri, usekernel = TRUE)
            probModel$train <- predict(probModel)$posterior
            probModel$x <- NULL
            probModel
        }
        cls <- class(out)
        ## spoof # obslevels so nothing is dropped in the predict function
        out$obsLevels <- ncol(y) + 1
        train <- spls:::predict.spls(out, as.matrix(tmpData$x))
        usedlevs <- 0
        probMod <- vector("list", N <- length(obsLevels))
        for (i in 1:N) {
            ysub  <- oldY[[i]]
            nlevs <- nlevels(ysub)
            traintemp <- train[, (i+usedlevs):(nlevs+usedlevs)]
            traintemp <- traintemp[, -length(obsLevels[[i]]), drop = FALSE]
            usedlevs  <- usedlevs + nlevels(ysub)
            probMod[[i]] <- makeModels(traintemp, y = ysub, pri = prior)
        }
        out$probModel <- probMod
        names(out$probModel) <- names(oldY)
    }
    else out$probModel <- NULL
    class(out) <- c("splsmda", "splsda")
    out
}


.spls <- function (x, y, K, eta, kappa = 0.5, select = "pls2", fit = "okernelpls", 
    scale.x = TRUE, scale.y = FALSE, eps = 1e-04, maxstep = 100, trace = FALSE, ...) {
    caret:::requireNamespaceQuietStop("spls")
    caret:::requireNamespaceQuietStop("pls")
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ip <- c(1:p)
    y <- as.matrix(y)
    q <- ncol(y)
    one <- matrix(1, 1, n)
    mu <- one %*% y/n
    y <- scale(y, drop(mu), FALSE)
    meanx <- drop(one %*% x)/n
    x <- scale(x, meanx, FALSE)
    if (scale.x) {
        normx <- sqrt(drop(one %*% (x^2))/(n - 1))
        if (any(normx < .Machine$double.eps)) {
            stop("Some of the columns of the predictor matrix have zero variance.")
        }
        x <- scale(x, FALSE, normx)
    }
    else {
        normx <- rep(1, p)
    }
    if (scale.y) {
        normy <- sqrt(drop(one %*% (y^2))/(n - 1))
        if (any(normy < .Machine$double.eps)) {
            stop("Some of the columns of the response matrix have zero variance.")
        }
        y <- scale(y, FALSE, normy)
    }
    else {
        normy <- rep(1, q)
    }
    betahat <- matrix(0, p, q)
    betamat <- list()
    x1 <- x
    y1 <- y
    type <- spls:::correctp(x, y, eta, K, kappa, "simpls", "simpls")
    eta <- type$eta
    K <- type$K
    kappa <- type$kappa
 #   select <- type$select
    if (is.null(colnames(x))) {
        xnames <- c(1:p)
    }
    else {
        xnames <- colnames(x)
    }
    new2As <- list()
    if (trace) {
        cat("The variables that join the set of selected variables at each step:\n")
    }
    for (k in 1:K) {
        Z <- t(x1) %*% y1
        what <- spls:::spls.dv(Z, eta, kappa, eps, maxstep)
        A <- unique(ip[what != 0 | betahat[, 1] != 0])
        new2A <- ip[what != 0 & betahat[, 1] == 0]
        xA <- x[, A, drop = FALSE]
        plsfit <- .plsr(y ~ xA, ncomp = min(k, length(A)), 
            method = fit, scale = FALSE, ...)
        betahat <- matrix(0, p, q)
        betahat[A, ] <- matrix(coef(plsfit), length(A), q)
        betamat[[k]] <- betahat
        pj <- plsfit$projection
        scores <- plsfit$scores
        if (select == "pls2") {
            y1 <- y - x %*% betahat
        }
        if (select == "simpls") {
            pw <- pj %*% solve(t(pj) %*% pj) %*% t(pj)
            x1 <- x
            x1[, A] <- x[, A, drop = FALSE] - x[, A, drop = FALSE] %*% pw
        }
        if (select == "okernelpls") {
            modbeta <- plsfit$modbeta
            y1 <- y - scores %*% modbeta
        }
        new2As[[k]] <- new2A
        if (trace) {
            if (length(new2A) <= 10) {
                cat(paste("- ", k, "th step (K=", k, "):\n", 
                  sep = ""))
                cat(xnames[new2A])
                cat("\n")
            }
            else {
                cat(paste("- ", k, "th step (K=", k, "):\n", 
                  sep = ""))
                nlines <- ceiling(length(new2A)/10)
                for (i in 0:(nlines - 2)) {
                  cat(xnames[new2A[(10 * i + 1):(10 * (i + 1))]])
                  cat("\n")
                }
                cat(xnames[new2A[(10 * (nlines - 1) + 1):length(new2A)]])
                cat("\n")
            }
        }
    }
    if (!is.null(colnames(x))) {
        rownames(betahat) <- colnames(x)
    }
    if (q > 1 & !is.null(colnames(y))) {
        colnames(betahat) <- colnames(y)
    }
    object <- list(x = x, y = y, betahat = betahat, A = A, betamat = betamat, 
        new2As = new2As, mu = mu, meanx = meanx, normx = normx,
        normy = normy, eta = eta, K = K, kappa = kappa, select = select, 
        fit = fit, projection = pj,  scores=plsfit$scores)
    class(object) <- "spls"
    object
}


predict.splsmda <- function (object, newdata, type = "class", ...) {
    caret:::requireNamespaceQuietStop("spls")
    tmpPred <- spls::predict.spls(object, newx = newdata)
    if (type == "raw") 
        return(tmpPred)
    if (is.null(object$probModel)) {
        # softmax
        usedlevs <- 1
        out <- vector("list", N <- length(object$obsLevels))
        for (i in 1:N) {
            levs <- object$obsLevels[[i]]
            nlevs     <- length(levs)
            tmpPred2  <- as.data.frame(tmpPred[, (usedlevs):(nlevs+usedlevs-1)])
            out[[i]]  <- switch(type, class = {
                classIndex <- levs[apply(tmpPred2, 1, which.max)]
                factor(classIndex, levels = levs)
             }, prob = t(apply(tmpPred2, 1, function(data) exp(data)/sum(exp(data))))
            )
            usedlevs  <- usedlevs + nlevs
        }
    }
    else {
        # Naive Bayes
        caret:::requireNamespaceQuietStop("klaR")
        usedlevs <- 1
        out <- vector("list", N <- length(object$obsLevels))
        for (i in 1:N) {
            levs <- object$obsLevels[[i]]
            nlevs     <- length(levs)
            tmpPred2  <- as.data.frame(tmpPred[, (usedlevs):(nlevs+usedlevs-1)])
            tmpPred2  <- tmpPred2[, -nlevs, drop = FALSE]
            pred      <- predict(object$probModel[[i]], tmpPred2)
            out[[i]]  <- switch(type, class = pred$class, prob = pred$posterior)
            usedlevs  <- usedlevs + nlevs
        }
    }
    names(out) <- names(object$obsLevels)
    out
}






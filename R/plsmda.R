#######################
# multilabel PLS-DA
#
#########################


plsmda <- function (x, y, ncomp = 2, probMethod = "softmax", prior = NULL, ...) {
## partial least squares multilabel discriminant analysis
    caret:::requireNamespaceQuietStop("pls")
    funcCall <- match.call(expand.dots = TRUE)
    if (!is.matrix(x)) 
        x <- as.matrix(x)
    if (length(ncomp) > 1) {
        ncomp <- max(ncomp)
        warning(paste("A value single ncomp must be specified.", 
            "max(ncomp) was used.", "Predictions can be obtained for values <= ncomp"))
    }
    if (probMethod == "softmax") {
        if (!is.null(prior)) 
            warning("Priors are ignored unless probMethod = \"Bayes\"")
    }
    if (is.factor(y)) {
        obsLevels <- levels(y)
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
            obsLevels <- colnames(y)
            if (is.null(obsLevels)) 
                stop("the y matrix must have column names")
            oldY <- obsLevels[apply(y, 1, which.max)]
        }
        else stop("y must be a matrix, data.frame or a factor")
    }
    tmpData <- data.frame(n = paste("row", 1:nrow(y), sep = ""))
    tmpData$y <- y
    tmpData$x <- x
    out <- .plsr(y ~ x, data = tmpData, ncomp = ncomp, ...)
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
        class(out) <- "mvr"
        train <- predict(out, as.matrix(tmpData$x), ncomp = 1:ncomp)
#       train <- train[, -length(obsLevels), , drop = FALSE]
        usedlevs <- 0
        probMod <- vector("list", N <- length(obsLevels))
        for (i in 1:N) {
            ysub  <- oldY[[i]]
            nlevs <- nlevels(ysub)
            traintemp <- train[, (i+usedlevs):(nlevs+usedlevs), ,drop=FALSE]
            traintemp <- traintemp[, -length(obsLevels[[i]]), , drop = FALSE]
            usedlevs  <- usedlevs + nlevels(ysub)
            probMod[[i]] <- apply(traintemp, 3, makeModels, y = ysub, pri = prior)
        }
        out$probModel <- probMod
        names(out$probModel) <- names(oldY)
    }
    else out$probModel <- NULL
    class(out) <- c("plsmda", "plsda", class(out))  
    out
}


.plsr <- function (..., method = "okernelpls") {
    cl <- match.call()
    cl$method <- match.arg(method, c("kernelpls", "widekernelpls", "okernelpls",
           "simpls", "oscorespls", "model.frame"))
    cl[[1]] <- quote(.mvr)
    res <- eval(cl, parent.frame())
    if (cl$method != "model.frame") 
        res$call[[1]] <- as.name(".plsr")
    if (missing(method)) 
        res$call$method <- NULL
    res
}


.mvr <- function (formula, ncomp, Y.add, data, subset, na.action, method = pls.options()$mvralg, 
    scale = FALSE, validation = c("none", "CV", "LOO"), model = TRUE, 
    x = FALSE, y = FALSE, ...) {
    caret:::requireNamespaceQuietStop("pls")
    ret.x <- x
    ret.y <- y
    mf <- match.call(expand.dots = FALSE)
    if (!missing(Y.add)) {
        Y.addname <- as.character(substitute(Y.add))
        mf$formula <- update(formula, paste("~ . +", Y.addname))
    }
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 
        0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    method <- match.arg(method, c("kernelpls", "widekernelpls", "okernelpls",
        "simpls", "oscorespls", "cppls", "svdpc", "model.frame"))
    if (method == "model.frame") 
        return(mf)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "numeric")
    if (is.matrix(Y)) {
        if (is.null(colnames(Y))) 
            colnames(Y) <- paste("Y", 1:dim(Y)[2], sep = "")
    }

    else {
        Y <- as.matrix(Y)
        colnames(Y) <- deparse(formula[[2]])
    }
    if (missing(Y.add)) {
        Y.add <- NULL
    }
    else {
        Y.add <- mf[, Y.addname]
        mt <- drop.terms(mt, which(attr(mt, "term.labels") == 
            Y.addname), keep.response = TRUE)
    }
    X <- pls:::delete.intercept(model.matrix(mt, mf))
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    if (length(attr(mt, "term.labels")) == 1 && !is.null(colnames(mf[[attr(mt, 
        "term.labels")]]))) 
        colnames(X) <- sub(attr(mt, "term.labels"), "", colnames(X))
    if (missing(ncomp)) {
        ncomp <- min(nobj - 1, npred)
        ncompWarn <- FALSE
    }
    else {
        if (ncomp < 1 || ncomp > min(nobj - 1, npred)) 
            stop("Invalid number of components, ncomp")
        ncompWarn <- TRUE
    }
    sdscale <- identical(TRUE, scale)
    if (is.numeric(scale)) 
        if (length(scale) == npred) 
            X <- X/rep(scale, each = nobj)
        else stop("length of 'scale' must equal the number of x variables")
    switch(match.arg(validation), CV = {
        val <- mvrCv(X, Y, ncomp, Y.add = Y.add, method = method, 
            scale = sdscale, ...)
    }, LOO = {
        segments <- as.list(1:nobj)
        attr(segments, "type") <- "leave-one-out"
        val <- mvrCv(X, Y, ncomp, Y.add = Y.add, method = method, 
            scale = sdscale, segments = segments, ...)
    }, none = {
        val <- NULL
    })
    if (identical(TRUE, ncomp > val$ncomp)) {
        ncomp <- val$ncomp
        if (ncompWarn) 
            warning("`ncomp' reduced to ", ncomp, " due to cross-validation")
    }
    fitFunc <- switch(method, kernelpls = pls:::kernelpls.fit, widekernelpls = pls:::widekernelpls.fit, 
        simpls = pls:::simpls.fit, oscorespls = pls:::oscorespls.fit, 
        cppls = pls:::cppls.fit, okernelpls=okernelpls.fit, svdpc = pls:::svdpc.fit)
    if (sdscale) {
        scale <- sqrt(colSums((X - rep(colMeans(X), each = nobj))^2)/(nobj - 
            1))
        if (any(abs(scale) < .Machine$double.eps^0.5)) 
            warning("Scaling with (near) zero standard deviation")
        X <- X/rep(scale, each = nobj)
    }
    start.time <- proc.time()[3]
    z <- fitFunc(X, Y, ncomp, Y.add = Y.add, ...)
    z$fit.time <- proc.time()[3] - start.time
    class(z) <- "mvr"
    z$na.action <- attr(mf, "na.action")
    z$ncomp <- ncomp
    z$method <- method
    if (is.numeric(scale)) 
        z$scale <- scale
    z$validation <- val
    z$call <- match.call()
    z$terms <- mt
    if (model) 
        z$model <- mf
    if (ret.x) 
        z$x <- X
    if (ret.y) 
        z$y <- Y
    z
}

print.plsmda <- function(x, ...) {
    if (x$method == "okernelpls") {
        ana = "Orthonormalized Partial Least Squares"
        alg = "kernel"
        kernfun = x$call$kernelfun
       cat(ana, ", fitted with the", alg, "algorithm and using a", kernfun, "kernel function")
    if (!is.null(x$validation)) 
        cat("\nCross-validated using", length(x$validation$segments), 
            attr(x$validation$segments, "type"), "segments.")
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    invisible(x)
    } else {
        pls:::print.mvr(x)
    }
}


predict.plsmda <- function (object, newdata = NULL, ncomp = NULL, type = "class",  ...)  {
    caret:::requireNamespaceQuietStop("pls")
    if (is.null(ncomp)) 
        if (!is.null(object$ncomp)) 
            ncomp <- object$ncomp
        else stop("specify ncomp")
    if (!is.null(newdata)) {
        if (!is.matrix(newdata)) 
            newdata <- as.matrix(newdata)
    }
    class(object) <- "mvr"
    tmpPred <- predict(object, newdata = newdata)[, , ncomp, drop = FALSE]
    if (type == "raw") 
        return(tmpPred)
    out <- vector('list', N <- length(object$obsLevels))
    if (is.null(object$probModel)) {
        ## softmax
        for (j in 1:N) {
            switch(type, class = {
                if (length(dim(tmpPred)) < 3) {
                    # only 1 latent component 
                    out[[j]] <- object$obsLevels[apply(tmpPred, 1, which.max)]
                    out[[j]] <- factor(out[[j]], levels = object$obsLevels[[j]])
                } else {
                    tmpOut <- matrix("", nrow = dim(tmpPred)[1], ncol = dim(tmpPred)[3])
                    for (i in 1:dim(tmpPred)[3]) {
                      tmpOut[, i] <- object$obsLevels[[j]][apply(tmpPred[, , i, drop = FALSE], 1, which.max)]
                    }
                    out[[j]] <- as.data.frame(tmpOut)
                    out[[j]] <- as.data.frame(lapply(out[[j]], function(x, y) factor(x, levels = y), y = object$obsLevels[[j]]))
                    names(out[[j]]) <- paste("ncomp", ncomp, sep = "")
                    rownames(out[[j]]) <- rownames(newdata)
                    if (length(ncomp) == 1) out[[j]] <- out[[j]][, 1]
                }
            }, prob = {
                if (length(dim(tmpPred)) < 3) {
                    out[[j]] <- t(apply(tmpPred, 1, function(data) exp(data)/sum(exp(data))))
                } else {
                    out[[j]] <- tmpPred * NA
                    for (i in 1:dim(tmpPred)[3]) {
                      out[[j]][, , i] <- t(apply(tmpPred[, , i, drop = FALSE], 1, function(data) exp(data)/sum(exp(data))))
                    }
                }
            })
        }
    }
    else {
        for (j in 1:N) {
            library(klaR)
            tmp <- vector(mode = "list", length = length(ncomp))
            for (i in seq(along = ncomp)) {
                tmp[[i]] <- predict(object$probModel[[j]][[ncomp[i]]], as.data.frame(tmpPred[, -length(object$obsLevels[[j]]), i]))
            }
            if (type == "class") {
                out[[j]] <- t(do.call("rbind", lapply(tmp, function(x) as.character(x$class))))
                rownames(out[[j]]) <- names(tmp[[1]]$class)
                colnames(out[[j]]) <- paste("ncomp", ncomp, sep = "")
                out[[j]] <- as.data.frame(out[[j]])
                out[[j]] <- as.data.frame(lapply(out[[j]], function(x, y) factor(x, levels = y), y = object$obsLevels[[j]]))
                if (length(ncomp) == 1) 
                    out[[j]] <- out[[j]][, 1]
            } else {
                out[[j]] <- array(dim = c(dim(tmp[[1]]$posterior), length(ncomp)), 
                    dimnames = list(rownames(tmp[[1]]$posterior), 
                       colnames(tmp[[1]]$posterior), paste("ncomp", ncomp, sep = "")))
                for (i in seq(along = ncomp)) out[[j]][, , i] <- tmp[[i]]$posterior
            }
        }
    }
    names(out) <- names(object$obsLevels)
    out
}


.class2ind <- function(y) {
    y <- caret:::class2ind(y)
    y[y==0] <- -1
    y
}

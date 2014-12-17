################################
# @author Zachary Kurtz
# @date 8/14/2013
## Methods for Partial Least Squares regression for compositional data




#' Partial Least Squares Discriminant Analysis
#' PLS regression to discriminate classes (via a logistic model)
#' basically this is a wrapper for the \code{plsda} function in the caret package,
#' but with default setup for dealing with uneven classes (via the priors option, see details)
#' see caret::plsda for implementation details
#' 
#' run this code if you don't need to fit paramaters by cross-validation
#'
#' @title plsDA partial least squares discriminant analysis
#' @param x data with samples in rows, features are columns (not necessarily compositional data)
#' @param grouping a numeric vector or factor with sample classes (length should equal \code{nrow(x)})
#' @param usePriors use priors for very biased sample size between groups (ie - put strong penalty on misclassifying small groups)
#' @param K number of components in the PLS model (default: number of classes - 1)
#' @return a plsda fitted model
#' @seealso \code{\link{plsDA_main}}, \code{\link{caret::plsda}}
#' @rdname plsDA
#' @export
plsDA <- function(data, grouping, K, usePriors=FALSE, ...) {
# wrapper function for plsda

    if (length(grouping) != nrow(data)) 
        stop('length of grouping vector does equal the number of samples')

    args <- list(...)
    if (!('probMethod' %in% names(args))) args <- c(args, list(probMethod='Bayes'))
    if (usePriors) {
        # set priors to 1-(class freq) 
        prior <- 1-(table(grouping)/length(grouping))
        args <- c(args, list(prior=as.vector(prior)))
    }

    y <- as.factor(grouping)

    if (missing(K)) {
        K <- length(levels(y))-1
    }
    args  <- c(args, ncomp=K)
    
    do.call(caret::plsda, c(list(data, y), args))
}


#' The main wrapper for full Partial Least Squares discriminant analysis,
#' performing cross-validation to tune model parameters (here, number of components)
#' and do permutation tests (ie bootstrapping) to get pseudo-pvals estimates for model coefficients
#'
#' @title plsDA_main partial least squares discriminant analysis
#' @param x data with samples in rows, features are columns (not necessarily compositional data)
#' @param grouping a numeric vector or factor with sample classes (length should equal \code{nrow(x)})
#' @param usePriors use priors for very biased sample size between groups (ie - put strong penalty on misclassifying small groups)
#' @param K numeric vector containing number of components in the PLS model
#' @param  fold number of partitions to randomly subsample for cross-validation
#' @param nboots number of bootstraps/permutations for estimating coefficient p-vals
#' @param n.core number of cores for paralellization of bootstraps
#' @param noise for very sparse components, some subsamples may have zero variance. Optionally, add some Gaussian noise to to avoid PLS errors
#' @param ... additional arguments passed to plsDA
#' 
#' @return a \code{plsDA} object that contains: the plsda model/object, \code{pvals}, the original data, \code{x}, and \code{groupings}
#' @seealso \code{\link{plsDA}}
#' @rdname plsDA_main
#' @export
plsDA_main  <- function(x, grouping, K, usePriors=FALSE, fold=5, nboots=999, n.core=4, noise=0, ...) {

    K   <- cv.plsDA.fit(x, grouping, eta=eta, K=K, noise=noise, fold=fold, usePriors, n.core=n.core)$K
    .bstat <- function(data, indices, ...)   plsDA(data[indices,,drop=FALSE], ...)$coefficients
    .pstat <- function(data, indices, ...)   plsDA(apply(data[indices,,drop=FALSE], 2, function(x) sample(x)), ...)$coefficients

    if (nboots > 1) {
        sboots <- plsdaboot(x, .bstat, .pstat, R=nboots, n.core, K=K, grouping=grouping, ...)
        pmat   <- pval(sboots)
#        keepInd <- which(suppressWarnings(apply(pmat, 1, min, na.rm=TRUE)) <= alpha)
    } else {
        keepInd <- 1:ncol(x)
        pmat    <- NULL
    }
    plsmod  <- plsDA(x, grouping, usePriors=usePriors, K=K, ...)

#    if (length(keepInd) > 0) {
#        splskmod <- plsDA(x[,keepInd], grouping, usePriors=usePriors, K=K, ...)	

#    } else splskmod <- NULL
    structure(list(plsda=plsmod, pvals=pmat, # keep=keepInd,
                  x=data, y=grouping), class="plsDA")


}


plsdaboot <- function(data, statisticboot, statisticperm, R, ncpus=1, ...) {
    res     <- boot::boot(data, statisticboot, R=R, parallel="multicore", ncpus=ncpus, noise=0, ...)
    null_av <- boot::boot(data, statisticperm, sim='permutation', R=R, parallel="multicore", ncpus=ncpus, ...)
    class(res) <- 'list'
    structure(c(res, list(null_av=null_av)), class='plsdaboot')
}



pval.plsdaboot <- function(x, sided='both', mar=2) {
# calculate 1 or 2 way pseudo p-val from boot object
# Args: a boot object
    if (sided != "both") stop("only two-sided currently supported")
    x$t0 <- as.matrix(x$t0)
    nparams  <- ncol(x$t)
    tmeans   <- colMeans(x$null_av$t)
#    check to see whether betas are unstable -- confirm 
    niters   <- nrow(x$t)
    ind95    <- max(1,round(.025*niters)):round(.975*niters)
    boot_ord <- apply(x$t, 2, sort)
    boot_ord95 <- boot_ord[ind95,]
    outofrange <- unlist(lapply(1:length(x$t0), function(i) {
            betas <- x$t0[i]
            range  <- range(boot_ord95[,i])
            range[1] > betas || range[2] < betas
        }))
    # calc whether center of mass is above or below the mean
    bs_above <- unlist(lapply(1:nparams, function(i) 
                    length(which(x$t[, i] > tmeans[i]))))
    is_above <- bs_above > x$R/2
    pvals    <- ifelse(is_above, 2*(1-bs_above/x$R), 2*bs_above/x$R)
    pvals[pvals > 1]  <- 1
    pvals[outofrange] <- NaN
    pmat <- matrix(pvals, ncol=ncol(x$t0))
    rownames(pmat) <- rownames(x$t0)
    colnames(pmat) <- colnames(x$t0)
    pmat
}

cv.plsDA.fit <- function(x, grouping, eta, K, noise=0, fold=5, usePriors=FALSE, n.core) {
# Find eta & k by cross validation
    x <- x + matrix(rnorm(prod(dim(x)), 0, noise), ncol=ncol(x))
    capture.output(out <- .cv.plsDA(x, grouping, fold=fold, K=K, plot.it=FALSE, n.core=n.core))
    out
}

#' @keywords internal
.cv.plsDA <- function (x, y, fold = 10, K, kappa = 0.5, classifier = c("lda", 
    "logistic"), scale.x = TRUE, n.core = 4, plot.it=FALSE, ...) 
{
    result.mat <- c()
    foldi <- .cv.split(y, fold)
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ip <- c(1:p)
    y <- as.matrix(y)
    q <- ncol(y)
#    eta.K.pair <- cbind(rep(1, each = length(K)), rep(K, 1))
#    eta.K.list <- split(eta.K.pair, c(1:nrow(eta.K.pair)))
    .fit.plsda <- function(K.val) {
        mspemati <- rep(0, fold)
        Ai <- rep(0, fold)
        for (k in 1:fold) {
            omit <- foldi[[k]]
            train.x <- x[-omit, ]
            train.y <- y[-omit, ]
            test.x <- x[omit, ]
            test.y <- y[omit, ]
            plsda.fit <- plsDA(train.x, train.y, K = K.val, 
                            scale.x = scale.x, classifier = classifier, ...)
            pred <- as.numeric(predict(plsda.fit, 
                newx = test.x))
            mspemati[k] <- mean(as.numeric(pred != test.y))
            Ai[k] <- mean(length(plsda.fit$A))
        }
        mspe.ij <- c(mean(mspemati), mean(Ai), 1, K.val)
        print(mspe.ij)
        return(mspe.ij)
    }
    if (.Platform$OS.type == "unix") {
        result.list <- parallel::mclapply(K, function(k) .fit.plsda(k), mc.cores = n.core)
    }
    else {
        result.list <- lapply(K, function(k) .fit.plsda(k))
    }
    result.mat <- c()
    for (i in 1:length(result.list)) {
        result.mat <- rbind(result.mat, result.list[[i]])
    }
    
    mspemat <- matrix(result.mat[, 1], length(K), 1)
    mspemat <- t(mspemat)
    rownames(mspemat) <- '1'
    colnames(mspemat) <- K
    cands <- result.mat[result.mat[, 1] == min(result.mat[, 1]), 
        , drop = FALSE]

    cands <- cands[cands[, 2] == min(cands[, 2]), , drop = FALSE]
    cands <- cands[cands[, 4] == min(cands[, 4]), , drop = FALSE]
    cands <- cands[cands[, 3] == max(cands[, 3]), , drop = FALSE]
    K.opt <- cands[, 4]
    eta.opt <- cands[, 3]

#    cat(paste("K = ", K.opt, "\n", sep = ""))
#    if (plot.it) {

#        spls::heatmap.spls(t(mspemat), xlab = "K", ylab = "eta", main = "CV MSPE Plot", 
#            coln = 16, as = "n")
#    }
#    rownames(mspemat) <- paste("eta=", eta)
    colnames(mspemat) <- paste("K =", K)
    cv <- list(mspemat = mspemat,  K.opt = K.opt)
    invisible(cv)
}

#' @keywords internal
.cv.split <- function (y, fold) {
    n <- length(y)
    group <- table(y)
    x <- c()
    for (i in 1:length(group)) {
        x.group <- c(1:n)[y == names(group)[i]]
        x <- c(x, sample(x.group))
    }
    foldi <- split(x, rep(1:fold, length = n))
    return(foldi)
}

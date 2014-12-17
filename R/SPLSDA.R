########################################
# SPLS discriminant analysis pipeline
#
#




#' sparse Partial Least Squares Discriminant Analysis
#' sPLS regression to discriminate classes (via a logistic model)
#' basically this is a wrapper for the \code{splsda} function in the caret package,
#' but with default setup for dealing with uneven classes (via the priors option, see details)
#' see caret::splsda for implementation details
#' 
#' run this code if you don't need to fit paramaters by cross-validation
#'
#' @title splsDA sparse partial least squares discriminant analysis
#' @param x x with samples in rows, features are columns (not necessarily compositional x)
#' @param grouping a numeric vector or factor with sample classes (length should equal \code{nrow(x)})
#' @param usePriors use priors for very biased sample size between groups (ie - put strong penalty on misclassifying small groups)
#' @param K number of components in the PLS model (default: number of classes - 1)
#' @param eta parameter that adjusts sparsity of the PLS model (between 0 and 1)
#' @return a plsda fitted model
#' @seealso \code{\link{plsDA_main}}, \code{\link{caret::plsda}}, \code{\link{caret::splsda}}
#' @rdname splsDA
#' @export
splsDA <- function(x, grouping, eta, K, usePriors=FALSE, ...) {
## splsda wrapper
 #   caret::splsda(...)

    if (length(grouping) != nrow(x)) 
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
        ncomp <- length(levels(y))-1
    }
    args  <- c(args, K=K, eta=eta)
    do.call(caret::splsda, c(list(x, y), args))
}



#' The main wrapper for full sparse Partial Least Squares discriminant analysis,
#' performing cross-validation to tune model parameters (here, number of components)
#' and do permutation tests (ie bootstrapping) to get pseudo-pvals estimates for model coefficients
#'
#' @title splsDA_main sparse partial least squares discriminant analysis
#' @param x data with samples in rows, features are columns (not necessarily compositional x)
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
splsDA_main <- function(x, grouping, eta, K, usePriors=FALSE, fold=5, nboots=999, n.core=4, noise=0, ...) {
    # find optimal eta & K by cross validation
    
    opt   <- cv.splsDA.fit(x, grouping, eta=eta, K=K, noise=noise, fold=fold, usePriors=usePriors, n.core=n.core)
    eta   <- opt$eta.opt
    K     <- opt$K.opt
    .bstat <- function(x, indices, ...)   splsDA(x[indices,,drop=FALSE], ...)$betahat
    .pstat <- function(x, indices, ...)   splsDA(apply(x[indices,], 2, function(x) sample(x)), ...)$betahat
    
    if (nboots > 1) {
        sboots <- splsdaboot(x, .bstat, .pstat, R=nboots, n.core, eta=eta, K=K, grouping=grouping, ...)
        pmat   <- pval(sboots)
#        keepInd <- which(suppressWarnings(apply(pmat, 1, min, na.rm=TRUE)) <= alpha)
    } else {
        keepInd <- 1:ncol(x)
        pmat    <- NULL
    }
    splsmod  <- splsDA(x, grouping, usePriors=usePriors, eta=eta, K=K, ...)

#    if (length(keepInd) > 0) {
#        splskmod <- splsDA(x[,keepInd], grouping, usePriors=usePriors, eta=.99, K=K, ...)

#    } else splskmod <- NULL
    structure(list(splsda=splsmod, pvals=pmat,
                eta=eta, x=x, y=grouping), class="splsDA")
}



predict.splsDA <- function(train.mod, test) {
    fullpred <- caret::predict.splsda(train.mod$splsda, test)
    
    if (!is.null(train.mod$splsdakeep)) {
        keeppred <- caret::predict.splsda(train.mod$splsdakeep, test[,train.mod$keep])
    } else {
        keeppred <- NULL
    }
    list(Full=fullpred, Keep=keeppred)
}

splsdaboot <- function(x, statisticboot, statisticperm, R, ncpus=1, ...) {
    res     <- boot::boot(x, statisticboot, R=R, parallel="multicore", ncpus=ncpus, ...)
    null_av <- boot::boot(x, statisticperm, sim='permutation', R=R, parallel="multicore", ncpus=ncpus, ...)
    class(res) <- 'list'
    structure(c(res, list(null_av=null_av)), class='splsdaboot')
}

pval <- function(x, ...) {
    UseMethod("pval")
}

pval.splsdaboot <- function(x, sided='both', mar=2) {
# calculate 1 or 2 way pseudo p-val from boot object
# Args: a boot object
    if (sided != "both") stop("only two-sided currently supported")
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


cv.splsDA.fit <- function(x, grouping, eta, K, noise=0, fold=5, usePriors=FALSE, n.core) {
# Find eta & k by cross validation
## wrapper for spls::cv.spls, need to add some random noise so variance isn't zero
    x <- x + matrix(rnorm(prod(dim(x)), 0, noise), ncol=ncol(x))
    capture.output(out <- spls::cv.splsda(x, grouping, eta=eta, K=K, fold=fold, plot.it=FALSE, n.core=n.core))
    out
}





#cross_validation.splsDA <- function(x, grouping, k=5, rep=1, allgroups=FALSE, ...) {
# # Perform n-fold cross validation of LDA decision rule by leaving out samples and returning prediction
# #Args:
# #  method -> a method string that takes x and grouping (or X & y) and has a predict method for resulting object
# #  x -> N by p x matrix
# #  k    -> divide x into n even fractions for cross validation
# #  rep  -> repeat randomized x-validation
# # allgroups -> all groups should be represented at least twice in the subsample (could mean doubling up)
# #  ...  -> additional arguments to lda_r
# # Returns:
# #  matrix of n X rep with prediction error

#    nrec <- function(N, k) {
#        r <- N %% k
#        b <- floor(N/k)
#        c(rep(b+1, r), rep(b, k-r))
#    }

#    ind       <- 1:nrow(x)
#    err       <- vector('list', length=rep)
#    errKeep   <- vector('list', length=rep)
#    subsizes  <- nrec(length(ind), k)
#    groupfact <- as.factor(grouping)
#    for (l in 1:rep) {
#        rind     <- sample(ind)
#        mmFull <- vector('list', length=length(subsizes))
#        mmKeep <- vector('list', length=length(subsizes))
#        
#        if (allgroups) {
#            memlist <- split(1:length(grouping), grouping[rind])
#            rindmat <- suppressWarnings(do.call('rbind', memlist))
#        }
#        
#        for (i in 1:length(subsizes)) {
#            rind.i    <- ((subsizes[i]*i)-(subsizes[i]-1)):(subsizes[i]*i)
#            
#            if (allgroups) {
#                rind.i    <- rindmat[rind.i]
#            }
#            
#            test.ind  <- rind[rind.i]
#            test      <- x[test.ind,]
#            train     <- x[-test.ind,]

#            args <- list(...)

#            vars <- apply(train, 2, var)
#            zind <- which(vars <= 1e-3)
#            vind  <- setdiff(1:ncol(train), zind)

#            train.mod <- do.call(splsDA_main, c(list(train[,vind], grouping[-test.ind]), args))

#            pred  <- predict.splsDA(train.mod, test[,vind])
#            if ('class' %in% names(pred)) pred <- pred$class  # eg if method is lda_r
#            mmFull[[i]]   <- table(factor(pred$Full, levels=levels(groupfact)), 
#                                    factor(grouping[test.ind], levels=levels(groupfact)))

#            if (!is.null(pred$Keep)) {
#                mmKeep[[i]]   <- table(factor(pred$Keep, levels=levels(groupfact)),
#                                    factor(grouping[test.ind], levels=levels(groupfact)))
#            }
#        }
#        mfind <- which(unlist(lapply(lapply(mmFull, dim), 
#                    function(x) all(x == rep(length(unique(grouping)), 2) ))))
#        mmFull  <- mmFull[mfind]
#        err[[l]] <- Reduce('+', mmFull) / length(mmFull)
#        err[[l]] <- err[[l]] / rowSums(err[[l]])
#        if (!all(unlist(lapply(mmKeep, is.null)))) {
#            mfind <- which(unlist(lapply(lapply(mmKeep, dim), 
#                            function(x) all(x == rep(length(unique(grouping)), 2) ))))
#            mmKeep  <- mmKeep[mfind]
#            errKeep[[l]] <- Reduce('+', mmKeep) / length(mmKeep)
#            errKeep[[l]] <- errKeep[[l]] / rowSums(errKeep[[l]])
#        } else errKeep <- NULL
#    }
#    structure(list(errFull = err, errKeep = errKeep), class='cv')
#}


#spls.dist <- function(x, ...) {
##    means <- irispls$meanx
##    X <- scale(xobj$x, center = means, scale = FALSE)
#    scores <- as.matrix(x$x[,x$A]) %*% as.matrix(x$projection)
#    dist(scores)
#}


#spls.sil <- function(dist, grouping, ...) {
#    x <- list()
#    x$clustering <- as.integer(grouping)
#    cluster:::silhouette.default(x, dist=dist)
#}


spls2pls <- function(x) {
    pls(x$x[,x$A])
}





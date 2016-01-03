
#' stability selection of sparse models 
#'' via stars
#' @export
spls.stars <- function (x, y, fold = 10, K, eta, kappa = 0.5, select = "pls2", 
    fit = "simpls", scale.x = TRUE, scale.y = FALSE,  stars.thresh = 0.05, ncores=2, rep.num=20,
    stars.subsample.ratio = NULL) {
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ip <- c(1:p)
    y <- as.matrix(y)
    q <- ncol(y)

     if (is.null(stars.subsample.ratio)) {
         if (n > 144) 
             stars.subsample.ratio = 10 * sqrt(n)/n
         if (n <= 144) 
             stars.subsample.ratio = 0.8
     }

   mergelist <- vector('list', length(eta))
    ## compute variability
    mergeArr <- array(0, dim=c(p, length(eta), length(K)))
    variability <- matrix(0, length(eta), length(K))

   for (j in 1:length(eta)) {
         cat(paste("eta =", eta[j], "\n"))
       merge <- parallel::mclapply(1:rep.num, function(i) {
           ind.sample = sample(c(1:n), floor(n * stars.subsample.ratio), replace = FALSE)
           object <- spls::spls(x[ind.sample, , drop = FALSE], y[ind.sample, 
                    , drop = FALSE], eta = eta[j], kappa = kappa, 
                    K = max(K), select = select, fit = fit, scale.x = scale.x, 
                    scale.y = scale.y, trace = FALSE)
              return(object$betamat)
         }, mc.cores=ncores)
        mergelist[[j]] <- merge
        for (k.ind in 1:length(K)) {
            k <- K[k.ind]
            mm <- sapply(mergelist[[j]][1:rep.num], function(x) sign(abs(x[[k]][,1])))
            mergeArr[,j,k.ind] <- apply(mm, 1, function(x) sum(x)/length(x))
            variability[j,k.ind] <- 4 * sum(mergeArr[,j,k.ind] * (1 - mergeArr[,j,k.ind]))/(p)
        }
    }
    colnames(variability) <- K
    rownames(variability) <- eta
    opt.index = max(which.max(variability >= 
                stars.thresh)[1] - 1, 1)
    return(list(merge=mergeArr, variability=variability, opt.index=opt.index))
}

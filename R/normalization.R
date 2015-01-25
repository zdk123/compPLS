#####################################################################
# Different normalization schemes for microbiome counts (real or fake)
#
# @author Zachary Kurtz
# @date 10/10/2013
#####################################################################


#' @export
norm_pseudo  <- function(x) norm_to_total(x+1)


#' @keywords internal
norm_diric   <- function(x, rep=1) {
    require(VGAM)
    dmat <- rdiric(rep, x+1)
    norm_to_total(colMeans(dmat))
}

#' @export
norm_to_total <- function(x) x/sum(x)


#' The centered log-ratio transformation for
#' compositional data (not necessarily closed/normalized!)
#'
#' The clr is computed as
#'  \code{x[i]} = log (\code{x[i]} / \code{exp(mean(log(x)))})
#'
#' @title clr The centered log-ratio transformation
#' @param x a numeric data vector containing components of a composition
#' @param ... additional arguments
#' @return clr transformed \code{x}
#' @examples
#' # vector examples:
#' clr(norm_to_total(1:10))
#' clr(1:10)
#'
#' # matrix examples:
#' dmat <- matrix(exp(rnorm(110)), 10)
# rows are samples/compositions, cols are features/components
#' clr(dmat, 1) 
# cols are samples/compositions, rows are features/components
#' clr(dmat, 2) 
#' @rdname clr
#' @export
clr <- function(x, ...) {
    UseMethod('clr', x)
}



#'
#' @param base base of log to use, default is natural log
#' @param tol machine tolerance for a zero count, default is machine tol (.Machine$double.eps)
#' @rdname clr
#' @method clr default
#' @export clr.default
clr.default <- function(x, base=exp(1), tol=.Machine$double.eps) {
    nzero <- (x >= tol)
    LOG <- log(ifelse(nzero, x, 1), base)
    ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}


#'
#' 
#' @rdname clr
#' @method clr matrix
#' @export clr.matrix
clr.matrix <- function(x, mar=2, base=exp(1), tol=.Machine$double.eps) {
#    apply(x, mar, clr, ...)
    if (!mar %in% c(1,2)) stop('mar (margin) must be 1 (compositions are rows) or 2 (compositions are columns)')
    
    if (mar == 1) x <- t(x)
    nzero <- (x >= tol)
    LOG <- log(ifelse(nzero, x, 1), base)
    means <- colMeans(LOG)
    x.clr <- LOG - rep(means/colMeans(nzero), each=nrow(x))
    x.clr[!nzero] <- 0.0
    if (mar == 1) x.clr <- t(x.clr)
    x.clr
}


#'
#' 
#' @rdname clr
#' @method clr data.frame
#' @export clr.data.frame
clr.data.frame <- function(x, mar=2, ...) {
    clr(as.matrix(x), mar, ...)
}



#' The additive log-ratio transformation for
#' compositional data (not necessarily closed/normalized!)
#' 
#'  The alr transformation is computed as:
#'  \code{x[i]} = log ( \code{x[i]} /  x[D] )
#'
#' @title alr The additive log-ratio transformation
#' @param x a numeric data vector containing components of a composition
#' @param ... additional arguments
#' @return alr transformed \code{x}
#' @examples
#' # vector examples:
#' alr(norm_to_total(1:10))
#' alr(1:10)
#'
#' # matrix examples:
#' dmat <- matrix(exp(rnorm(110)), 10)
# rows are samples/compositions, cols are features/components
#' alr(dmat, 1) 
# cols are samples/compositions, rows are features/components
#' alr(dmat, 2) 
#' @rdname alr
#' @export
alr <- function(x, ...) {
    UseMethod("alr", x)
}


#' @param divcomp index of the divisor component
#' @param removeDivComp remove divisor component from the resulting data
#' @param base base of log to use, default is natural log
#' @param tol machine tolerance for a zero count, default is machine tol (.Machine$double.eps)
#' @rdname alr
#' @method alr default
#' @export alr.default
alr.default <- function(x, divcomp=1, base=exp(1), removeDivComp=TRUE,
                        tol=.Machine$double.eps) {
    zero <- (x >= tol)
    LOG <- log(ifelse(zero, x, 1), base)
    x.alr <- ifelse(zero, LOG - LOG[divcomp], 0.0)
    if (removeDivComp) x.alr[-divcomp]
    else x.alr
}

#'
#' 
#' @rdname alr
#' @method alr matrix
#' @export alr.default
alr.matrix <- function(x, mar=2, divcomp=1, base=exp(1), removeDivComp=TRUE,
                        tol=.Machine$double.eps) {
    if (mar == 1) x <- t(x)
    zero <- (x >= tol)
    LOG <- log(ifelse(zero, x, 1), base)
    x.alr <- ifelse(zero, LOG - rep(LOG[divcomp,], each=nrow(x)), 0.0)
    if (removeDivComp) x.alr <- x.alr[-divcomp,]

    if (mar ==1) t(x.alr)
    else x.alr
}

#'
#' 
#' @rdname alr
#' @method alr data.frame
#' @export alr.data.frame
alr.data.frame <- function(x, mar=2, ...) {
    alr(as.matrix(x), mar, ...)
}

#' @keywords internal
ilr <- function(x.f, V, ...) {
    UseMethod('ilr')
}

#' @keywords internal
ilr.default <- function(x.f, ...) {

}



#' @keywords internal
CSS <- function(x, ...) {
    UseMethod("CSS")
}
#' @keywords internal
CSS.default <- function(x, p=0.05, sl=1000) {
# Cumulative sum scaling Normalization Paulson et al 2013 (Nature Methods)
    xx <- x
    xx[x==0] <- NA
    qs <- quantile(xx, p=p, na.rm=TRUE)
    xx <- x - .Machine$double.eps
    normFactor <- sum(xx[xx <= qs])
    (x/normFactor)*sl
}
#' @keywords internal
CSS.matrix <- function(x, p=CSSstat(x), sl=1000, mar=2) {
    apply(x, mar, CSS, p=p, sl=sl)
}
#' @keywords internal
CSSstat <- function(mat, rel=0.1) {
    smat <- sapply(1:ncol(mat), function(i) {
        sort(mat[, i], decreasing = FALSE)
    })
    ref <- rowMeans(smat)
    yy  <- mat
    yy[yy == 0] <- NA
    ncols <- ncol(mat)
    refS  <- sort(ref)
    k     <- which(refS > 0)[1]
    lo    <- (length(refS) - k + 1)
    diffr <- sapply(1:ncols, function(i) {
            refS[k:length(refS)] - quantile(yy[, i], p = seq(0, 
                1, length.out = lo), na.rm = TRUE)})
    diffr2 <- apply(abs(diffr), 1, median, na.rm = TRUE)
    x <- which(abs(diff(diffr2))/diffr2[-1] > rel)[1]/length(diffr2)
    names(x) <- NULL
    x
}

#' @keywords internal
DESeq <- function(x, ...) {
    UseMethod("DESeq")
}
#' @keywords internal
DESeq.matrix <- function(mat, c) {
    # compute geometric mean along columns
    matt  <- mat
    matt[matt == 0] <- NA
    k_ref <- apply(matt, 1, function(x) exp(mean(log(x), na.rm=TRUE)))
    krefmat <- matrix(rep(k_ref,ncol(mat)), nrow=nrow(mat))
    s_hat   <- apply(matt/krefmat, 2, median, na.rm=TRUE)
    if (missing(c)) {
        fn <- function(c, s_hat) abs(sum(log(c*s_hat)))
        c  <- optimize(fn, interval=0:10, s_hat=s_hat)$minimum
    }
    s <- c * s_hat
    smat <- matrix(rep(s,nrow(mat)), ncol=ncol(mat), byrow=TRUE)
    mat/smat
}




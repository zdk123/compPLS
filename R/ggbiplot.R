################################################################
#  Methods for making pretty biplots with ggplot2 package
#
#  based on code by http://www.vince.vu/software/#ggbiplot


#' Pretty biplots using ggplots
#' @title methods for making biplots from various projection & classification models
#' @param xobj The object to be plotted
#' @rdname ggbiplot
#' @export
ggbiplot <- function(xobj, ...) {
    UseMethod("ggbiplot")
}

#' @rdname ggbiplot
#' @method ggbiplot princomp
#' @export ggbiplot.princomp
ggbiplot.princomp <- function(xobj, ...) {
    nobs.factor <- sqrt(xobj$n.obs)
    d <- xobj$sdev
    scores <- sweep(xobj$scores, 2, 6/(d * nobs.factor), FUN = '*')
    ggbiplot.default(list(scores=scores, loadings=xobj$loadings), ...)
}


#' @rdname ggbiplot
#' @method ggbiplot prcomp
#' @export ggbiplot.prcomp
ggbiplot.prcomp <- function(xobj, ...) {
    nobs.factor <- sqrt(nrow(xobj$x) - 1)
    d <- xobj$sdev
    scores   <- sweep(xobj$x, 2, 1/(d * nobs.factor), FUN = '*')
    loadings <- xobj$rotation
    ggbiplot.default(list(scores=scores, loadings=loadings), ...)
}


#' @rdname ggbiplot
#' @method ggbiplot lda
#' @export ggbiplot.lda
ggbiplot.lda <- function(xobj, ...) {
    xname <- xobj$call$x
    gname <- xobj$call[[3L]]
    X <- eval.parent(xname)
    g <- eval.parent(gname)
    means <- colMeans(xobj$means)
    X <- scale(X, center = means, scale = FALSE)
    x <- as.matrix(X) %*% xobj$scaling
    nobs.factor <- sqrt(nrow(x))
    d <- apply(x, 2, sd)
    x <- sweep(x, 2, 25/(d * nobs.factor), FUN = '*')
    ggbiplot.default(list(scores=x, loadings=xobj$scaling), ...)
}


#' @rdname ggbiplot
#' @method ggbiplot plsda
#' @export ggbiplot.plsda
ggbiplot.plsda <- function(xobj, Yplot=FALSE, ...) {
    if (Yplot) {
        scores <- xobj$Yscores
        loadings <- xobj$Yloadings
    } else {
        scores <- xobj$scores
        loadings <- xobj$loadings
    }
    nobs.factor <- sqrt(nrow(scores))
    d     <- apply(scores, 2, sd)
    means  <- colMeans(scores)
    scores <- scale(scores, center = means, scale = FALSE)
    scores <- sweep(scores, 2, 6/(d * nobs.factor), FUN = '*')
    ggbiplot.default(list(scores=scores, loadings=loadings), ...)

}


#' @rdname ggbiplot
#' @method ggbiplot splsda
#' @export ggbiplot.splsda
ggbiplot.splsda <- function(xobj, ...) {
#    means <- irispls$meanx
#    X <- scale(xobj$x, center = means, scale = FALSE)
    X <- xobj$x
    if (is.null(xobj$scores))
      xobj$scores <- as.matrix(X[,xobj$A]) %*% as.matrix(xobj$projection)
    xobj$loadings <- xobj$projection
    ggbiplot.plsda(xobj, ...)
}


#' @rdname ggbiplot
#' @method ggbiplot matrix
#' @export ggbiplot.matrix
ggbiplot.matrix <- function(xobj, ...) {
    scores   <- xobj
    loadings <- matrix(NA, nrow(xobj), ncol(xobj))
    ggbiplot.default(list(scores=scores, loadings=loadings), plot.loadings=FALSE, 
                          equalcoord=FALSE, ...)
}


#' @param grouping an optional grouping vector (ie - for coloring points)
#' @param select index of components to be plotted (must be length 2)
#' @param circle enclose points in a circle
#' @param circle.prob controls circle diameter (scales data std dev) if \code{circle = TRUE}
#' @param plot.loadings should loading vectors be plotted
#' @param label.loadings text of loadings labels, taken from rownames of loadings (depends on class of \code{xobj})
#' @param label.offset absolute offset for loading labels, so labels don't cover loadings vectors
#' @param scale.loadings scale length of loading vectors for plotting purposes
#' @param col.loadings a single value of vector for color of loadings
#' @param alpha controls relative transparency of various plot features
#' @param col color factor for points
#' @param group.ellipse enclose within-group points in an covariance ellipse
#' @param scale.ellipse scale \code{group.ellipse} to 1 standard deviation
#' @param group.cloud connect within-group points to a group mean point with a straight edge
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @param equalcoord equal coordinates, ie should the plot area be square?
#' @param size point size
#' @param size.loadings line width of loading vectors
#'
#' @details additional plotting attributes (eg colors, themes, etc) can be chained on in the usual way for ggplots
#' @examples
#'  # an LDA example with iris data
#' ldamod <- lda(iris[,1:4], grouping=iris[,5])
#' ggbiplot(ldamod, grouping=iris[,5], alpha=.7, group.cloud=TRUE) + theme_bw()

#' @rdname ggbiplot
#' @method ggbiplot default
#' @export ggbiplot.default
ggbiplot.default <- function(xobj, grouping, select=1:2, circle = FALSE, circle.prob = 0.69,
                     plot.loadings=TRUE, label.loadings=FALSE, label.offset=0, label.size=4.5,
                     scale.loadings = 1, col.loadings=scales::muted("red"), alpha = 1, col=grouping, 
                     group.ellipse=FALSE, scale.ellipse = 1, group.cloud = FALSE, 
                     xlab="", ylab="", equalcoord=TRUE, size=3, size.loadings=1) {
    ## get scores and loadings from xobj
    if (length(select) > 2) stop("Error: only 2d plots supported")
    if (length(select) < 2) stop("Error: need at least 2 coordinates/components")
    scores   <- data.frame(xvar=xobj$scores[,select[1]], yvar=xobj$scores[,select[2]])
    loadings <- data.frame(xvar=xobj$loadings[,select[1]], yvar=xobj$loadings[,select[2]])

    # standardize scores (?)
    # Base plot
    g <- ggplot(data = scores, aes(x = xvar, y = yvar))
    
    if (plot.loadings) {
         g <- g +
          geom_segment(data = loadings*scale.loadings,
                       aes(x = 0, y = 0, xend = xvar, yend = yvar),
                       arrow = grid::arrow(length = grid::unit(1/2, 'picas')),
                       color = col.loadings, size = size.loadings)
    }
    
    if (label.loadings) {
        labs <- rownames(loadings)
        # compute angles from orig.
        ang <- atan2(loadings$yvar*scale.loadings, loadings$xvar*scale.loadings)
        hyp <- sqrt((loadings$yvar*scale.loadings)^2 + (loadings$xvar*scale.loadings)^2)
        
        labdat <- data.frame(newx=(hyp + label.offset)*cos(ang),
                             newy=(hyp + label.offset)*sin(ang),
                             label=labs)
        g <- g +
           geom_text(aes(x=newx, y=newy, label=label), data=labdat, size=label.size)
    }
    
    if (!missing(grouping)) {
        gind <- order(grouping)
        grouping <- grouping[gind]
        scores   <- scores[gind,]
        g <- g + geom_point(data = data.frame(scores, grouping=grouping),
                       aes(color = grouping), alpha = alpha, size=size) 

    } else {
        if (!missing(col)) {
            g <- g + geom_point(data=data.frame(scores, col), aes(color = col), 
                                alpha = alpha, size=size)
                     
        } else {
        g <- g + geom_point(alpha = alpha, size=size)
        }
    }
    
    
    if (group.ellipse && !missing(grouping)) {
        l   <- 200
        group.scores  <- split(scores[,1:2], grouping)
        group.centers <- lapply(group.scores, colMeans)
        group.cov     <- lapply(group.scores, cov)
        group.RR      <- lapply(group.cov, chol)
        angles   <- seq(0, 2*pi, length.out=l)
        ell.list <- lapply(group.RR, function(RR) 
                     scale.ellipse * cbind(cos(angles), sin(angles)) %*% RR)
        ellCntr  <- lapply(1:length(ell.list), function(i)
                        sweep(ell.list[[i]], 2, group.centers[[i]], "+"))
        names(ellCntr) <- names(ell.list)
        ell.df   <- as.data.frame(do.call("rbind", ellCntr))
        ell.df$grouping <- factor(rep(names(ellCntr), each=l), levels=names(ellCntr))
        g <- g + geom_path(data = ell.df, aes(color = grouping, group = grouping))
    }
    
    if (group.cloud && !missing(grouping)) {
        group.scores  <- split(scores[,1:2], grouping)
        group.centers <- lapply(group.scores, colMeans)
        centers.df    <- do.call('rbind', rep(group.centers, table(grouping)))
        rownames(centers.df) <- rownames(scores)
        colnames(centers.df) <- c("xcntr", "ycntr")
        loadCntr.df <- cbind(scores, centers.df, grouping)
        g <- g +
          geom_segment(data = loadCntr.df,
                       aes(x = xcntr, y = ycntr, xend = xvar, yend = yvar,
                       color = grouping), alpha=10^log(alpha/1.4))
    }

    
    if (circle) {
    # scale circle radius
        r1 <- sqrt(qchisq(circle.prob, df = 2)) * max(scores$xvar^2)^(1/2)
        r2 <- sqrt(qchisq(circle.prob, df = 2)) * max(scores$yvar^2)^(1/2)

        theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
        circdat <- data.frame(xvar = r1 * cos(theta), yvar = r2 * sin(theta))
        g <- g + geom_path(aes(x=xvar, y=yvar), data = circdat, color = scales::muted('black'), 
                           size = 0.5, alpha = alpha/3)
    }

    if (equalcoord) {
        if (circle) {
            xrange <- range(circdat$xvar)
            yrange <- range(circdat$yvar)
        } else {
            xrange <- c(-max(abs(scores$xvar)), max(abs(scores$xvar)))
            yrange <- c(-max(abs(scores$yvar)), max(abs(scores$yvar)))
        }
        
        g <- g + scale_x_continuous(xlab, limits=xrange) +
                 scale_y_continuous(ylab, limits=yrange)
    } else {
        g <- g + scale_x_continuous(xlab) +
                 scale_y_continuous(ylab)
    }
    g
}

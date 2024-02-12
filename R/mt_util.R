
#' ========================================================================
#' lwc-29-03-2013: PCA outlier plot by lattice
#' wll-29-11-2015: Examples of 'panel.elli' and 'panel.outl' give more
#' general information about ellipses and outliers. If you *only* want to
#' plot outliers based on PCA in a general way, for example outliers in
#' different groups or in conditional panel, you can write an wrapper
#' function combining with 'pca.comp', 'panel.elli' and 'panel.oult'. The
#' example is 'pca_plot_wrap'.
pca.outlier <- function(x, center = TRUE, scale = TRUE,
                        conf.level = 0.975, ...) {

  #' lwc-29-03-2013: Lattice panel for plotting outliers with ellipse
  #' lwc-03-04-2013: To avoid error: formal argument "***" matched by
  #' multiple actual arguments, use: ..., type,col, lty, lwd.
  #' wll-29-11-2015: More general panel function for outlier 'panel.outl'
  panel.outlier <- function(x, y, groups = NULL, elli, labs, id, ...,
                            type, col, lty, lwd) {
    #' dots <- list(...)
    if (is.null(groups)) {
      panel.xyplot(x, y, ...)
    } else {
      panel.superpose(x, y, groups, ...)
    }

    panel.abline(h = 0, v = 0, col = c("gray"), lty = 2)
    #' overall ellipse line
    panel.points(elli[, 1], elli[, 2], type = "l", col = "red", lwd = 2, ...)
    #' labelling outliers
    if (any(id)) {
      ltext(x[id], y[id], labs[id], ...)
      #' cex = dots$cex, adj = dots$adj)
    }
  }

  #' argument list
  dots <- list(...)
  if (length(dots) > 0) {
    args <- dots
  } else {
    args <- list()
  }

  #' calculate PCA
  pcs <- 1:2 #' only use PC1 and PC2
  pca <- prcomp(x, center = center, scale. = scale) #' strip off dot arguments
  vars <- pca$sdev^2
  vars <- vars / sum(vars) #' Proportion of Variance
  names(vars) <- colnames(pca$rotation)
  vars <- round(vars * 100, 2)
  dfn <- paste(names(vars), " (", vars[names(vars)], "%)", sep = "")
  x <- data.frame(pca$x)
  names(x) <- dfn
  x <- x[, pcs]

  #' outlier detection by Mahalanobis distances
  cent <- colMeans(x)
  cova <- cov(x)
  dist <- sqrt(mahalanobis(x, center = cent, cov = cova))
  cuto <- sqrt(qchisq(conf.level, ncol(x)))
  id <- dist > cuto

  #' get ellipse point
  elli <- ellipse(var(x), centre = cent, level = conf.level)

  #' handle args
  labs <- rownames(x)
  args <- c(list(x = x[, 2] ~ x[, 1], data = x), args)

  if (is.null(args$xlab)) args$xlab <- names(x)[1]
  if (is.null(args$ylab)) args$ylab <- names(x)[2]
  if (F) {
    xlim <- c(min(x[, 1], elli[, 1]), max(x[, 1], elli[, 1]))
    ylim <- c(min(x[, 2], elli[, 2]), max(x[, 2], elli[, 2]))
    if (is.null(args$xlim)) args$xlim <- xlim
    if (is.null(args$ylim)) args$ylim <- ylim
  }

  args <- c(args, panel = panel.outlier)

  #' arguments for panel.outlier
  args$elli <- elli
  args$labs <- labs
  args$id <- id

  p <- do.call("xyplot", args)

  ret <- list(
    plot = p, outlier = which(id), conf.level = conf.level,
    mah.dist = dist, cutoff = cuto
  )
  return(ret)
}

#' ========================================================================
#' lwc-03-06-2010: PCA plot with outlier detection
#' lwc-01-09-2010: Add group info.
#' To-Do:
#'  1.) Display group text inside the ellipse
pca.outlier.1 <- function(x, center = TRUE, scale = TRUE, conf.level = 0.975,
                          group = NULL, main = "PCA", cex = 0.7, ...) {
  #' calculate PCA
  pcs <- 1:2 #' only use PC1 and PC2
  pca <- prcomp(x, center = center, scale. = scale) #' strip off dot arguments
  vars <- pca$sdev^2
  vars <- vars / sum(vars) #' Proportion of Variance
  names(vars) <- colnames(pca$rotation)
  vars <- round(vars * 100, 2)
  dfn <- paste(names(vars), " (", vars[names(vars)], "%)", sep = "")
  x <- data.frame(pca$x)
  names(x) <- dfn
  x <- x[, pcs]

  #' outlier detection by Mahalanobis distances
  cent <- colMeans(x)
  cova <- cov(x)
  dis <- sqrt(mahalanobis(x, center = cent, cov = cova))
  cutoff <- sqrt(qchisq(conf.level, ncol(x)))
  outlier <- which(dis > cutoff)

  #' Plot PCA with ellipse and outliers
  z <- ellipse(x = cova, center = cent, level = conf.level)
  x1 <- c(min(x[, 1], z[, 1]), max(x[, 1], z[, 1]))
  y1 <- c(min(x[, 2], z[, 2]), max(x[, 2], z[, 2]))

  if (is.null(group)) {
    plot(x, xlim = x1, ylim = y1, main = main, ...)
  } else {
    col <- unclass(group)
    pch <- unclass(group)
    plot(x, xlim = x1, ylim = y1, main = main, col = col, pch = pch, ...)
    legend("topright",
      legend = sort(unique(levels(group))),
      cex = cex, col = sort(as.vector(unique(col))),
      pch = sort(as.vector(unique(pch)))
    )
  }
  lines(z, type = "l", col = "red") #' plot ellipse

  #' plot the outliers
  if (length(outlier) > 0) {
    xrange <- par("usr")
    xrange <- xrange[2] - xrange[1] #' control offset of text position

    #' display names
    txt <- names(dis[outlier])
    if (is.null(txt)) txt <- outlier
    text(x[outlier, 1], x[outlier, 2] + xrange / 50, txt,
      col = "blue",
      cex = cex
    )
  }

  ret <- list(
    outlier = outlier, conf.level = conf.level, mah.dist = dis,
    cutoff = cutoff
  )
  return(ret)
}

#' =======================================================================
#' lwc-13-11-2007: Plot columns of matrix-like object by group
#' lwc-18-12-2007: Major changes. Call .grpplot.
#' lwc-28-09-2008: Add ocall. For details, see plot.shingle in lattice
#' lwc-11-02-2010: Change the point of median in boxplot as line
#' lwc-21-02-2010: change name from gplot to grpplot.
#' lwc-15-07-2015: remove 'ep'
#' Note: Some examples of auto.key:
#'         auto.key=list(columns=nlevels(x.1$.y)),
#'         auto.key = list(space = "right"),
#' Usages
#'  data(iris)
#'  x <- iris[,1:4]
#'  y <- iris[,5]
#'  grpplot(x[,1:2],y, scale=T, pcs=c(2,1),ep=2)
#'  grpplot(x,y, scale=T, pcs=c(2,1),ep=1)
grpplot <- function(x, y, plot = "pairs", ...) {
  ocall <- sys.call(sys.parent())
  ocall[[1]] <- quote(grpplot)
  x <- as.data.frame(x)
  #' lwc-19-12-07: data.frame(x) will change names of columns in some
  #' situations.
  y <- factor(y)

  plot <- match.arg(plot, c("strip", "box", "density", "pairs"))

  #' reshape data set for "strip", "boxplot" and "density".
  x.1 <- stack(x)
  x.1$ind <- factor(x.1$ind, levels = unique.default(x.1$ind))
  #' lwc-21-11-2007: Original ind of stack is a sorted-level factor. Lattice
  #'   will use this factor level to arrange panel order. To be consist with
  #'   feature rank descent order, the factor levels are ordered by the
  #'   feature rank from to to bottom. Therefore, no sort used inside factor
  #'   function.
  x.1$.y <- rep(y, ncol(x))

  grpplot <-
    switch(tolower(plot),
      strip = stripplot(values ~ .y | ind,
        data = x.1, as.table = T, ,
        ...
      ),
      box = bwplot(values ~ .y | ind,
        data = x.1, as.table = T,
        pch = "|", ...
      ),
      density = densityplot(~ values | ind,
        data = x.1,
        groups = x.1$.y, plot.points = F, #' kern = "rect",
        as.table = T, ...
      ),
      pairs = .grpplot(x, y, ...)
    )
  grpplot$call <- ocall
  return(grpplot)
}

#' =======================================================================
#' lwc-13-12-2007: Scatter plot by group
#' lwc-17-12-2007: deal with 1 column problem. This idea is from
#' ldahist
#' lwc-09-01-2008: argument of default of auto.key.
#' lwc-12-01-2008-note:
#'   Colours in supervose.symbol re-cycle after 7 by default. I did not
#'   change the default number inside the function by: par.settings =
#'   list(superpose.symbol=list(pch=1:nlevels(y), col=1:nlevels(y))),
#'   For convenient, user can change the color scheme outside .grpplot, such as
#'       superpose.symbol <- trellis.par.get("superpose.symbol")
#'       superpose.symbol$col <- rainbow(16)
#'       trellis.par.set("superpose.symbol",superpose.symbol)
#'  Then call .grpplot.  A list of colour's names is produced by colors().
#' lwc-17-01-2008: pch range from 1 to 25. So I recycle the symbols if number
#'     of symbol exceed the limits.
#' TO-DO: 1.) Check validity of arguments
#'        2.) Better way to process ep(Currently, ep will be 0,1 and 2).
#'        3.) panel.xyplot doesn't know anything about xlab, ylab, etc., and
#'           you can specify cex as part of the top level call. (claimed by
#'           Deepayan Sarkar, author of package lattice). So
#'           trellis.par.get(), trellis.par.set() or par.settings will be
#'           help for global or local parameters setting.
#' lwc-15-07-2015: remove 'ep' and call panel.elli.1
.grpplot <- function(x, y, auto.key = list(space = "right"),
                     par.settings = list(superpose.symbol = list(pch = rep(1:25))),
                     xlab, ylab, ...) {
  ocall <- sys.call(sys.parent())
  ocall[[1]] <- quote(.grpplot)
  if (!is.matrix(x) && !is.data.frame(x)) {
    x <- as.matrix(x)
  }
  y <- factor(y)

  if (ncol(x) == 2) {
    if (missing(xlab)) xlab <- names(x)[2]
    if (missing(ylab)) ylab <- names(x)[1]
    p <-
      xyplot(x[, 1] ~ x[, 2],
        groups = y, as.table = T,
        #' xlab=names(x)[2], ylab=names(x)[1],  #' lwc-07-08-14
        xlab = xlab, ylab = ylab,
        auto.key = auto.key,
        par.settings = par.settings,

        #' par.settings =
        #'   list(superpose.symbol=list(pch = rep(1:25, len = nlevels(y)))),

        scales = list(cex = 0.8), #' for axis font
        panel = function(x, y, ...) {
          panel.xyplot(x, y, ...)
          panel.elli.1(x, y, ...)
          #' panel.outl(x,y, ...)
        }, ...
      )
  } else if (ncol(x) > 2) {
    p <-
      splom(~x,
        groups = y, as.table = T, xlab = "",
        auto.key = auto.key,
        par.settings = par.settings,

        #' par.settings =
        #'   list(superpose.symbol=list(pch = rep(1:25, len = nlevels(y))),
        #'        axis.text=list(cex=0.7)),

        #' varname.cex = 1.0, cex=0.6, #' pscales = 0,
        panel = function(x, y, ...) {
          panel.xyplot(x, y, ...)
          panel.elli.1(x, y, ...)
          #' panel.outl(x,y, ...)
        }, ...
      )
  } else {
    p <- stripplot(x[, 1] ~ y, groups = y, as.table = T,
                   ylab = colnames(x)[1], ...)
    #' p <- stripplot(x ~ y, groups=y, as.table=T, ylab= "", ...)
    #' p <- densityplot(~ x, groups = y, as.table=T,
    #'                  auto.key = auto.key,
    #'                  par.settings = list(superpose.line = list(lty=c(1:7))),
    #'                  plot.points = FALSE, ref = TRUE,...)
    if (F) {
      p <- histogram(~ x | y,
        xlab = "", type = "density",
        layout = c(1, nlevels(y)),
        panel = function(x, ...) {
          panel.histogram(x, ...)
          panel.mathdensity(
            dmath = dnorm, col = "black",
            args = list(mean = mean(x), sd = sd(x))
          )
        },
        ...
      )
    }
  }
  p$call <- ocall
  p
}

#' ========================================================================
#' lwc-12-13-2007: plot PCA using lattice package
#' lwc-15-07-2015: remove 'ep'
#' TO-DO: 1). check validity of PCs used for plotting.
#' Usage:
#'   data(iris)
#'   x <- iris[,1:4]
#'   y <- iris[,5]
#'   pcaplot(x,y, scale=T, pcs=c(2,1),ep=2)
pcaplot <- function(x, y, scale = TRUE, pcs = 1:2, ...) {
  #' pca  <- prcomp(x, scale.=scale, ...)
  pca <- prcomp(x, scale. = scale)
  vars <- pca$sdev^2
  vars <- vars / sum(vars) #' Proportion of Variance
  names(vars) <- colnames(pca$rotation)
  vars <- round(vars * 100, 2)
  dfn <- paste(names(vars), " (", vars[names(vars)], "%)", sep = "")

  x <- data.frame(pca$x)
  names(x) <- dfn
  x <- x[, pcs]

  p <- grpplot(x, y, plot = "pairs", ...)
  p
}

#' =======================================================================
#' lwc-15-07-2015: ellipse panel function which support individual and
#'  combined group plotting. It is the extension of panel.elli.
#' Usage: Under ep=2, there are three options to plot ellipse.
#'        com.grp: control which combination of groups to be plotted.
#'        no.grp:  control which individual group not to be plotted. Note
#'                 it will be overridden by com.grp.
#'        If no com.grp and no.grp, the each individual group ellipse should
#'        be plotted.
panel.elli.1 <- function(x, y, subscripts, groups = NULL, conf.level = 0.975,
                         ep = 0, com.grp = NULL, no.grp = NULL,
                         ell.grp = NULL, ...) {
  plot.elli <- function(x, y, ...) { #' plot ellipse
    Var <- var(cbind(x, y))
    Mean <- cbind(mean(x), mean(y))
    Elli <- ellipse(Var, centre = Mean, level = conf.level)
    #' panel.xyplot(x, y,...)
    #' panel.xyplot(Elli[,1], Elli[,2],...)
    panel.points(Elli[, 1], Elli[, 2], ...)
  }

  #' lwc-14-07-2015: do NOT plot x and y inside the panel function
  if (FALSE) {
    if (!is.null(groups)) {
      panel.superpose(x, y, subscripts, groups, ...)
    } else {
      panel.xyplot(x, y, ...)
    }
    panel.abline(h = 0, v = 0, col = c("gray"), lty = 2)
  }

  if (!is.null(ell.grp)) { #' ellipse based on other group info
    grp <- ell.grp[subscripts]
    tmp <- data.frame(x = x, y = y, grp = grp)

    #' wl-05-11-2021, Fri: use base R 'by'
    by(tmp, tmp$grp, function(x) {
      plot.elli(x$x, x$y, ..., type = "l", lty = 2, col = "cyan")
    })
    ## plyr::ddply(tmp, .(grp), function(x) {
    ##   plot.elli(x$x, x$y, ..., type = "l", lty = 2, col = "cyan")
    ## })

  } else if (ep == 1) { #' over-line ellipse
    plot.elli(x, y, type = "l", col = "red", ...) #' lwd=2
    #' ellipse based on groups, individual or combination.
  } else if (ep == 2) { #' plot group ellipse
    if (!is.null(com.grp)) { #' plot combined groups
      grp <- groups[subscripts]
      for (i in names(com.grp)) {
        id <- grp %in% com.grp[[i]]
        plot.elli(x[id], y[id], ..., type = "l", col = "gray")
      }
    } else if (!is.null(no.grp)) { #' plot remained groups
      grp <- groups[subscripts]
      for (i in levels(grp)) {
        id <- i == grp
        if (!(i %in% no.grp)) {
          plot.elli(x[id], y[id], ..., type = "l", col = "gray")
        }
      }
    } else { #' plot all groups
      panel.superpose(x, y, subscripts, groups, ...,
        panel.groups = plot.elli,
        type = "l", lty = 2
      )
    }
  }
}

#' =======================================================================
#' lwc-02-04-2013: Panel function for plotting ellipse  used by lattice.
#' Note: For details, see panel.ellipse of package latticeExtra
panel.elli <- function(x, y, groups = NULL, conf.level = 0.975, ...) {
  if (!is.null(groups)) {
    panel.superpose(
      x = x, y = y, groups = groups, conf.level = conf.level,
      panel.groups = panel.elli, ...
    )
  } else {
    Var <- var(cbind(x, y))
    Mean <- cbind(mean(x), mean(y))
    Elli <- ellipse(Var, centre = Mean, level = conf.level)
    panel.xyplot(Elli[, 1], Elli[, 2], ...)
  }
}

#' =======================================================================
#' lwc-02-04-2013: Panel function for plotting outliers with ellipse in
#' lattice
#' To-Do: How to keep colour of text as consistent with groups?
panel.outl <- function(x, y, subscripts, groups = NULL,
                       conf.level = 0.975, labs, ...) {
  if (!is.null(groups)) {
    panel.superpose(
      x = x, y = y, groups = groups, subscripts = subscripts,
      conf.level = conf.level, labs = labs,
      panel.groups = panel.outl, ...
    )
  } else {
    #' outlier detection by Mahalanobis distances
    mat <- cbind(x, y)
    #' row.names(mat) <- labs[subscripts]
    cent <- colMeans(mat)
    cova <- cov(mat)
    dist <- sqrt(mahalanobis(mat, center = cent, cov = cova))
    cuto <- sqrt(qchisq(conf.level, ncol(mat)))
    id <- dist > cuto

    if (any(id)) {
      panel.text(x[id], y[id], labs[subscripts][id], ...)
      #' ltext(x[id], y[id], labs[subscripts][id],...)
      #' from lattice: panel.text <- function(...) ltext(...)
    }
  }
}

#' =======================================================================
#' lwc-30-07-2013: group stats of column of a matrix/data frame including
#'  fold changes, auc and p-values.
#' lwc-08-01-2014: tidy up and minor changes.
#' wll-24-11-2015: add the adjusted p-values. Beware that 'stats.vec' has no
#' such thing.
stats.mat <- function(x, y, method = "mean", test.method = "wilcox.test",
                      padj.method = "fdr", fc = TRUE, ...) {
  #' function for calculation based on column vector.
  x <- as.data.frame(x, stringsAsFactors = F)
  res <- t(sapply(x, function(i) stats.vec(i, y, method, test.method, fc, ...)))
  res <- as.data.frame(res, stringsAsFactors = FALSE)

  #' get adjusted p-values
  padj <- round(p.adjust(res$pval, method = padj.method), digits = 4)
  res <- cbind(res, padj) #' or res <- data.frame(res, padj)

  return(res)
}

#' =======================================================================
#' lwc-30-07-2013: group stats for vector
#' lwc-09-01-2014: lack of error handling, such as limits of 'method' and
#' two groups.
#' wll-11-08-2014: add overall mean
#' wll-24-11-2015: add adjusted p-values
#' wll-01-12-2015: add an argument for fold-change. fc is only for positive
#'   values of 'x'. If not, the results are useless.
#' wll-26-01-2016: drop off change direction so the results are numeric, not
#'   the character. Note that the fold change indicates the changing
#'   direction.
stats.vec <- function(x, y, method = "mean", test.method = "wilcox.test",
                      fc = TRUE, ...) {
  #' overall mean
  omn <- do.call(method, list(x, na.rm = TRUE))
  names(omn) <- method
  #' group mean
  gmn <- tapply(x, y, method, na.rm = TRUE)
  names(gmn) <- paste(names(gmn), method, sep = ".")

  auc <- round(cl.auc(x, y), digits = 2)
  p.val <- round(.pval(x, y, test.method = test.method, ...), digits = 4)
  #' p.val  <- wilcox.test(x ~ y,correct = FALSE)$"p.value"

  if (F) {
    direc <- if (gmn[1] > gmn[2]) {
      "Down"
    } else if (gmn[1] < gmn[2]) {
      "Up"
    } else {
      "No change"
    }
  }

  if (fc) {
    fc <- round(.foldchange(gmn[2], gmn[1]), digits = 2)
    #' lwc-23-08-2013: gmn[1] is baseline
    names(fc) <- NULL
    log2.fc <- round(.foldchange2logratio(fc), digits = 2)

    res <- c(omn, gmn, #' direction=direc,
      fold.change = fc, log2.fold.change = log2.fc, auc = auc, pval = p.val
    )
  } else {
    res <- c(omn, gmn, #' direction=direc,
      auc = auc, pval = p.val
    )
  }
  return(res)
}

#' =======================================================================
#' lwc-30-07-2013: wrapper functions for p-values from test
.pval <- function(x, y, test.method = "oneway.test", ...) {
  test.method <- if (is.function(test.method)) {
    test.method
  } else if (is.character(test.method)) {
    get(test.method)
  } else {
    eval(test.method)
  }
  pval <- test.method(x ~ y, ...)$p.value
  return(pval)
}

#' =======================================================================
#' lwc-16-07-2013: Fold change from gtools
#' Usage:
#' a <- 1:21
#' b <- 21:1
#' f <- .foldchange(a, b)
#' cbind(a, b, f)
.foldchange <- function(num, denom) {
  ifelse(num >= denom, num / denom, -denom / num)
}

.foldchange2logratio <- function(foldchange, base = 2) {
  retval <- ifelse(foldchange < 0, 1 / -foldchange, foldchange)
  retval <- log(retval, base)
  retval
}

.logratio2foldchange <- function(logratio, base = 2) {
  retval <- base^ (logratio)
  retval <- ifelse(retval < 1, -1 / retval, retval)
  retval
}

#' ========================================================================
#' lwc-24-08-2011: Summary function for data vector.
#' lwc-26-08-2011: add error checking (All NAs)
#' Usage:
#'   x <- iris[,1]
#'   vec.summ.1(x)
vec.summ.1 <- function(x) {
  if (sum(!is.na(x)) < 2) {
    #' if (all(is.na(x))) {
    mean <- median <- sd <- iqr <- CI.L <- CI.H <- NA
  } else {
    mean <- mean(x, na.rm = T)
    median <- median(x, na.rm = T)
    sd <- sd(x, na.rm = T)
    iqr <- IQR(x, na.rm = T)
    conf <- t.test(x)$conf.int
    CI.L <- conf[1]
    CI.H <- conf[2]
  }

  res <- c(
    N = sum(!is.na(x)), Mean = mean, Median = median,
    "95% CI.l" = CI.L, "95% CI.u" = CI.H,
    IQR = iqr, Std = sd
  )

  #' res <- format(res,digits=3)
  res <- round(res, digits = 3)
  return(res)
}

#' =======================================================================
#' lwc-03-03-2010: Summary function for vector data
#' lwc-11-11-2011: Change Nval as N.
vec.summ <- function(x) {
  res <- c(
    N = sum(!is.na(x)), Min = min(x, na.rm = T),
    Mean = mean(x, na.rm = T),
    Median = median(x, na.rm = T),
    Max = max(x, na.rm = T),
    Std = sd(x, na.rm = T)
  )
  res <- round(res, digits = 3)
  return(res)
}

#' =======================================================================
#' lwc-03-03-2010: Summary function for data frame/matrix by column.
#' lwc-24-08-2011: Summary function for data matrix (wrapper function of
#' vec.summ).
#' lwc-22-05-2013: add dots for method's arguments.
#' Usage:
#'  data(abr1)
#'  dat <- (abr1$pos)[,110:150]
#'  dat <- mv.zene(dat)
#'  summ   <- df.summ(dat, method=vec.summ)
#'  summ.1 <- df.summ(dat, method=vec.summ.1)
df.summ <- function(dat, method = vec.summ, ...) {
  method <-
    if (is.function(method)) {
      method
    } else if (is.character(method)) {
      get(method)
    } else {
      eval(method)
    }

  dat <- as.data.frame(dat, stringsAsFactors = F)
  #' only numeric, not categorical data
  dat <- Filter(is.numeric, dat)
  #' lwc-11-10-2011: dat must be data frame here.

  res <- t(sapply(dat, function(i) method(i, ...)))

  res <- as.data.frame(res, stringsAsFactors = FALSE)
  #' res <- cbind(Variable=rownames(res),res)
  #' Do we need to add one column?
  return(res)
}

#' =========================================================================
#' lwc-11-10-2011: replace zero/negative with NA.
mv.zene <- function(dat) {
  vec.func <- function(x) {
    x <- ifelse(x < .Machine$double.eps, NA, x)
  }

  dat <- as.data.frame(dat, stringsAsFactors = F)
  res <- sapply(dat, function(i) vec.func(i))
  return(res)
}

#' ========================================================================
#' lwc-23-04-2010: Fill the zero/NA values by the mean of vector.
#' lwc-15-06-2011: minor changes.
#' lwc-22-06-2011: Replace ifelse(x < .Machine$double.eps, m, x) with
#'   ifelse(x < .Machine$double.eps, NA, x) and change its line position.
#' wl-27-11-2021, Sat: bug. if dat have minus values, do not use 'ze_ne = T'
#' Usage
#'   data(abr1)
#'   dat <- abr1$pos[,1970:1980]
#'   dat.1 <- mv.fill(dat,method="mean",ze_ne = TRUE)
mv.fill <- function(dat, method = "mean", ze_ne = FALSE) {
  method <-
    if (is.function(method)) {
      method
    } else if (is.character(method)) {
      get(method)
    } else {
      eval(method)
    }

  vec.func <- function(x) {
    if (ze_ne) {
      x <- ifelse(x < .Machine$double.eps, NA, x)
      #' vectorisation of ifelse
    }
    m <- method(x, na.rm = TRUE)

    x[is.na(x)] <- m
    #' 10-10-2011: more general for multiple filing points
    x
    #' x <- ifelse(is.na(x), m, x)  #' for missing values
    #' x <- ifelse(is.nan(x), m, x)  #' for missing values
  }

  dat <- as.data.frame(dat, stringsAsFactors = F)
  res <- sapply(dat, function(i) vec.func(i))
  return(res)
}

#' =======================================================================
#' lwc-07-09-2010: Test codes for missing values processing
#' wll-11-12-2007: Statistics and plot for missing values
#' lwc-03-03-2010: Get number of missing values by column.
#' lwc-15-06-2011: re-write
#' lwc-10-10-2011: major change. remove stats based on row.
#' Usage:
#'  data(metaboliteData, package="pcaMethods")
#'  dat <- t(metaboliteData)
#'  colnames(dat) <- paste("V", 1:ncol(dat), sep="")
#'  cls <- rownames(dat)
#'  cls <- sub("[.].*", "", cls)
#'  cls <- factor(cls)
#'  tmp <- mv.stats(dat, grp=cls)
#'  tmp <- mv.fill(dat, method = "median")
#'  tmp <- mt:::mv.pattern(dat)
mv.stats <- function(dat, grp = NULL, ...) {
  #' overall missing values rate
  mv.all <- sum(is.na(as.matrix(dat))) / length(as.matrix(dat))

  #' MV stats function for vector
  vec.func <-
    function(x) round(sum(is.na(x) | is.nan(x)) / length(x), digits = 3)
  #' vec.func  <- function(x)  sum(is.na(x)|is.nan(x)) #' number of MV
  #' sum(is.na(x)|is.nan(x)|(x==0))

  #' get number of Na, NaN and zero in each of feature/variable
  #' mv.rep <- apply(dat, 1, vec.func)
  mv.var <- apply(dat, 2, vec.func)

  ret <- list(mv.overall = mv.all, mv.var = mv.var)

  if (!is.null(grp)) {
    #' MV rate with respect of variables and class info
    mv.grp <- sapply(levels(grp), function(y) {
      idx <- (grp == y)
      mat <- dat[idx, ]
      mv <- apply(mat, 2, vec.func)
    })

    #' lwc-10-10-2011: Use aggregate. Beware that values pased in the
    #' function is vector(columns).
    #' mv.grp <- aggregate(dat, list(cls), vec.func)
    #' rownames(mv.grp) <- mv.grp[,1]
    #' mv.grp <- mv.grp[,-1]
    #' mv.grp <- as.data.frame(t(mv.grp),stringsAsFactors=F)

    #' reshape matrix for lattice
    mv.grp.1 <- data.frame(mv.grp)
    mv.grp.1$all <- mv.var #' Combine all

    var <- rep(1:nrow(mv.grp.1), ncol(mv.grp.1))
    mv.grp.1 <- stack(mv.grp.1)
    mv.grp.1$ind <- factor(mv.grp.1$ind,
      levels = unique.default(mv.grp.1$ind)
    )
    mv.grp.1$var <- var

    mv.grp.plot <-
      xyplot(values ~ var | ind,
        data = mv.grp.1, groups = mv.grp.1$ind, as.table = T,
        layout = c(1, nlevels(mv.grp.1$ind)), type = "l",
        auto.key = list(space = "right"),
        #' main="Missing Values Percentage With Respect of Variables",
        xlab = "Index of variables", ylab = "Percentage of missing values",
        ...
      )

    ret$mv.grp <- mv.grp
    ret$mv.grp.plot <- mv.grp.plot
  }

  ret
}

#' ========================================================================
#' wll-05-12-2007: Calculate the pattern of missing values.
#' Value:
#'   A matrix with (nrow(x)+1, ncol(x)+1) dimension. Except the last row and
#'   column, each row corresponds to a missing data pattern
#'   (1=observed, 0=missing). The row names shows the number of pattern.
#'   The last row contains the number of missing values
#'   with respect to each column and the last column represent the counts of
#'   each row.
#' See Also:
#'   md.pattern in package mice and prelim.norm in package norm.
#' NOTE: 1.The motivation of the function is that Ted Harding mentioned that
#'       that prelim.norm can only encode NA-patterns in an R integer for up
#'       to 31 columns. More than that, and it will not work properly or at
#'       all. (http://article.gmane.org/gmane.comp.lang.r.general/55185).
#'       Function md.pattern has also this problem since it modified from
#'       prelim.norm. 2. The function is not sorted at current stage.
#' Usage:
#'   library(mice)
#'   data(nhanes)
#'   md.pattern(nhanes)     #' from mice
#'   mv.pattern(nhanes)
mv.pattern <- function(x) {
  "%all.==%" <-
    function(a, b) apply(b, 2, function(x) apply(t(a) == x, 2, all))

  if (!(is.matrix(x) | is.data.frame(x))) {
    stop("Data should be a matrix or data frame")
  }

  #' get the pattern of missing values
  mat <- 1 * !is.na(x)
  pattern <- unique(mat)
  counts <- colSums(mat %all.==% t(unique(mat)))
  rownames(pattern) <- counts

  #' number of missing values with respect to column (variable)
  nmis <- apply(1 * is.na(x), 2, sum)
  #' number of missing values in the pattern
  pmis <- ncol(pattern) - apply(pattern, 1, sum)

  pattern <- rbind(pattern, c(nmis)) #' a trick to take off the row name
  pattern <- cbind(pattern, c(pmis, sum(nmis)))
  pattern
}

#' =========================================================================
#' lwc-24-11-2010: Get heatmap colours
#' Note: compare this:
#'  col.regions = colorRampPalette(c("green", "black", "red"))
#'  in lattice levelplot
hm.cols <- function(low = "green", high = "red", n = 123) {
  low <- col2rgb(low) / 255
  if (is.character(high)) {
    high <- col2rgb(high) / 255
  }
  col <- rgb(
    seq(low[1], high[1], len = n),
    seq(low[2], high[2], len = n),
    seq(low[3], high[3], len = n)
  )
  return(col)
}

#' =========================================================================
#' wll-23-10-2008: Wrapper function for plotting PCA. The first two PCs are
#'                 fixed in this routine.
#' wll-12-01-2008: add dot arguments for lattice plot arguments
#' Arguments:
#'   data.list - A two-layer list structure, in which the second layer
#'               include a data frame and a factor of class label. It should
#'               be noted the names of the first layer of data.list must be
#'               given.
#'   title     - A part of title string for plotting
pca_plot_wrap <- function(data.list, title = "plotting", ...) {
  if (is.null(names(data.list))) {
    names(data.list) <-
      paste(deparse(substitute(data.list)), 1:length(data.list), sep = ":")
  }
  dn <- names(data.list)
  pca <- lapply(dn, function(x) {
    pca <- pca.comp(data.list[[x]]$dat, scale = F, pcs = 1:2)
    scores <- cbind(pca$scores, cls = data.list[[x]]$cls, type = x)
    list(scores = scores, vars = pca$vars)
  })
  names(pca) <- dn

  pca.scores <- do.call(rbind, lapply(pca, function(x) x$scores))
  pca.vars <- do.call(rbind, lapply(pca, function(x) x$vars))

  pca.p <-
    xyplot(PC2 ~ PC1 | type,
      data = pca.scores, groups = pca.scores$cls, as.table = T,
      xlab = "PC1", ylab = "PC2", main = paste(title, ": PCA", sep = ""),
      auto.key = list(space = "right"),
      par.settings = list(superpose.symbol = list(pch = rep(1:25))),
      #' par.settings = list(superpose.symbol=list(pch=rep(1:25),
      #'                                           col=c("black","brown3"))),
      #' lwc-15-12-2010: check R_colour_card

      #' scales = "free",
      panel = function(x, y, ...) {
        panel.xyplot(x, y, ...)
        panel.elli.1(x, y, ...)
        #' panel.outl(x,y, ...)  #' wll-29-11-15: need to provide 'labs'
      }, ...
    )

  #' plot the PCA proportion of variance (#' reverse pca.vars for dotplot)
  pca.p.1 <- dotplot(pca.vars[nrow(pca.vars):1, , drop = F],
    groups = T, as.table = T,
    auto.key = list(space = "right"),
    par.settings = list(superpose.symbol = list(pch = rep(1:25))),
    xlab = "Percentage",
    main = paste(title, ": PCA proportion of variance", sep = ""),
    ...
  )
  list(pca.p = pca.p, pca.p.1 = pca.p.1, pca.vars = pca.vars)
}

#' =========================================================================
#' wll-01-06-2015: Wrapper function for plotting MDS. Only the first two
#'   dimensions are  plotted.
mds_plot_wrap <- function(data.list, method = "euclidean",
                          title = "plotting", ...) {
  if (is.null(names(data.list))) {
    names(data.list) <- paste(deparse(substitute(data.list)),
                              1:length(data.list), sep = ":")
  }
  dn <- names(data.list)

  #' MDS
  METHODS <- c(
    "euclidean", "maximum", "manhattan", "canberra",
    "binary", "minkowski"
  )
  meth <- pmatch(method, METHODS)

  mds <- lapply(dn, function(x) {
    dis <- dist(data.list[[x]]$dat, method = METHODS[meth])
    mds <- cmdscale(dis) #' only consider 2 dimension
    mds <- as.data.frame(mds)
    names(mds) <- c("Coord_1", "Coord_2")
    mds <- cbind(mds, cls = data.list[[x]]$cls, type = x)
  })
  names(mds) <- dn
  mds <- do.call(rbind, lapply(mds, function(x) x))

  #' MDS plot
  mds.p <-
    xyplot(Coord_2 ~ Coord_1 | type,
      data = mds, groups = mds$cls, as.table = T,
      xlab = "Coordinate 1", ylab = "Coordinate 2",
      main = paste(title, ": MDS Plot", sep = ""),
      auto.key = list(space = "right"),
      par.settings = list(superpose.symbol = list(pch = rep(1:25))),
      panel = function(x, y, ...) {
        panel.xyplot(x, y, ...)
        panel.elli.1(x, y, ...)
        #' panel.outl(x,y, ...)#' wll-29-11-15: need to provide 'labs'
      }, ...
    )
}

#' =========================================================================
#' wll-23-10-2008: Wrapper function for plotting PCALDA
#' lwc-11-02-2010: replace nlda with pcalda and change correspondingly, e.g.
#'                 DF to LD.
#' lwc-19-10-2010: handle with 2-class and more than 3-class problem. For
#'                 2-class, DF2 is a dummy variable, identical to LD1 for
#'                 general plotting reason.
#' Arguments:
#'   data.list - A two-layer list structure, in which the second layer
#'               include a data frame and a factor of class label. It should
#'               be noted the names of the first layer of data.list must be
#'               given.
#'   title     - A part of title string for plotting
lda_plot_wrap <- function(data.list, title = "plotting", ...) {
  if (is.null(names(data.list))) {
    names(data.list) <- paste(deparse(substitute(data.list)),
                              1:length(data.list), sep = ":")
  }
  dn <- names(data.list)
  lda <- lapply(dn, function(x) { #' x=dn[1]
    res <- pcalda(data.list[[x]]$dat, data.list[[x]]$cls)
    dfs <- as.data.frame(res$x)
    eig <- res$lda.out$svd

    if (ncol(dfs) > 2) {
      dfs <- dfs[, 1:2, drop = F]
      eig <- eig[1:2]
    }
    if (ncol(dfs) == 1) {
      dfs$LD2 <- dfs$LD1
      eig <- c(eig, eig)
    }
    dfs <- cbind(dfs, cls = data.list[[x]]$cls, type = x)

    names(eig) <- c("LD1", "LD2")
    list(dfs = dfs, eig = eig)
  })
  names(lda) <- dn

  lda.dfs <- do.call(rbind, lapply(lda, function(x) x$dfs))
  lda.eig <- do.call(rbind, lapply(lda, function(x) x$eig))

  lda.p <- xyplot(LD2 ~ LD1 | type,
    data = lda.dfs, groups = lda.dfs$cls, as.table = T,
    xlab = "DF1", ylab = "DF2", main = paste(title, ": LDA", sep = ""),
    #' auto.key = list(columns=nlevels(cl)),
    auto.key = list(space = "right"),
    par.settings = list(superpose.symbol = list(pch = rep(1:25))),
    #' scales = "free",
    panel = function(x, y, ...) {
      panel.xyplot(x, y, ...)
      panel.elli.1(x, y, ...)
      #' panel.outl(x,y, ...)#' wll-29-11-15: need to provide 'labs'
    }, ...
  )

  #' plot LDA eigenvales
  lda.p.1 <- dotplot(lda.eig[nrow(lda.eig):1, , drop = F],
    groups = T, as.table = T,
    auto.key = list(space = "right"),
    par.settings = list(superpose.symbol = list(pch = rep(1:25))),
    xlab = "Eigenvalues",
    main = paste(title, ": LDA Eigenvalues", sep = ""), ...
  )
  list(lda.p = lda.p, lda.p.1 = lda.p.1, lda.eig = lda.eig)
}

#' =========================================================================
#' wll-23-10-2008: Wrapper function for plotting PCALDA
#' Note: Will plot 2-class problem differently with stripplot.
lda_plot_wrap.1 <- function(data.list, title = "plotting", ...) {
  if (is.null(names(data.list))) {
    names(data.list) <- paste(deparse(substitute(data.list)),
                              1:length(data.list), sep = ":")
  }
  dn <- names(data.list)
  lda <- lapply(dn, function(x) {
    res <- pcalda(data.list[[x]]$dat, data.list[[x]]$cls)
    dfs <- as.data.frame(res$x)
    dfs <- cbind(dfs,
      cls = data.list[[x]]$cls,
      type = rep(x, nrow(data.list[[x]]$dat))
    )
    #' list(dfs=dfs, eig=res$Tw)
    eig <- res$lda.out$svd
    names(eig) <- colnames(res$x)
    list(dfs = dfs, eig = eig)
  })
  names(lda) <- dn

  lda.dfs <- do.call(rbind, lapply(lda, function(x) x$dfs))
  lda.eig <- do.call(rbind, lapply(lda, function(x) x$eig))

  if (length(grep("LD2", colnames(lda.dfs))) > 0) {
    lda.p <-
      xyplot(LD2 ~ LD1 | type,
        data = lda.dfs, groups = lda.dfs$cls, as.table = T,
        xlab = "LD1", ylab = "LD2", main = paste(title, ": LDA", sep = ""),
        #' auto.key = list(columns=nlevels(cl)),
        auto.key = list(space = "right"),
        par.settings = list(superpose.symbol = list(pch = rep(1:25))),
        #' scales = "free",
        panel = function(x, y, ...) {
          panel.xyplot(x, y, ...)
          panel.elli.1(x, y, ...) #' panel.outl(x,y, ...)
        }, ...
      )
  } else {
    lda.p <-
      stripplot(LD1 ~ cls | type,
        data = lda.dfs, as.table = T, groups = lda.dfs$cls,
        auto.key = list(space = "right"),
        par.settings = list(superpose.symbol = list(pch = rep(1:25))),
        #' scales = "free",
        main = paste(title, ": LDA", sep = ""), ...
      )
    #' wll-07-05-2009: I have added the auto.key and par.settings
    #'   only large number of sub-figures. Otherwise, this two
    #'   line should be removed.
  }

  #' plot LDA eigenvales
  lda.p.1 <- dotplot(lda.eig[nrow(lda.eig):1, , drop = F],
    groups = T, as.table = T,
    auto.key = list(space = "right"),
    par.settings = list(superpose.symbol = list(pch = rep(1:25))),
    xlab = "Eigenvalues",
    main = paste(title, ": LDA Eigenvalues", sep = ""), ...
  )
  list(lda.p = lda.p, lda.p.1 = lda.p.1, lda.eig = lda.eig)
}

#' =========================================================================
#' wll-23-10-2008: Wrapper function for plotting PLSDA. Only the first two
#'                 components are plotted.
#' Note: Use plsc instead of plslda. You can call it PLSDA if PLS is
#'       employed for discrimination.
#' Arguments:
#'   data.list - A two-layer list structure, in which the second layer
#'               include a data frame and a factor of class label. It should
#'               be noted the names of the first layer of data.list must be
#'               given.
#'   title     - A part of title string for plotting
pls_plot_wrap <- function(data.list, title = "plotting", ...) {
  if (is.null(names(data.list))) {
    names(data.list) <- paste(deparse(substitute(data.list)), 1:length(data.list),
      sep = ":"
    )
  }
  dn <- names(data.list)
  pls <-
    lapply(dn, function(x) {
      res <- plsc(data.list[[x]]$dat, data.list[[x]]$cls)
      scores <- as.data.frame(res$x)[, 1:2] #' The first two components
      scores <-
        cbind(scores,
          cls = data.list[[x]]$cls, type = rep(x, nrow(data.list[[x]]$dat))
        )
      vars <- round((res$pls.out$Xvar / res$pls.out$Xtotvar) * 100, 2)[1:2]
      names(vars) <- c("LC1", "LC2")
      list(scores = scores, vars = vars)
    })
  names(pls) <- dn

  pls.scores <- do.call(rbind, lapply(pls, function(x) x$scores))
  pls.vars <- do.call(rbind, lapply(pls, function(x) x$vars))

  pls.p <-
    xyplot(LC2 ~ LC1 | type,
      data = pls.scores, groups = pls.scores$cls, as.table = T,
      xlab = "LC1", ylab = "LC2", main = paste(title, ": PLS", sep = ""),
      auto.key = list(space = "right"),
      par.settings = list(superpose.symbol = list(pch = rep(1:25))),
      #' scales = "free",
      panel = function(x, y, ...) {
        panel.xyplot(x, y, ...)
        panel.elli.1(x, y, ...)
        #' panel.outl(x,y, ...)
      }, ...
    )

  #' plot PLS proportion of variance
  pls.p.1 <- dotplot(pls.vars[nrow(pls.vars):1, , drop = F],
    groups = T, as.table = T,
    auto.key = list(space = "right"),
    par.settings = list(superpose.symbol = list(pch = rep(1:25))),
    xlab = "Percentage",
    main = paste(title, ": PLS proportion of variance", sep = ""), ...
  )
  list(pls.p = pls.p, pls.p.1 = pls.p.1, pls.vars = pls.vars)
}

#' =======================================================================
#' wll-01-10-2009: Misc function for splitting PCA/LDA/PLS plot.
#' Note:  1.) To deal with data.frame and matrix in the same way, should use
#'           a.) colnames instead of names;
#'           b.) index should be tmp[,x] instead of tmp[x] or tmp[[x]].
#'           c.) keep [,,drop=F].
#' TO-DO: 1.) How to extract strip text?
#'        2.) How to keep the legend being consistent with sub-figure's
#'            symbol /colour? It means how to remove the irrelevant
#'            symbol/colours in legend to keep consistent with
#'            sub-figure's.
#' Note: Internal function.
#' Usage:
#'  data(iris)
#'  x <- subset(iris, select = -Species)
#'  y <- iris$Species
#'  #' generate data list by dat.sel
#'  iris.pw <- dat.sel(x,y,choices=NULL)
#'  res <- pls_plot_wrap(iris.pw)
#'  ph  <- plot_wrap.split(res[[1]], res[[3]], perc=F)
#'  win.metafile(filename = paste("pls_plot","%02d.emf",sep="_"))
#'  for(i in 1:length(ph)) plot(ph[[i]])
#'  dev.off()
plot_wrap.split <- function(plot.handle, plot.lab, perc = T) {
  n <- dim(plot.handle)
  pca.ph <- lapply(1:n, function(x) {
    ph <- plot.handle[x]
    xylab <- round(plot.lab[x, , drop = F], digits = 2)
    xylab <- lapply(colnames(xylab), function(y) {
      if (perc) {
        paste(y, " (", xylab[, y], "%)", sep = "")
      } else {
        paste(y, " (", xylab[, y], ")", sep = "")
      }
    })
    ph$ylab <- xylab[[1]]
    if (length(xylab) > 1) ph$xlab <- xylab[[2]]
    ph
  })
}

#' ========================================================================
#' lwc-09-06-2015: plot MDS using lattice package
#' Usage:
#' data(iris)
#' x <- iris[,1:4]
#' y <- iris[,5]
#' mdsplot(x,y, dimen = c(1,2),ep = 2)
#' mdsplot(x,y, dimen = c(2,1),ep = 1)
mdsplot <- function(x, y, method = "euclidean", dimen = c(1, 2), ...) {
  METHODS <- c(
    "euclidean", "maximum", "manhattan", "canberra",
    "binary", "minkowski"
  )
  meth <- pmatch(method, METHODS)

  dis <- dist(x, method = METHODS[meth])

  #' mds <- cmdscale(dis)      #' only consider 2 dimension
  #' mds <- as.data.frame(mds)
  mds <- cmdscale(dis, k = 2, eig = TRUE) #' Classical MDS
  #' mds <- isoMDS(dis, k=2)                  #' Non-metric MDS
  mds <- as.data.frame(mds$points)
  names(mds) <- c("Coordinate 1", "Coordinate 2")

  #' want to change order?
  mds <- mds[, dimen]

  #' call group plot
  p <- grpplot(mds, y, plot = "pairs", ...)
  p
}

#' ========================================================================
#' lwc-29-05-2008: Another version of PCA plot with proportion indication.
#'   Use text indicate different groups.
#' lwc-13-09-2010: major changes and add ellipse plot
pca.plot <- function(x, y, scale = TRUE, abbrev = FALSE, ep.plot = FALSE, ...) {

  #' lwc-12-12-2008: Plot ellipse
  elli.plot <- function(x, y, ...) {
    Var <- var(cbind(x, y))
    Mean <- cbind(mean(x), mean(y))
    Elli <- ellipse(Var, centre = Mean, level = 0.975)
    lines(Elli[, 1], Elli[, 2], ...)
  }

  x <- as.matrix(x)
  y <- factor(y)

  pca <- pca.comp(x, scale = scale, pcs = 1:2, ...)
  val <- pca$scores
  val <- val[c("PC2", "PC1")]

  if (abbrev) levels(y) <- abbreviate(levels(y), abbrev)

  plot(val,
    type = "n", cex = 1.0, cex.lab = 1.0, cex.axis = 1.0, cex.main = 1.0,
    ylab = paste("PC1", " (", pca$vars[1], "%)", sep = ""),
    xlab = paste("PC2", " (", pca$vars[2], "%)", sep = ""), ...
  )

  text(val[, 1], val[, 2], as.character(y), cex = 0.7, col = unclass(y), ...)
  if (ep.plot) {
    tmp <- as.factor(as.numeric(unclass(y)))
    for (i in levels(tmp)) {
      idx <- tmp == i
      elli.plot(val[idx, 1], val[idx, 2], col = i)
    }
    #' Note: I think that it is not a good way to keep the colors
    #'   of ellipse consistent with group text colors.
  }
  invisible(NULL)
}

#' ========================================================================
#' wll-29-03-2008: Compute the PCA scores and proportion of variance
#' TO-DO: 1.) check validity of argument pcs.
pca.comp <- function(x, scale = FALSE, pcs = 1:2, ...) {
  pca <- prcomp(x, scale. = scale, ...)
  vars <- pca$sdev^2 #' i.e. eigenvalues/variance
  vars <- vars / sum(vars) #' Proportion of Variance
  names(vars) <- colnames(pca$rotation)
  vars <- round(vars * 100, 2)
  #' dfn  <- paste(names(vars)," (",vars[names(vars)],"%)",sep="")
  dfn <- paste(names(vars), ": ", vars[names(vars)], "%", sep = "")
  x <- data.frame(pca$x)
  x <- x[, pcs]
  vars <- vars[pcs]
  dfn <- dfn[pcs]

  return(list(scores = x, vars = vars, varsn = dfn))
}

#' ========================================================================
#' wll-17-11-2008: Get number of rejected hypotheses for several multiple
#'                 testing procedures based on Type I error rates.
#' WLL-21-07-2014: Add na.rm = TRUE in sum.
#' Arguments:
#'   adjp  - a matrix-like p-values of simultaneously testing
#'   alpha - a vector of cut-off of p-values or Type I error rate.
#' Note: See mt.reject in package multest.
pval.reject <- function(adjp, alpha) {
  adjp <- as.data.frame(adjp)
  tmp <- sapply(alpha, function(y) {
    p.num <- sapply(adjp, function(x) sum(x <= y, na.rm = TRUE))
  })
  colnames(tmp) <- alpha
  tmp <- t(tmp)
  return(tmp)
}

#' ========================================================================
#' lwc-20-01-2009: Calculate the p-values for columns of data matrix
#'   with respect to group information. Support multiple categorical data.
#' lwc-16-06-2010: Provide argument method. Support user defined test method
#'   which has formula format and returned p.value.
#' Arguments:
#'   x      - data frame or matrix
#'   y      - categorical data
#'   method - hypothesis test such as t.test and wilcox.test.
pval.test <- function(x, y, method = "oneway.test", ...) {
  method <-
    if (is.function(method)) {
      method
    } else if (is.character(method)) {
      get(method)
    } else {
      eval(method)
    }

  pval <- sapply(as.data.frame(x), function(x) {
    method(x ~ y, ...)$p.value

    #' Normality test. Note that H0 is normal distribution!
    #' shapiro.test(x)$p.value
    #' t.test(x ~ y,var.equal=F)$p.value
    #' oneway.test(x ~ y,var.equal=F)$p.value
  })

  #' return(list(pval=pval, method=method))
  return(pval)
}

#' ========================================================================
#' lwc-23-03-2010: corrgram with ellipse
#' lwc-27-04-2012: put scales in argument list
#' Arguments:
#'  co    - Correlation matrices
#'  lable - Logical value indicating whether the correlation coefficient (x
#'          100) should be displayed.
#'  \references{
#'    Michael Friendly (2002).
#'    \emph{Corrgrams: Exploratory displays for correlation matrices}.
#'    The American Statistician, 56, 316--324.
#'    D.J. Murdoch, E.D. Chow (1996).
#'    \emph{A graphical display of large correlation matrices}.
#'    The American Statistician, 50, 178--180.
#'  }
#' Usages:
#'   tmp <- iris[,1:4]
#'   co  <- cor(tmp)
#'   corrgram.ellipse(co,label=T)
#'   corrgram.circle(co)
corrgram.ellipse <- function(co, label = FALSE,
                             col.regions =
                               colorRampPalette(c("red", "white", "blue")),
                             scales = list(x = list(rot = 90)), ...) {
  ph <-
    levelplot(co,
      xlab = NULL, ylab = NULL, at = do.breaks(c(-1.01, 1.01), 20),
      colorkey = list(space = "right"),
      #' col.regions = heat.colors,#terrain.colors,#cm.colors,
      #' col.regions = colorRampPalette(c("yellow", "red")),
      col.regions = col.regions, scales = scales,
      panel = panel.corrgram.ell, label = label, ...
    )
  return(ph)
}

#' ========================================================================
#' lwc-23-03-2010: corrgram with circle/pie
#' lwc-27-04-2012: put scales in argument list
#' Arguments:
#'  co    - Correlation matrices
corrgram.circle <- function(co,
                            col.regions =
                              colorRampPalette(c("red", "white", "blue")),
                            scales = list(x = list(rot = 90)), ...) {
  ph <-
    levelplot(co,
      xlab = NULL, ylab = NULL,
      colorkey = list(space = "right"),
      at = do.breaks(c(-1.01, 1.01), 101),
      col.regions = col.regions, scales = scales,
      panel = panel.corrgram.cir, ...
    )
  return(ph)
}

#' ========================================================================
#' lwc-23-03-2010: Panel function for ellipse corrgram.
#' From Lattice book chapter 13. Internal function.
panel.corrgram.ell <- function(x, y, z, subscripts, at, level = 0.9,
                               label = FALSE, ...) {
  #' require("ellipse", quietly = TRUE)
  x <- as.numeric(x)[subscripts]
  y <- as.numeric(y)[subscripts]
  z <- as.numeric(z)[subscripts]
  zcol <- level.colors(z, at = at, ...)
  for (i in seq(along = z)) {
    ell <- ellipse(z[i],
      level = level, npoints = 50, scale = c(.2, .2),
      centre = c(x[i], y[i])
    )
    panel.polygon(ell, col = zcol[i], border = zcol[i], ...)
  }
  if (label) {
    panel.text(
      x = x, y = y, lab = 100 * round(z, 2),
      cex = 0.8, col = ifelse(z < 0, "white", "black")
    )
  }
}

#' =========================================================================
#' lwc-23-03-2010:Panel function for partially filled circles corrgram.
#' From Lattice book chapter 13. Internal function.
panel.corrgram.cir <- function(x, y, z, subscripts, at = pretty(z),
                               scale = 0.8, ...) {
  x <- as.numeric(x)[subscripts]
  y <- as.numeric(y)[subscripts]
  z <- as.numeric(z)[subscripts]
  zcol <- level.colors(z, at = at, ...)
  for (i in seq(along = z)) {
    lims <- range(0, z[i])
    tval <- 2 * base::pi * seq(from = lims[1], to = lims[2], by = 0.01)
    grid.polygon(
      x = x[i] + .5 * scale * c(0, sin(tval)),
      y = y[i] + .5 * scale * c(0, cos(tval)),
      default.units = "native", gp = gpar(fill = zcol[i])
    )
    grid.circle(
      x = x[i], y = y[i], r = .5 * scale,
      default.units = "native"
    )
  }
}

#' ========================================================================
#' lwc-15-04-2010: pairwise combination of categorical data set
dat.sel <- function(dat, cls, choices = NULL) {
  #' get the index of pairwise combination
  idx <- combn.pw(cls, choices = choices)
  #' construct data set consisting of data matrix and its class info
  dat.pair <-
    lapply(idx, function(x) {
      cls.pw <- factor(cls[x]) #' force drop factor levels

      dat.pw <- dat[x, , drop = F]
      list(dat = dat.pw, cls = cls.pw)
    })
  return(dat.pair)
}

#' ========================================================================
#' lwc-13-04-2010: Index of pairwise combination for categorical vectors.
combn.pw <- function(cls, choices = NULL) {
  .combn.pw <- function(choices, lev) {
    choices <- if (is.null(choices)) lev else unique(choices)
    pm <- pmatch(choices, lev)
    if (any(is.na(pm))) {
      stop("'choices' should be one of ", paste(lev, collapse = ", "))
    }

    #' Get the binary combinations using combn (core package utils)
    if (length(choices) == 1) {
      if (F) { #' simple implementation
        lev.1 <- setdiff(lev, choices)
        com <- cbind(choices, lev.1)
        dimnames(com) <- NULL
      } else { #' Keep comparable with dat.sel.1
        com <- t(combn(lev, 2))
        idx <- sapply(1:nrow(com), function(x) {
          if (match(choices, com[x, ], nomatch = 0) > 0) {
            return(T)
          } else {
            (F)
          }
        })
        com <- com[idx, , drop = F] #' lwc-01-12-2009: fix a bug
      }
    } else {
      com <- t(combn(choices, 2))
    }
    return(com)
  }

  cls <- as.factor(cls)
  lev <- levels(cls)
  if (is.list(choices)) {
    com <- lapply(choices, function(x) .combn.pw(x, lev))
    com <- do.call("rbind", com)
    com <- unique(com)
  } else {
    com <- .combn.pw(choices, lev)
  }
  idx <- apply(com, 1, function(x) {
    ind <- cls %in% x
    #' ind <- which(ind)
    #' comment: don't use this otherwise you have to switch
    #' data frame to list.
  })
  colnames(idx) <- apply(com, 1, paste, collapse = "~")

  idx <- as.data.frame(idx) #' for easy manipulation.
  return(idx)
}

#' ========================================================================
#' lwc-13-08-2006: Generates the pairwise data set based on the class label.
#' History:
#'   18-09-2006: Fix a bug.
#'   31-05-2007: Major changes
#'   01-12-2009: fix a bug when cl is two-class.
#'   21-02-2010: Change name from dat.set to .dat.sel and treat as internal
#'               function
#' NOTE: Using drop=F to keep the format of matrix even the matrix has one
#'       element.
.dat.sel <- function(dat, cl, choices = NULL) {

  #' lwc-29-10-2006: combinations is from package gtools.
  #' $Id: mt_util_1.r,v 1.16 2009/07/27 10:23:41 wll Exp $
  #' From email by Brian D Ripley <ripley@stats.ox.ac.uk> to r-help
  #' dated Tue, 14 Dec 1999 11:14:04 +0000 (GMT) in response to
  #' Alex Ahgarin <datamanagement@email.com>.  Original version was
  #' named "subsets" and was Written by Bill Venables.
  combinations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) {
    if (mode(n) != "numeric" || length(n) != 1
    || n < 1 || (n %% 1) != 0) {
      stop("bad value of n")
    }
    if (mode(r) != "numeric" || length(r) != 1
    || r < 1 || (r %% 1) != 0) {
      stop("bad value of r")
    }
    if (!is.atomic(v) || length(v) < n) {
      stop("v is either non-atomic or too short")
    }
    if ((r > n) & repeats.allowed == FALSE) {
      stop("r > n and repeats.allowed=FALSE")
    }
    if (set) {
      v <- unique(sort(v))
      if (length(v) < n) stop("too few different elements")
    }
    v0 <- vector(mode(v), 0)

    #' Inner workhorse
    if (repeats.allowed) {
      sub <- function(n, r, v) {
        if (r == 0) {
          v0
        } else
        if (r == 1) {
          matrix(v, n, 1)
        } else
        if (n == 1) {
          matrix(v, 1, r)
        } else {
          rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n - 1, r, v[-1]))
        }
      }
    } else {
      sub <- function(n, r, v) {
        if (r == 0) {
          v0
        } else
        if (r == 1) {
          matrix(v, n, 1)
        } else
        if (r == n) {
          matrix(v, 1, n)
        } else {
          rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), Recall(n - 1, r, v[-1]))
        }
      }
    }

    sub(n, r, v[1:n])
  }

  func <- function(choices) {
    if (is.null(choices)) {
      choices <- g
    } else {
      choices <- unique(choices)
    }

    i <- pmatch(choices, g)
    if (any(is.na(i))) {
      stop("'choices' should be one of ", paste(g, collapse = ", "))
    }

    #' Get the binary combinations based on the class labels (package GTOOLS)
    if (length(choices) == 1) {
      com <- combinations(length(g), 2, v = g)
      idx <- sapply(1:nrow(com), function(x) {
        if (match(choices, com[x, ], nomatch = 0) > 0) {
          return(T)
        } else {
          (F)
        }
      })
      com <- com[idx, , drop = F] #' lwc-01-12-2009: fix a bug
    } else {
      com <- combinations(length(choices), 2, v = choices)
    }
    return(com)
  }

  if (missing(dat) || missing(cl)) {
    stop(" The data set and/or class label are missing!")
  }
  cl <- as.factor(cl)
  g <- levels(cl)

  if (is.list(choices)) {
    com <- lapply(choices, function(x) func(x))
    com <- do.call("rbind", com)
    com <- unique(com)
  } else {
    com <- func(choices)
  }

  #' process the data set labels being selected
  dat.sub <- list()
  cl.sub <- list()
  for (i in (1:nrow(com))) {
    idx <- (cl == com[i, ][1]) | (cl == com[i, ][2])
    cl.sub[[i]] <- cl[idx]
    cl.sub[[i]] <- cl.sub[[i]][, drop = T] #' drop the levels
    dat.sub[[i]] <- dat[idx, , drop = F]
  }

  #' get comparison names
  com.names <- apply(com, 1, paste, collapse = "~")
  names(dat.sub) <- names(cl.sub) <- com.names

  return(list(dat = dat.sub, cl = cl.sub, com = com))
}

#' ========================================================================
#' lwc-28-07-2009: panel function for plotting regression line with red
#' color.
#' lwc-29-03-2010: add dots arguments for panel.xyplot.
panel.smooth.line <- function(x, y, ...) {
  panel.grid(h = -1, v = -1)
  panel.xyplot(x, y, ...)
  #' panel.xyplot(x, y, type="p")
  if (sd(y) > 0.001) { # .Machine$double.eps)
    panel.loess(x, y, span = 1, col = "red")
  } else {
    panel.lmline(x, y, col = "red")
  }
}

#' ========================================================================
#' lwc-04-12-2006: Pre-process Data Set
#' lwc-27-03-2007: support multiple methods
#' lwc-27-06-2007: 'rescaler' function in package 'reshape' provides a R
#'   standard to deal with vector, matrix and data.frame using S3 method.
#'   Also another version of range method. Need to check source code to hack
#'   something.
preproc <- function(x, y = NULL, method = "log", add = 1) {

  #' TIC normalisation
  TICnorm <- function(x, y = NULL) {
    scale <- apply(x, 1, function(x) sum(x, na.rm = T))

    if (!is.null(y)) {
      grpm <- as.vector(by(scale, y, mean))
      grpm <- grpm - mean(scale)
      for (k in 1:nlevels(y)) {
        scale[y == levels(y)[k]] <- scale[y == levels(y)[k]] - grpm[k]
      }
    }
    x <- sweep(x, 1, scale, "/")
  }

  me <- function(x) mean(x, na.rm = T)
  se <- function(x) sd(x, na.rm = T)
  mx <- function(x) max(x, na.rm = T)
  mn <- function(x) min(x, na.rm = T)
  #' sm  <- function(x) sum(x,na.rm=T)

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("x must be a matrix or data frame.")
  }
  x <- as.data.frame(x)
  if (!is.null(y)) {
    y <- as.factor(y)
  }

  for (i in method) {
    i <- match.arg(i, c(
      "center", "auto", "range", "pareto", "vast", "level",
      "log", "log10", "sqrt", "asinh", "TICnorm"
    ))

    x <- switch(i,
      #' by colmns
      "center"  = sapply(x, function(x) (x - me(x))),
      "auto"    = sapply(x, function(x) (x - me(x)) / se(x)),
      "range"   = sapply(x, function(x) (x - me(x)) / (mx(x) - mn(x))),
      "pareto"  = sapply(x, function(x) (x - me(x)) / sqrt(se(x))),
      "vast"    = sapply(x, function(x) (x - me(x)) * me(x) / se(x)^2),
      "level"   = sapply(x, function(x) (x - me(x)) / me(x)),
      #' by all
      "log"     = log(x + add),
      "log10"   = log10(x + add),
      "sqrt"    = sqrt(x),
      "asinh"   = asinh(x),
      #' by row
      "TICnorm" = TICnorm(x, y)
    )
  }

  rownames(x) <- 1:nrow(x)
  return(x)
}

#' =========================================================================
#' Remove variables which has (near) zero S.D with/without respect to class.
#' lwc-18-01-2007:  For more details, ?.Machine
#' lwc-15-03-2008:  Fix a bug
#' lwc-01-03-2010:  add na.rm
preproc.sd <- function(x, y = NULL, na.rm = FALSE) {
  if (is.null(y)) {
    #' take off the columns with the same values.
    id <- which(apply(x, 2, sd, na.rm = na.rm) > .Machine$double.eps)
    x <- x[, id]
    return(x)
  } else {
    y <- factor(y)
    #' group s.d. with respect to features
    z <- sapply(data.frame(x), function(i) tapply(i, y, sd, na.rm = na.rm))
    #' minimal s.d.
    z.1 <- sapply(data.frame(z), function(i) min(i))

    #' which one is zero within group?
    if (any(z.1 <= .Machine$double.eps)) {
      z.2 <- which(z.1 <= .Machine$double.eps)
      x <- x[, -z.2, drop = F]
    }
    return(x)
  }
}

#' =========================================================================
#' Remove variables appear to be constant within groups/class
#' lwc-24-01-2007: The function is hacked from lda.default in MASS package
#' lwc-25-01-2007: Fix a bug
preproc.const <- function(x, y, tol = 1.0e-4) {
  if (is.null(dim(x))) stop("'x' is not a matrix or data frame")
  x <- as.matrix(x) #' lwc-04-03-2008: must be matrix
  n <- nrow(x)
  p <- ncol(x)
  if (n != length(y)) {
    stop("nrow(x) and length(y) are different")
  }
  g <- as.factor(y)

  group.means <- tapply(x, list(rep(g, p), col(x)), mean)

  f1 <- sqrt(diag(var(x - group.means[g, ])))

  #' which one is constant within group?
  if (any(f1 < tol)) {
    const <- (1:p)[f1 < tol]
    x <- x[, -const, drop = F]
  }
  x
}

#' ========================================================================
#' lwc-19-06-2008: Correlation analysis of data set and extract the
#'                 pairs with correlation coefficient larger than cutoff
#' lwc-21-10-2010: Minor modify in Manchester: 1.) add abs.f=FALSE; 2). Add
#'                 dot arguments for passing additional info for function
#'                 cor.
#' lwc-23-06-2015: fix a minor bug
#' Arguments:
#'   mat     - A matrix-like data set
#'   cutoff  - A scalar value of threshold
#' Returns:
#'   A data frame with three columns, in which the first and second columns
#'   are variable names and their correlation (lager than cutoff) are
#'   given in the third column.
cor.cut <- function(mat, cutoff = 0.75, abs.f = FALSE,
                    use = "pairwise.complete.obs", method = "pearson", ...) {
  #' co <- cor(mat,use=use, method=method)
  co <- cor(x = mat, use = use, method = method, ...) #' added on 23-06-2015
  co[upper.tri(co)] <- NA
  diag(co) <- NA
  co <- co[-1, -ncol(co), drop = F]
  #' extract items above the cutoff value
  if (abs.f) {
    idx <- which(abs(co) >= cutoff, arr.ind = T)
  } else {
    idx <- which(co >= cutoff, arr.ind = T)
  }

  if (length(idx) != 0) {
    #' tow-columns correlation
    fs1 <- rownames(co)[idx[, 1]]
    fs2 <- colnames(co)[idx[, 2]]
    res <- data.frame(
      com1 = fs1, com2 = fs2, cor = co[idx],
      stringsAsFactors = FALSE
    )
  } else {
    res <- NA
  }
  res
}

#' ========================================================================
#' lwc-16-04-2008: Hierarchical cluster analysis based on correlation
#' analysis.
#' lwc-19-05-2008: Fix a tiny bug
#' lwc-21-05-2008: Check extreme situation
#' lwc-14-10-2009: Change name from my.fs.cor to fs.cor.bas
#' lwc-17-02-2010: change name from fs.cor.bas to cor.hcl.
#' lwc-18-02-2010: add use and method for function cor
#' LWC-13-02-2012: plot cluster use plot method for dendrogram instead of
#' for hclust.
#' Arguments:
#'   mat       - A matrix-like data set
#'   cutoff    - A vector of cutoff (should be in increasing-order)
#'   fig.f     - A logical value for plotting clustering
#' Returns:
#'   A list including all clustering information.
cor.hcl <- function(mat, cutoff = 0.75, use = "pairwise.complete.obs",
                    method = "pearson", fig.f = TRUE, hang = -1,
                    horiz = FALSE, main = "Cluster Dendrogram",
                    ylab = ifelse(!horiz, "1 - correlation", ""),
                    xlab = ifelse(horiz, "1 - correlation", ""), ...) {
  co <- cor(mat, use = use, method = method)
  res <- list()
  if (ncol(co) <= 1) { #' if number of FS is less than 2, no clustering.
    res$all <- co
  } else {
    hc <- hclust(as.dist(1 - co))
    if (fig.f && ncol(co) > 2) { #' 14-10-09: change & to &&
      #' plot(hc, hang=-1,sub="", ylab="1 - correlation", xlab="Variables",
      #'      cex=0.6,...)
      #' lwc-13-02-2012: Not plot hc directly.

      den.hc <- as.dendrogram(hc, hang = hang)
      plot(den.hc, main = main, ylab = ylab, xlab = xlab, horiz = horiz, ...)
      if (horiz) {
        abline(v = 1 - cutoff, col = "red")
      } else {
        abline(h = 1 - cutoff, col = "red")
      }
    }
    id <- cutree(hc, h = 1 - cutoff)
    res <- lapply(unique(id), function(x) {
      cent <- mat[, id == x, drop = FALSE]
      res <- if (ncol(cent) < 2) NA else cor(cent, use = use, method = method)
    })
    #' names(res) <- paste("Cluster",unique(id), sep="_")

    #' shrink the list
    id <- sapply(res, function(x) {
      if (!any(is.na(x))) TRUE else FALSE
    })
    if (all(id == FALSE)) { #' lwc-21-05-2008: Dont't use if (!all(id)) !!!
      res$all <- co
    } else {
      res <- res[id]
      names(res) <- paste("Cluster", 1:length(res), sep = "_")
      res$all <- co
    }
  }
  return(res)
}

#' ========================================================================
#' lwc-18-02-2010: Correlation heatmap using lattice
#' lwc-17-03-2010: add dendrogram.
cor.heat <- function(mat, use = "pairwise.complete.obs", method = "pearson",
                     dend = c("right", "top", "none"), ...) {
  dend <- match.arg(dend)

  co <- cor(mat, use = use, method = method)
  #' Prepare for dendrogram
  dd <- as.dendrogram(hclust(as.dist(1 - co))) #' for correlation only
  ord <- order.dendrogram(dd)

  co.p <-
    switch(dend,
      right =
        levelplot(co[ord, ord],
          aspect = "fill",
          scales = list(x = list(rot = 90), cex = 0.6),
          colorkey = list(space = "bottom"),
          legend = list(right = list(
            fun = dendrogramGrob,
            args = list(
              x = dd, ord = ord,
              side = "right", size = 6
            )
          )), ...
        ),
      top =
        levelplot(co[ord, ord],
          aspect = "fill",
          scales = list(x = list(rot = 90), cex = 0.6),
          colorkey = list(space = "bottom"),
          legend = list(top = list(
            fun = dendrogramGrob,
            args = list(
              x = dd, ord = ord,
              side = "top", size = 6
            )
          )), ...
        ),
      none =
        levelplot(co[ord, ord], #' still want to order them by HCL.
          aspect = "fill",
          scales = list(x = list(rot = 90), cex = 0.6),
          colorkey = list(space = "bottom"), ...
        )
    )
  #' Fix-Me: Is there any efficient and simple way to do switch inside
  #' levelplot?
  return(co.p) #' must use return for lattice object.
}

#' ========================================================================
#' lwc-09-03-2010: Correlation analysis between two data sets
#' lwc-14-09-2010: convert into function from scratch.
#' lwc-05-04-2011: Since the correlation matrix here is not squared, its
#'   order methods are limited. If the similarity matrix is squared, the
#'   functions for ordering objects using hierarchical clustering in package
#'   gclus can be used. These functions are order.single, order.endlink and
#'   order.hclust.
#' Usages
#' x1  <-rnorm(20,40,1)
#' x2  <-rnorm(20,40,2.5)
#' df1 <-data.frame(x1,x2)
#' y1  <-rnorm(20,1,0.47)
#' y2  <-rnorm(20,1,0.59)
#' y3  <-rnorm(20,1,0.75)
#' df2 <-data.frame(y1,y2,y3)
#' cor(df1, df2)
#' cor.heat.gram(df1, df2)
cor.heat.gram <- function(mat.1, mat.2, use = "pairwise.complete.obs",
                          method = "pearson", main = "Heatmap of correlation",
                          cex = 0.75, ...) {
  co <- cor(mat.1, mat.2, use = use, method = method)
  co <- co[complete.cases(co), ]

  if (F) {
    ph <- levelplot(co,
      scales = list(x = list(rot = 90), cex = cex),
      xlab = "", ylab = "", main = main, ...
    )
    #' heatmap.2(co, Rowv=T, Colv=T, col=rev(heat.colors(16)),
    #'          #' distfun=function(x) as.dist(1 - x),
    #'          trace="none", dendrogram = "both", density.info="none")
  }

  #' The heatmap need to be ordered by some rules so we can easily spot some
  #' patterns.
  row.dd <- as.dendrogram(hclust(dist(co)))
  #' not as.dist since co is not squared.
  row.ord <- order.dendrogram(row.dd)
  col.dd <- as.dendrogram(hclust(dist(t(co))))
  col.ord <- order.dendrogram(col.dd)
  ph <-
    levelplot(co[row.ord, col.ord],
      aspect = "fill",
      scales = list(x = list(rot = 60), cex = cex),
      xlab = "", ylab = "", main = main,
      #' main=list(paste("Heatmap of correlation between data - ",des, sep=""),
      #'           cex=cex),
      #' wll-10-09-2015: Use panel.fill() to fill the background with
      #' your 'NA' color.From Deepayan Sarkar
      panel = function(...) {
        panel.fill(col = "black")
        panel.levelplot(...)
      },
      colorkey = list(space = "bottom"),
      #' colorkey = list(space = "bottom", labels=list(cex=cex)),
      legend =
        list(
          right =
            list(
              fun = dendrogramGrob,
              args =
                list(
                  x = col.dd, ord = col.ord,
                  side = "right",
                  size = 10
                )
            ),
          top =
            list(
              fun = dendrogramGrob,
              args =
                list(
                  x = row.dd, ord = row.ord,
                  side = "top",
                  type = "rectangle", size = 5
                )
            )
        ), ...
    )
  ph

  ph.1 <- corrgram.circle(co[row.ord, col.ord], main = main, ...)
  ph.1 <- update(ph.1, scales = list(x = list(rot = 60), cex = cex))
  ph.1

  #' convert short format to long format
  co.1 <- my_melt(co)
  #' co.1 <- reshape::melt(co)
  co.1 <- co.1[complete.cases(co.1), ] #' 17-03-2010: in case NA
  #' co.max  <- max(co.1[,3], na.rm=T)
  #' co.thre <- co.1[co.1[,3] >= 0.4,] #' lwc-09-03-2010: Very specific

  #' lwc-23-06-2015: If co is a symmetric matrix,
  if (F) {
    co.1 <- co
    co.1[upper.tri(co.1)] <- NA
    diag(co.1) <- NA
    co.1 <- co.1[-1, -ncol(co.1), drop = F]

    co.1 <- my_melt(co.1)
    #' co.1 <- reshape::melt(co.1)
    co.1 <- co.1[complete.cases(co.1), ]
  }

  res <- list(cor.heat = ph, cor.gram = ph.1, 
              cor.short = co, cor.long = co.1)

  return(res)
}

#' ========================================================================
#' wl-05-11-2021, Fri: Convert matrix into long format
#' Internal format.  It is used to replace reshape::melt(x). 
my_melt <- function(x) { 
  res <- matrix(x, dimnames = list(t(outer(colnames(x), rownames(x), 
                                           FUN = paste)), NULL))
  res <- as.data.frame(res)
  rn <- rownames(res)
  rn <- do.call("rbind", sapply(rn, strsplit, " "))
  res <- cbind(rn, res)
  dimnames(res) <- list(1:nrow(res), c("X2", "X1", "value"))
  res <- res[c("X1", "X2", "value")]
  return(res)
}

#' =======================================================================
#' lwc-26-04-2008: save a list into a table file.
save.tab <- function(x, filename = "temp.csv", firstline = "\n") {
  # options(warn = -1) #' disable warning message
  write(firstline, file = filename)
  if (is.list(x) && !is.data.frame(x)) { #' lwc-18-02-2010: fix
    for (i in names(x)) {
      write(paste("\n", i, sep = ""), file = filename, sep = ",", append = T)
      write.table(x[[i]],
        file = filename, append = T, sep = ",", na = "",
        quote = F, row.names = T, col.names = NA
      )
    }
  } else {
    write(paste("\n", sep = ""), file = filename, sep = ",", append = T)
    write.table(x,
      file = filename, append = T, sep = ",", na = "",
      quote = F, row.names = T, col.names = NA
    )
  }
  # options(warn = 0) #' restore to default value
  invisible()
}

#' =======================================================================
#' lwc-13-12-2008: Convert a list with components of vector to a data frame
#'  for writing into an Excel file. Shorter vector will be filled with NA.
list2df <- function(x) {
  len <- max(sapply(x, length))
  df <- sapply(x, function(y) c(y, rep(NA, len - length(y))))
  #' lwc-26-06-2008: bug fix. Convert into matrix if fs.order is a vector.
  if (is.vector(df)) df <- t(df)
  return(df)
}

#' ========================================================================
#' lwc-26-04-2008: my version of unlist, which collapse the higher-depths
#' list to 1-depth list. This function uses recursive programming skill to
#' tackle any depths of list.
un.list <- function(x, y = "") {
  res <- list()
  for (i in names(x)) {
    id <- if (y == "") i else paste(y, i, sep = "_")
    if (is.list(x[[i]]) && !is.data.frame(x[[i]])) {
      #' Since data frame has also property of list
      tmp <- un.list(x[[i]], y = id)
      res <- c(res, tmp)
    } else {
      res[[id]] <- x[[i]]
    }
  }
  res
}

#' ========================================================================
#' lwc-16-09-2010: Remove all NULL or NA entries from a list.
#'  Hacked from function compact of package plyr.
#' wll-15-09-2015: Has some problem in new version of R.
shrink.list <- function(x) {
  tmp <- Filter(Negate(is.null), x)
  tmp <- Filter(Negate(is.na), tmp)
  #' Note-16-09-2010: Get a warning if swapping the above two lines.
  return(tmp)
}

#' ========================================================================
#' Generates Class Indicator Matrix from a Factor.
#' A matrix which is zero except for the column corresponding to the class.
#' Internal function.  From package NNET
class.ind <- function(cl) {
  n <- length(cl)
  cl <- as.factor(cl)
  x <- matrix(0, n, length(levels(cl)))
  x[(1:n) + n * (unclass(cl) - 1)] <- 1
  dimnames(x) <- list(names(cl), levels(cl))
  x
}

#'  1) pca.outlier
#'  2) pca.outlier.1
#'  3) grpplot
#'  4) .grpplot
#'  5) pcaplot
#'  6) panel.elli.1
#'  7) panel.elli
#'  8) panel.outl
#'  9) stats.mat
#' 10) stats.vec
#' 11) .pval
#' 12) .foldchange
#' 13) .foldchange2logratio
#' 14) .logratio2foldchange
#' 15) vec.summ.1
#' 16) vec.summ
#' 17) df.summ
#' 18) mv.zene
#' 19) mv.fill
#' 20) mv.stats
#' 21) mv.pattern
#' 22) hm.cols
#' 23) pca_plot_wrap
#' 24) mds_plot_wrap
#' 25) lda_plot_wrap
#' 26) lda_plot_wrap.1
#' 27) pls_plot_wrap
#' 28) plot_wrap.split
#' 29) mdsplot
#' 30) pca.plot
#' 31) pca.comp
#' 32) pval.reject
#' 33) pval.test
#' 34) corrgram.ellipse
#' 35) corrgram.circle
#' 36) panel.corrgram.ell
#' 37) panel.corrgram.cir
#' 38) dat.sel
#' 39) combn.pw
#' 40) .dat.sel
#' 41) panel.smooth.line
#' 42) preproc
#' 43) preproc.sd
#' 44) preproc.const
#' 45) cor.cut
#' 46) cor.hcl
#' 47) cor.heat
#' 48) cor.heat.gram
#' 49) my_melt
#' 50) save.tab
#' 51) list2df
#' 52) un.list
#' 53) shrink.list
#' 54) class.ind

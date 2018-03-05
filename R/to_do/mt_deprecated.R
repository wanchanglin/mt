## =========================================================================
## lwc-15-07-2015: some changes
## Changed functions:
##   1). pca.plot.wrap
##   2). lda.plot.wrap
##   3). lda.plot.wrap.1
##   4). pls.plot.wrap
##   5). mds.plot.wrap
##   6). .grpplot
##   8). grpplot (remove ep)
##   7). pcaplot (remove ep)
## Change RD files:
##   1). grpplot.rd
##   2). pcaplot.rd 
##   3). data.visualisation.rd
##   4). panel.elli.rd
##   5). plot.pcalda.rd
##   6). plot.plsc.rd
##   7). plot.plslda.rd

## =========================================================================
## lwc-05-10-2009: ellipse panel function which support individual and
##  combined group plotting. It is the extension of panel.elli.
## ------------------------------------------------------------------------
## Usage: Under ep=2, there are three options to plot ellipse.
##        com.grp: control which combination of groups to be plotted.
##        no.grp:  control which individual group not to be plotted. Note
##                 it will be overided by com.grp.
##        If no com.grp and no.grp, the each individual group ellipse should
##        be plotted.
##  -----------------------------------------------------------------------
##  load("lda_plot.RData") # it includes lda.dfs and lda.eig
##  p <- xyplot(DF1~DF2|type, data=lda.dfs,groups = y, as.table=T,
##         xlab="DF2", ylab="DF1", main="LDA", auto.key=list(space="right"),
##         par.settings = list(superpose.symbol=list(pch=rep(1:25))),
##         scales = "free", panel=panel.elli.1,
##         #com.grp=list(pre="PR",pos=c("1.5","3","4.5"),voi="0"),
##         #no.grp=c("1.5","3","4.5"),
##         ep=2)
##  p
## ------------------------------------------------------------------------
## Note:It is difficult to manupulate groups and subscripts in
## Lattice package. subscripts is more useful than groups. It seems that
## subscripts have direct relationship with x and y. But groups is not.
## Therefore the groups can be manupulate by the subscripts.
##
panel.elli.1 <- function(x, y, groups, subscripts, ep=0, conf.level = 0.975,
                         com.grp=NULL,no.grp=NULL, ...) {
  ## ------------------------------------------------------------------
  plot.elli <- function(x,y,subscripts,...){ ## plot ellipse
    Var  <- var(cbind(x,y))
    Mean <- cbind(mean(x),mean(y))
    Elli <- ellipse(Var, centre = Mean, level = conf.level)
    ## here conf.level comes from outside

    if (missing(subscripts)){
      panel.points(Elli[,1], Elli[,2],...)
    }
    ## dealing with individual group plotting (called by panel.superpose).
    else {
      grp <- levels(factor(groups[subscripts]))
      ## here groups comes from outside
      if (!(grp %in% no.grp)){
        panel.points(Elli[,1], Elli[,2],...)
      }
    }
    ## label the centre
    ## ltext(x=Mean[1],y=Mean[2], labels="Centre",col="red")
  }
  ## ------------------------------------------------------------------
  panel.superpose(x,y,groups, subscripts,...)
  panel.abline(h=0, v=0,col=c("gray"), lty=2)

  ## overall ellipse based on all data
  if (ep == 1){
    plot.elli(x,y,type="l",col="red",lwd=2,...)
  }
  ## ellipse based on groups, individual or combination.
  else if (ep == 2){
    if (is.null(com.grp)) { ## individual group ellipse
       panel.superpose(x,y,groups=groups, subscripts=subscripts,...,
                       panel.groups = plot.elli, type="l", lty=2)
    } else {                ## combined group ellipse
      grp <- groups[subscripts]
      for (i in names(com.grp)){
        id <- grp %in% com.grp[[i]]
        plot.elli(x[id],y[id],type="l",col="gray",lwd=2,...)
      }
    }
  }
}

## =========================================================================
## lwc-14-02-2010: Panel function for plotting ellipse  used by lattice.
panel.elli <- function(x, y, ep=0,conf.level = 0.975, ...) {
  plot.elli <- function(x,y,...){ ## plot ellipse
    Var  <- var(cbind(x,y))
    Mean <- cbind(mean(x),mean(y))
    Elli <- ellipse(Var, centre = Mean, level = conf.level)
    ## panel.xyplot(x, y,...)
    ## panel.xyplot(Elli[,1], Elli[,2],...)
    panel.points(Elli[,1], Elli[,2],...)
  }

  panel.superpose(x,y,...)
  panel.abline(h=0, v=0,col=c("gray"), lty=2)
  if (ep == 1){
    ## overline ellipse
    plot.elli(x,y,type="l",col="red",lwd=2,...)
  } else if (ep == 2){
    ## group ellipse
    panel.superpose(x,y,..., panel.groups = plot.elli, type="l", lty=2)
  }
}

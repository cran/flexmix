#
#  Copyright (C) 2004-2005 Friedrich Leisch
#  $Id: plot.R 1848 2005-10-06 11:21:12Z leisch $
#

setGeneric("plot")

setMethod("plot", signature(x="flexmix", y="missing"),
function(x, y, mark=NULL, markcol=NULL, col=NULL, 
         eps=1e-4, root=TRUE, ylim=TRUE,
         main=NULL, mfrow=NULL, ...){

    k = length(x@prior)
    opar = par(c("mfrow", "mar", "oma"))
    on.exit(par(opar))

    if(is.null(markcol)) markcol <- FullColors[5]
    if(is.null(col)) col <- LightColors[4]

    if(is.null(main)){
        main = ifelse(root,
                      "Rootogram of posterior probabilities",
                      "Histogram of posterior probabilities")
        main = paste(main, ">", eps)
    }

    if(is.null(mfrow)){
        n = floor(sqrt(k))
        mfrow=c(ceiling(k/n),n)
    }
    
    par(mfrow=mfrow)
    par(mar=c(3,3,1,1))
    par(oma=c(0,0,3,0))

    h = list()
    h1 <- list()
    for(n in 1:k){
        z = x@posterior$scaled[,n]
        ok = z>eps
        h[[n]] = hist(z[ok], breaks=seq(0,1,length=log2(length(z[ok]))+1),
                 plot=FALSE)
        
        if(!is.null(mark)){
            h1[[n]] <- hist(z[ok & (x@cluster==mark)],
                            breaks=h[[n]]$breaks,
                            plot=FALSE)
        }
                        
        
        if(root){
            h[[n]]$counts = sqrt(h[[n]]$counts)
            if(!is.null(mark)){
                h1[[n]]$counts = sqrt(h1[[n]]$counts)
            }
        }
    }

    if(is.logical(ylim)){
        if(ylim)
            ylim = c(0, maxcount = max(sapply(h, function(x) max(x$count))))
        else
            ylim = NULL
    }

    for(n in 1:k){
        plot(h[[n]], xlab="", ylab="", ylim=ylim, main="", col=col, ...)
        title(main=paste("Component", n), line=0)
        if(!is.null(mark)){
            lines(h1[[n]], col=markcol)
        }

    }
    
    title(main=main, outer=TRUE)
    invisible(list(hists=h, marks=h1))
})
        
    
          
###**********************************************************

plotEll <- function(object, data, which=1:2,
                    project=NULL, points=TRUE,
                    eqscale=TRUE, col=NULL,
                    number = TRUE, cex=1.5, numcol="black", 
                    pch=NULL, ...)
{
    if(is.null(col)) col <- rep(FullColors, length=object@k)
    
    if(!is.null(project)){
        data <- predict(project, data)
    }

    type=ifelse(points, "p", "n")

    if(is.null(pch)){
        pch <- (object@cluster %% 10)
        pch[pch==0] <- 10
    }
    else if(length(pch)!=nrow(data)){
        pch <- rep(pch, length=object@k)
        pch <- pch[object@cluster]
    }
    
    if(eqscale)
        MASS::eqscplot(data[,which], col=col[object@cluster],
                       pch=pch, type=type, ...)
    else
        plot(data[,which], col=col[object@cluster], ,
             pch=pch, type=type, ...)
        
    
    for(k in 1:length(object@components)){
        p = parameters(object, k)
        if(!is.null(project)){
            p <- projCentCov(project, p)
        }
        lines(ellipse::ellipse(p$cov[which,which],
                               centre=p$center[which], level=0.5),
              col=col[k], lwd=2)

        lines(ellipse::ellipse(p$cov[which,which],
                               centre=p$center[which], level=0.95),
              col=col[k], lty=2)
    }
    
    ## und nochmal fuer die zentren und nummern (damit die immer oben sind)
    for(k in 1:length(object@components)){
        p = parameters(object, k)
        if(!is.null(project)){
            p <- projCentCov(project, p)
        }
        if(number){
            rad <- ceiling(log10(object@k)) + 1.5
            points(p$center[which[1]],
                   p$center[which[2]],
                   col=col[k], pch=21, cex=rad*cex, lwd=cex,
                   bg="white")
            text(p$center[which[1]],
                 p$center[which[2]], k, cex=cex, col=numcol)
        }
        else{
            points(p$center[which[1]],
                   p$center[which[2]],
                   pch=16, cex=cex, col=col[k])
        }
    }
}

projCentCov <- function(object, p) UseMethod("projCentCov")

projCentCov.default <- function(object, p)
    stop(paste("Cannot handle projection objects of class",
               sQuote(class(object))))

projCentCov.prcomp <- function(object, p)
{
    cent <- matrix(p$center, ncol=length(p$center))
    cent <- scale(cent, object$center, object$scale)  %*% object$rotation

    cov <- p$cov
    if(length(object$scale)>1)
        cov <- outer(object$scale, object$scale, "*") * cov
    cov <- object$rotation %*% cov %*% t(object$rotation)
    
    list(center=cent, cov=cov)
}
    

#
#  FlexMix: Flexible mixture modelling in R
#  Copyright (C) 2004 Friedrich Leisch
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

setGeneric("plot")

setMethod("plot", signature(x="flexmix", y="missing"),
function(x, y, mark=NULL, markcol="red",
         eps=1e-4, root=TRUE, ylim=TRUE,
         main=NULL, mfrow=NULL, ...){

    k = length(x@prior)
    opar = par(c("mfrow", "mar", "oma"))
    on.exit(par(opar))

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
        plot(h[[n]], xlab="", ylab="", ylim=ylim, main="", ...)
        title(main=paste("Component", n), line=0)
        if(!is.null(mark)){
            lines(h1[[n]], col=markcol)
        }

    }
    
    title(main=main, outer=TRUE)
    invisible(list(hists=h, marks=h1))
})
        
    
          
###**********************************************************

plotEll <- function(object, data, ...)
{
    require("MASS")
    require("ellipse")
    
    eqscplot(data, col=object@cluster, ...)
    for(k in 1:length(object@components)){
        p = parameters(object, k)
        lines(ellipse(p$cov, centre=p$center), col=k)
    } 
}


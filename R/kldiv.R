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

setGeneric("KLdiv", function(object, ...) standardGeneric("KLdiv"))

setMethod("KLdiv", "matrix",
function(object, eps=1e-4, ...)
{
    if(!is.numeric(object))
        stop("object must be a numeric matrix\n")
    
    z <- matrix(NA, nrow=ncol(object), ncol=ncol(object))
    colnames(z) <- rownames(z) <- colnames(object)

    for(k in 1:(ncol(object)-1)){
        for(l in 2:ncol(object)){
            ok <- (object[,k] > eps) & (object[,l] > eps)
            if(any(ok)){
                z[k,l] <- sum(object[ok,k] *
                              (log(object[ok,k]) - log(object[ok,l])))
                z[l,k] <- sum(object[ok,l] *
                              (log(object[ok,l]) - log(object[ok,k])))
            }
        }
    }
    z <- sweep(z, 1, colSums(object), "/")
    diag(z)<-0
    z
})


setMethod("KLdiv", "flexmix",
function(object, ...) KLdiv(object@posterior$unscaled, ...))

###**********************************************************



#
#  Copyright (C) 2004-2005 Friedrich Leisch
#  $Id: kldiv.R 1664 2005-06-13 06:11:03Z leisch $
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



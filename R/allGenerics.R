setGeneric("flexmix",
           function(formula, data=list(), k=NULL,
                    cluster=NULL, model=NULL, concomitant=NULL, control=NULL,
                    weights = NULL)
           standardGeneric("flexmix"))

setGeneric("FLXfit",
           function(model, concomitant, control,
                    postunscaled=NULL, groups, weights)
           standardGeneric("FLXfit"))


###**********************************************************

setGeneric("FLXgetModelmatrix",
           function(model, data, ...) standardGeneric("FLXgetModelmatrix"))

setGeneric("FLXfillConcomitant",
           function(concomitant, ...) standardGeneric("FLXfillConcomitant"))

###**********************************************************

setGeneric("logLik")

setGeneric("clogLik", function(object, ...) standardGeneric("clogLik"))

###**********************************************************

setGeneric("FLXcheckComponent", function(model, ...) standardGeneric("FLXcheckComponent"))

setGeneric("FLXgetK", function(model, ...) standardGeneric("FLXgetK"))

setGeneric("FLXgetObs", function(model) standardGeneric("FLXgetObs"))

setGeneric("FLXmstep", function(model, ...) standardGeneric("FLXmstep"))

setGeneric("FLXremoveComponent", function(model, ...) standardGeneric("FLXremoveComponent"))

setGeneric("FLXdeterminePostunscaled", function(model, ...) standardGeneric("FLXdeterminePostunscaled"))

## just to make sure that some S3 generics are available in S4
setGeneric("fitted")
setGeneric("predict")
setGeneric("summary")
setGeneric("unique")


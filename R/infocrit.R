setGeneric("logLik")

setMethod("logLik", signature(object="flexmix"),
function(object, ...){
    z <- object@logLik
    attr(z, "df") <- object@df
    attr(z, "nobs") <- nrow(object@posterior$scaled)
    class(z) <- "logLik"
    z
})

setGeneric("AIC")

setMethod("AIC", signature(object="flexmix"),
function(object, ..., k=2){
    AIC(logLik(object), k=k)
})

if (!isGeneric("BIC")) {
    setGeneric("BIC", function(object, ...)
               standardGeneric("BIC"))
}

setMethod("BIC", "logLik",
          function(object, ...)
          -2 * (c(object) - attr(object, "df") * log(attr(object, "nobs"))/2)
          )

setMethod("BIC", signature(object="flexmix"),
function(object, ...){
    BIC(logLik(object), k=k)
})


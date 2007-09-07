setClass("FLXcontrol",
         representation(iter.max="numeric",
                        minprior="numeric",
                        tolerance="numeric",
                        verbose="numeric",
                        classify="character",
                        nrep="numeric"),
         prototype(iter.max=200,
                   minprior=0.05,
                   tolerance=10e-7,
                   verbose=0,
                   classify="auto",
                   nrep=1))

setAs("list", "FLXcontrol",
function(from, to){
    z = list2object(from, to)
    z@classify = match.arg(z@classify,
                           c("auto", "weighted", "hard", "random", "SEM", "CEM"))
    z
})

setAs("NULL", "FLXcontrol",
function(from, to){
    new(to)
})


###**********************************************************

setClass("FLXM",
         representation(fit="function",
                        defineComponent="expression",
                        weighted="logical",
                        name="character",
                        formula="formula",
                        fullformula="formula",
                        terms="ANY",
                        contrasts="ANY",
                        xlevels="ANY",
                        x="matrix",
                        y="matrix",
                        preproc.x="function",
                        preproc.y="function",
                        "VIRTUAL"),
         prototype(formula=.~.,
                   fullformula=.~.,
                   preproc.x = function(x) x,
                   preproc.y = function(x) x))

## model-based clustering
setClass("FLXMC",
         representation(dist="character"),
         contains = "FLXM")

## regression
setClass("FLXMR",
         contains = "FLXM")

setMethod("show", "FLXM",
function(object){
    cat("FlexMix model of type", object@name,"\n\nformula: ")
    print(object@formula)
    cat("Weighted likelihood possible:", object@weighted,"\n\n")
    if(nrow(object@x)>0){
        cat("Regressors:\n")
        print(summary(object@x))
    }
    if(nrow(object@y)>0){
        cat("Response:\n")
        print(summary(object@y))
    }
    cat("\n")
})

setClass("FLXcomponent",
         representation(df="numeric",
                        logLik="function",
                        parameters="list",
                        predict="function"))

setMethod("show", "FLXcomponent",
function(object){
    if(length(object@parameters)>0)
        print(object@parameters)
})
    


###**********************************************************

setClass("FLXP",
         representation(name="character",
                        formula="formula",
                        x="matrix",
                        fit="function",
                        refit="function",
                        df="function"),
         prototype(formula=~1, df = function(x, k, ...) (k-1)*ncol(x)))

setClass("FLXPmultinom",
         representation(coef="matrix"),
         contains="FLXP")

setMethod("show", "FLXP",
function(object){
    cat("FlexMix concomitant model of type", object@name,"\n\nformula: ")
    print(object@formula)
    if(nrow(object@x)>0){
        cat("\nRegressors:\n")
        print(summary(object@x))
    }
    cat("\n")
})

###**********************************************************

setClass("FLXdist",
         representation(model="list",
                        prior="numeric",
                        components="list",
                        concomitant="FLXP",
                        formula="formula",
                        call="call",
                        k="integer"),
         validity=function(object) {
           (object@k == length(object@prior))
         },
         prototype(formula=.~.))

setClass("flexmix",
         representation(posterior="ANY",
                        weights="ANY",
                        iter="numeric",
                        cluster="integer",
                        logLik="numeric",
                        df="numeric",
                        control="FLXcontrol",
                        group="factor",
                        size="integer",
                        converged="logical",
                        k0="integer"),
         prototype(group=(factor(integer(0))),
                   formula=.~.),
         contains="FLXdist")

setMethod("show", "flexmix",
function(object){
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\nCluster sizes:\n")
    print(object@size)
    cat("\n")
    if(!object@converged) cat("no ")
    cat("convergence after", object@iter, "iterations\n")
})


###**********************************************************

setClass("summary.flexmix",
         representation(call="call",
                        AIC="numeric",
                        BIC="numeric",
                        logLik="logLik",
                        comptab="ANY"))

setMethod("show", "summary.flexmix",
function(object){
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\n")
    print(object@comptab, digits=3)
    cat("\n")
    print(object@logLik)
    cat("AIC:", object@AIC, "  BIC:", object@BIC, "\n")
    cat("\n")    
})

###**********************************************************

setClass("FLXMRglm",
         representation(family="character",
                        refit="function"),
         contains="FLXMR")

setClass("FLXR",
         representation(k="integer",
                        refit="ANY",
                        call="call"))

setMethod("show", signature(object="FLXR"), function(object) {
  cat("\nCall:", deparse(object@call,0.75*getOption("width")),
      sep="\n")
  cat("\nNumber of components:", object@k, "\n\n")
  ## <fixme> for 2.5.0
  if (!length(object@refit)) {
    show(object@refit)
  }
  else if (is(object@refit, "list")) {
    if (any(!sapply(object@refit, function(x) is.null(x@summary)))) show(object@refit)
  }
  else show(object@refit)
})

setClass("FLXRM",
         representation(fitted="ANY",
                        summary="ANY"))

setMethod("show", signature(object="FLXRM"), function(object) {
  if (!is.null(object@summary)) show(object@summary)
})

setClass("FLXRMRglm",
         contains="FLXRM")

setClass("Coefmat",
         extends("matrix"))

setMethod("show", signature(object="Coefmat"), function(object) {
  printCoefmat(object, signif.stars = getOption("show.signif.stars"))
})

###**********************************************************

setClass("FLXnested",
         representation(formula = "list",
                        k = "numeric"),
         validity = function(object) {
           length(object@formula) == length(object@k)
         })

setAs("numeric", "FLXnested",
      function(from, to) {
        new("FLXnested", formula = ~0, k = from)
      })

setAs("list", "FLXnested",
      function(from, to) {
        z = list2object(from, to)
      })

setAs("NULL", "FLXnested",
      function(from, to) {
        new(to)
      })

setMethod("initialize", "FLXnested", function(.Object, formula=~0, k = numeric(0), ...) {
  if (is(formula, "formula")) .Object@formula <- lapply(rep(1,length(k)),
                                                             function(i) formula)
  else .Object@formula <- formula
  .Object@k <- k
  .Object <- callNextMethod()
  .Object
})

###**********************************************************

setClass("FLXMRglmfix",
         representation(design = "matrix",
                        nestedformula = "FLXnested",
                        fixed = "formula",
                        segment = "matrix",
                        variance = "vector"),
         contains="FLXMRglm")

setClass("FLXRMRglmfix",
         representation(design = "matrix"),
         contains = "FLXRMRglm")


###**********************************************************

setClass("FLXRP",
         contains="FLXRM")

setClass("FLXRPmultinom",
         contains="FLXRP")

###**********************************************************

setClassUnion("listOrdata.frame", c("list", "data.frame"))


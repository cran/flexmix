setClass("FLXMRlmer",
         representation(random = "formula", 
                        mc = "call",
                        fr = "list",
                        FL = "list",
                        preproc.z = "function"),
         prototype(preproc.z = function(x, ...) x),
         contains = "FLXMRglm")

defineComponent_lmer <- expression({
    predict <- function(x, ...) x%*%coef
    logLik <- function(x, y, fr, FL, ...) {
      Ztl <- lapply(FL$trms, "[[", "Zt")
      z <- lapply(Ztl, as, "matrix")
      grouping <- FL$fl[[1]]
      llh <- vector(length=nrow(x))
      for (i in seq_len(nlevels(grouping))) {
        index1 <- which(grouping == levels(grouping)[i])
        V <- Reduce("+", c(lapply(seq_along(sigma2$Random), function(u) {
          index2 <- rownames(z[[u]]) %in% levels(grouping)[i]
          t(z[[u]][index2,index1,drop=FALSE]) %*% sigma2$Random[[u]] %*% z[[u]][index2,index1,drop=FALSE]
        }), list(diag(length(index1)) * sigma2$Residual)))
        llh[index1] <- mvtnorm::dmvnorm(y[index1,], mean=predict(x[index1,,drop=FALSE], ...), sigma = V, log=TRUE)/length(index1)
      }
      llh
    }
      
    new("FLXcomponent",
        parameters=list(coef=coef, sigma2=sigma2),
        logLik=logLik, predict=predict,
        df=df)
  })

FLXMRlmer <- function(formula = . ~ ., random, weighted = FALSE,
                      control = list(), eps = .Machine$double.eps)
{
  mc <- match.call()
  mc$weights <- as.symbol("w")
  random <- if (length(random) == 3) random else formula(paste(".", paste(deparse(random), collapse = "")))
  object <- new("FLXMRlmer", formula = formula, random = random,
                family = "gaussian", weighted = weighted, name = "FLXMRlmer:gaussian")
  cv <- do.call(lme4:::lmerControl, control)
  if (weighted) object@preproc.z <- function(FL, x) { 
    if (length(unique(names(FL[["fl"]]))) != 1) stop("only a single variable for random effects is allowed")
    for (i in seq_along(FL[["fl"]])) {
      DIFF <- t(sapply(levels(FL$fl[[i]]), function(id) {
        index1 <- which(FL$fl[[i]] == id)
        index2 <- rownames(FL$trms[[i]]$Zt) == id
        sort(apply(FL$trms[[i]]$Zt[index2, index1, drop=FALSE], 1, paste, collapse = ""))
      }))
      if (length(unique(table(FL[["fl"]][[i]]))) != 1 || nrow(unique(DIFF)) != 1)
        stop("FLXMRlmer does only work correctly if the covariates of the random effects are the same for all observations")
    }
    FL
  }

  lmer.wfit <- function(x, y, w, fr, FL, mc) {
    zero.weights <- any(w < eps)
    if (zero.weights) {
      ok <- w >= eps
      w <- w[ok]
      fr$Y <- fr$Y[ok]
      fr$X <- fr$X[ok,,drop = FALSE]
      if (length(fr$wts) > 0) fr$wts <- fr$wts[ok]
      if (length(fr$off) > 0) fr$off <- fr$off[ok]
      fr$mf <- fr$mf[ok,,drop = FALSE]
      for (i in seq_along(FL$trms)) {
        index <- rownames(FL$trms[[i]]$A) %in% as.character(unique(FL$fl[[i]][ok]))
        FL$fl[[i]] <- factor(FL$fl[[i]][ok])
        FL$trms[[i]]$Zt <- FL$trms[[i]]$Zt[index, ok, drop = FALSE]
        FL$trms[[i]]$A <- FL$trms[[i]]$A[index, ok, drop = FALSE]
      }
      FL$dims["n"] <- sum(ok)
      FL$dims["q"] <- sum(index)
    }
    wts <- sqrt(w)
    fr$Y <- fr$Y * wts
    fr$X <- fr$X * wts
    mer <- lme4:::lmer_finalize(fr, FL, start = NULL, REML = FALSE, verbose = FALSE)
    rm(fr, FL)
    sigma_res <- mer@deviance["sigmaML"]/sqrt(mean(w))
    names(sigma_res) <- NULL
    ans <- mer@ST
    for (i in seq_along(ans)) {
      ai <- ans[[i]]
      dm <- dim(ai)
      ans[[i]] <- if (dm[1] < 2) (ai * sigma_res)^2 else {
        dd <- diag(ai)
        diag(ai) <- rep(1, dm[1])
        sigma_res^2 * crossprod(dd * t(ai))
      }
    }
    list(coefficients = lme4::fixef(mer),
         sigma2 = list(Random  = ans,
           Residual = sigma_res^2),
         df = mer@dims["p"] + 1 + prod(dim(ai) + c(0, 1))/2)
  }
  
  object@defineComponent <- defineComponent_lmer
  
  object@fit <- function(x, y, w, fr, FL, mc){
    fit <- lmer.wfit(x, y, w, fr, FL, mc)
    with(list(coef = coef(fit),
              df = fit$df,
              sigma2 =  fit$sigma2),
         eval(object@defineComponent))
  }
  object
}

setMethod("FLXgetModelmatrix", signature(model="FLXMRlmer"),
          function(model, data, formula, lhs=TRUE, contrasts = NULL, ...)
{
  mc <- match.call()
  formula_nogrouping <- RemoveGrouping(formula)
  if (identical(paste(deparse(formula_nogrouping), collapse = ""), paste(deparse(formula), collapse = ""))) stop("please specify a grouping variable")
  model <- callNextMethod(model, data, formula, lhs)
  random_formula <- update(model@random,
                           paste(".~. |", .FLXgetGroupingVar(formula)))
  fullformula <- model@fullformula
  if (!lhs) fullformula <- fullformula[c(1,3)]
  fullformula <- update(fullformula,
                        paste(ifelse(lhs, ".", ""), "~. + ", paste(deparse(random_formula[[3]]), collapse = "")))
  model@fullformula <- update(model@fullformula,
                              paste(ifelse(lhs, ".", ""), "~. |", .FLXgetGroupingVar(formula)))
  fr <- lme4:::lmerFrames(mc, fullformula, contrasts)
  FL <- lme4:::lmerFactorList(random_formula, fr, 0L, 0L)
  model@fr <- fr
  model@FL <- model@preproc.z(FL, model@x)
  model@mc <- mc
  model
})

setMethod("FLXmstep", signature(model = "FLXMRlmer"),
          function(model, weights, ...)
{
  apply(weights, 2, function(w) model@fit(model@x, model@y, w, model@fr, model@FL, model@mc))
})

setMethod("FLXdeterminePostunscaled", signature(model = "FLXMRlmer"), function(model, components, ...) {
  sapply(components, function(x) x@logLik(model@x, model@y, model@fr, model@FL))
})


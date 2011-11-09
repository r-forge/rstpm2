## package.skeleton(name="rstpm2", path="c:/usr/src/R", force=T, namespace=T, code_files="pm2-3.R")
## Local Windows setup:
## Rtools.bat
## R CMD INSTALL --html "c:/usr/src/R/rstpm2/pkg"
## R CMD build "c:/usr/src/R/rstpm2/pkg"
## R CMD build --binary "c:/usr/src/R/rstpm2/pkg"
##
## Local Ubuntu setup:
## R CMD INSTALL --html ~/src/R/rstpm2/pkg --library=~/R/x86_64-pc-linux-gnu-library/2.12
## R CMD build ~/src/R/rstpm2/pkg
## R CMD build --binary ~/src/R/rstpm2/pkg
##
## testPackage <- TRUE
## if (testPackage) {
##   require(splines)
##   require(survival)
##   require(bbmle)
## }

## extension of ns() to include different boundary derivatives,
## centering and cure
nsx <- 
function (x, df = NULL, knots = NULL, intercept = FALSE,
          Boundary.knots = range(x),
          derivs = if (cure) c(2,1) else c(2,2),
          log=FALSE, 
          centre = FALSE, cure = FALSE, stata.stpm2.compatible=FALSE) 
{
    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) 
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
            Boundary.knots[2L])
    }
    else outside <- FALSE
    if (!missing(df) && missing(knots)) {
        nIknots <- df - 1 - intercept
        if (nIknots < 0) {
            nIknots <- 0
            warning("'df' was too small; have used ", 1 + intercept)
        }
        knots <- if (nIknots > 0) {
          knots <- if (!cure)
            seq.int(0, 1, length.out = nIknots + 2L)[-c(1L, 
                            nIknots + 2L)]
          else c(seq.int(0, 1, length.out = nIknots + 1L)[-c(1L, 
                                 nIknots + 1L)], 0.95)
          if (!stata.stpm2.compatible)
            stats::quantile(x[!outside], knots)
          else stats::quantile(x[!outside], round(knots,2), type=2)
        }
    }
    else nIknots <- length(knots)
    Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
    if (any(outside)) {
        basis <- array(0, c(length(x), nIknots + 4L))
        if (any(ol)) {
            k.pivot <- Boundary.knots[1L]
            xl <- cbind(1, x[ol] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
                1))$design
            basis[ol, ] <- xl %*% tt
        }
        if (any(or)) {
            k.pivot <- Boundary.knots[2L]
            xr <- cbind(1, x[or] - k.pivot)
            tt <- spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
                1))$design
            basis[or, ] <- xr %*% tt
        }
        if (any(inside <- !outside)) 
            basis[inside, ] <- spline.des(Aknots, x[inside], 
                4)$design
    }
    else basis <- spline.des(Aknots, x, 4)$design
    const <- spline.des(Aknots, Boundary.knots, 4, derivs)$design
    if (!intercept) {
        const <- const[, -1, drop = FALSE]
        basis <- basis[, -1, drop = FALSE]
    }
    qr.const <- qr(t(const))
    basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L), 
        drop = FALSE])
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    if (centre) {
      centreBasis <- nsx(centre,
                         knots=if (is.null(knots)) numeric(0) else knots,
                         Boundary.knots=Boundary.knots, 
                         intercept=intercept, derivs=derivs, centre=FALSE, log=log)
      oldAttributes <- attributes(basis)
      basis <- t(apply(basis,1,function(x) x-centreBasis))
      attributes(basis) <- oldAttributes
    }
    a <- list(degree = 3, knots = if (is.null(knots)) numeric(0) else knots, 
        Boundary.knots = Boundary.knots, intercept = intercept, derivs=derivs,
              centre=centre, log=log)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("nsx", "basis", "matrix")
    basis
}
makepredictcall.nsx <- 
function (var, call) 
{
    if (as.character(call)[1L] != "nsx") 
        return(call)
    at <- attributes(var)[c("knots", "Boundary.knots", "intercept",
                            "derivs", "centre", "log")]
    xxx <- call[1L:2]
    xxx[names(at)] <- at
    xxx
}
predict.nsx <- 
function (object, newx, ...) 
{
    if (missing(newx)) 
        return(object)
    a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots", 
        "intercept", "derivs", "centre", "log")])
    do.call("nsx", a)
}
## first derivative of the nsx() function
nsxDeriv<- 
  function (x, dorder=1, ...) {
    stopifnot(dorder %in% 1:2)
    basis <- nsx(x, ...)
    if (dorder==1) {
      h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
      basisD <- apply(predict(basis,x+h)-predict(basis,x-h),2,
                      function(x) x/(2*h))
      if( attr(basis,"log"))
        basisD <- apply(basisD,2,function(y) y/exp(x))
    }
    if (dorder==2) {
      h <- .Machine$double.eps^(1/6)*ifelse(abs(x)>1,abs(x),1)
      basisD <- apply(predict(basis,x+h)-2*predict(basis,x+h)+predict(basis,x-h),2,
                      function(x) x/(h^2))
      if( attr(basis,"log")) {
        h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
        basisD1 <- apply(predict(basis,x+h)-predict(basis,x-h),2,
                      function(x) x/(2*h))
        basisD <- apply(basisD-basisD1,2,function(y) y/exp(2*x))
      }
    }    
    attributes(basisD) <- attributes(basis)
    class(basisD) <- c("nsxDeriv", "basis")
    basisD
  }
## nsxDerivOld <- 
##   function (x, df = NULL, knots = NULL, intercept = FALSE,
##             Boundary.knots = range(x),
##             derivs = if (cure) c(2,1) else c(2,2),
##             centre = FALSE, cure = FALSE) {
##     mf <- match.call(expand.dots = FALSE)
##     mf[[1L]] <- as.name("nsx")
##     basis <- eval(mf, parent.frame())
##     h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
##     basisD <- apply(predict(basis,x+h)-predict(basis,x-h),2,
##                       function(x) x/(2*h))
##     attributes(basisD) <- attributes(basis)
##     class(basisD) <- c("nsxDeriv", "basis")
##     basisD
##   }
makepredictcall.nsxDeriv <- 
  function (var, call) 
{
  if (as.character(call)[1L] != "nsxDeriv") 
    return(call)
  at <- attributes(var)[c("knots", "Boundary.knots", "intercept",
                          "derivs", "centre", "log")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
}
predict.nsxDeriv <- 
  function (object, newx, ...) 
{
  if (missing(newx)) 
    return(object)
  a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots", 
                "intercept", "derivs", "centre", "log")])
  do.call("nsxDeriv", a)
}
##
## nslog <- 
## function (x, df = NULL, knots = NULL, intercept = FALSE,
##           Boundary.knots = range(log(x)),
##           derivs = if (cure) c(2,1) else c(2,2),
##           log=TRUE, # placemarker only (has no effect)
##           centre = FALSE, cure = FALSE, stata.stpm2.compatible=FALSE) 
## {
##     mf <- match.call(expand.dots = FALSE)
##     mf[[1L]] <- as.name("nsx")
##     mf[["x"]] <- call("log",substitute(x))
##     if (!is.null(knots))
##       mf[["knots"]] <- call("log",substitute(knots))
##     if (centre)
##       mf[["centre"]] <- call("log",substitute(centre))
##     mf[["log"]] <- TRUE
##     basis <- eval(mf, parent.frame())
##     class(basis) <- c("nslog", "basis", "matrix")
##     basis
## }
## makepredictcall.nslog <- 
## function (var, call) 
## {
##     if (as.character(call)[1L] != "nslog") 
##         return(call)
##     at <- attributes(var)[c("knots", "Boundary.knots", "intercept",
##                             "derivs", "centre", "log")]
##     xxx <- call[1L:2]
##     xxx[names(at)] <- at
##     xxx
## }
## predict.nslog <- 
## function (object, newx, ...) 
## {
##     if (missing(newx)) 
##         return(object)
##     a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots", 
##         "intercept", "derivs", "centre", "log")])
##     do.call("nslog", a)
## }
## ## first derivative of the nslog() function
## nslogDeriv<- 
##   function (x, ...) {
##     basis <- nslog(x, ...)
##     h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
##     basisD <- apply(predict(basis,x+h)-predict(basis,x-h),2,
##                       function(u) u/(2*h))
##     ## basisD <- apply(basisD,2,function(y) y/x)
##     attributes(basisD) <- attributes(basis)
##     class(basisD) <- c("nslogDeriv", "basis")
##     basisD
##   }
## makepredictcall.nslogDeriv <- 
##   function (var, call) 
## {
##   if (as.character(call)[1L] != "nslogDeriv") 
##     return(call)
##   at <- attributes(var)[c("knots", "Boundary.knots", "intercept",
##                           "derivs", "centre", "log")]
##   xxx <- call[1L:2]
##   xxx[names(at)] <- at
##   xxx
## }
## predict.nslogDeriv <- 
##   function (object, newx, ...) 
## {
##   if (missing(newx)) 
##     return(object)
##   a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots", 
##                 "intercept", "derivs", "centre", "log")])
##   do.call("nslogDeriv", a)
## }
## ## first derivative of bs() function
## bsDeriv<- 
##   function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x)) {
##     mf <- match.call(expand.dots = FALSE)
##     mf[[1L]] <- as.name("bs")
##     basis <- eval(mf, parent.frame())
##     h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
##     ## basis <- ns(x, df=df, knots=knots, intercept=intercept,
##     ##             Boundary.knots=Boundary.knots)
##     basisD <- apply(predict(basis,x+h)-predict(basis,x-h),2,
##                       function(x) x/(2*h))
##     attributes(basisD) <- attributes(basis)
##     class(basisD) <- c("bsDeriv", "basis")
##     basisD
##   }
## makepredictcall.bsDeriv <- 
##   function (var, call) 
## {
##   if (as.character(call)[1L] != "bsDeriv") 
##     return(call)
##   at <- attributes(var)[c("knots", "Boundary.knots", "intercept")]
##   xxx <- call[1L:2]
##   xxx[names(at)] <- at
##   xxx
## }
## predict.bsDeriv <- 
##   function (object, newx, ...) 
## {
##   if (missing(newx)) 
##     return(object)
##   a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots", 
##                 "intercept")])
##   do.call("bsDeriv", a)
## }
## first derivative of the ns() function
## nsDeriv<- 
##   function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x)) {
##     mf <- match.call(expand.dots = FALSE)
##     mf[[1L]] <- as.name("ns")
##     basis <- eval(mf, parent.frame())
##     h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
##     ## basis <- ns(x, df=df, knots=knots, intercept=intercept,
##     ##             Boundary.knots=Boundary.knots)
##     basisD <- apply(predict(basis,x+h)-predict(basis,x-h),2,
##                       function(x) x/(2*h))
##     attributes(basisD) <- attributes(basis)
##     class(basisD) <- c("nsDeriv", "basis")
##     basisD
##   }
## nsDeriv2<- 
##   function (x, ...) {
##     basis <- ns(x, ...)
##     h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
##     basisD <- apply(predict(basis,x+h)-predict(basis,x-h),2,
##                       function(x) x/(2*h))
##     attributes(basisD) <- attributes(basis)
##     class(basisD) <- c("nsDeriv", "basis")
##     basisD
##   }
## nsDeriv3<- 
## function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x)) 
## {
##     nx <- names(x)
##     x <- as.vector(x)
##     nax <- is.na(x)
##     if (nas <- any(nax)) 
##         x <- x[!nax]
##     if (!missing(Boundary.knots)) {
##         Boundary.knots <- sort(Boundary.knots)
##         outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
##             Boundary.knots[2L])
##     }
##     else outside <- FALSE
##     if (!missing(df) && missing(knots)) {
##         nIknots <- df - 1 - intercept
##         if (nIknots < 0) {
##             nIknots <- 0
##             warning("'df' was too small; have used ", 1 + intercept)
##         }
##         knots <- if (nIknots > 0) {
##             knots <- seq.int(0, 1, length.out = nIknots + 2L)[-c(1L, 
##                 nIknots + 2L)]
##             stats::quantile(x[!outside], knots)
##         }
##     }
##     else nIknots <- length(knots)
##     Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
##     if (any(outside)) {
##         basis <- array(0, c(length(x), nIknots + 4L))
##         basisD <- array(0, c(length(x), nIknots + 4L))
##         if (any(ol)) {
##             k.pivot <- Boundary.knots[1L]
##             xl <- cbind(1, x[ol] - k.pivot)
##             xlD <- cbind(0, rep(1,sum(ol)))
##             tt <- splines:::spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
##                 1))$design
##              basis[ol, ] <- xl %*% tt
##             basisD[ol, ] <- xlD %*% tt
##         }
##         if (any(or)) {
##             k.pivot <- Boundary.knots[2L]
##             xr <- cbind(1, x[or] - k.pivot)
##             xrD <- cbind(0, rep(1,sum(or)))
##             tt <- splines:::spline.des(Aknots, rep(k.pivot, 2L), 4, c(0, 
##                 1))$design
##              basis[or, ] <- xr %*% tt
##             basisD[or, ] <- xrD %*% tt
##         }
##         if (any(inside <- !outside)) {
##              basis[inside, ] <- splines:::spline.des(Aknots, x[inside], 
##                  4)$design
##             basisD[inside, ] <- splines:::spline.des(Aknots, x[inside], 
##                 4, rep(1,sum(inside)))$design
##           }
##     }
##     else {
##       basis <- splines:::spline.des(Aknots, x, 4)$design
##       basisD <- splines:::spline.des(Aknots, x, 4, rep(1,length(x)))$design
##     }
##     const <- splines:::spline.des(Aknots, Boundary.knots, 4, c(2, 2))$design
##     if (!intercept) {
##         const <- const[, -1, drop = FALSE]
##         basis <- basis[, -1, drop = FALSE]
##         basisD <- basisD[, -1, drop = FALSE]
##     }
##     qr.const <- qr(t(const))
##     basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L), 
##         drop = FALSE])
##     basisD <- as.matrix((t(qr.qty(qr.const, t(basisD))))[, -(1L:2L), 
##         drop = FALSE])
##     n.col <- ncol(basis)
##     if (nas) {
##         nmat <- matrix(NA, length(nax), n.col)
##         nmat[!nax, ] <- basisD
##         basisD <- nmat
##     }
##     dimnames(basisD) <- list(nx, 1L:n.col)
##     a <- list(degree = 3, knots = if (is.null(knots)) numeric(0) else knots, 
##         Boundary.knots = Boundary.knots, intercept = intercept)
##     attributes(basisD) <- c(attributes(basisD), a)
##     class(basisD) <- c("nsDeriv", "basis")
##     basisD
## }
## makepredictcall.nsDeriv <- 
##   function (var, call) 
## {
##   if (as.character(call)[1L] != "nsDeriv") 
##     return(call)
##   at <- attributes(var)[c("knots", "Boundary.knots", "intercept")]
##   xxx <- call[1L:2]
##   xxx[names(at)] <- at
##   xxx
## }
## predict.nsDeriv <- 
##   function (object, newx, ...) 
## {
##   if (missing(newx)) 
##     return(object)
##   a <- c(list(x = newx), attributes(object)[c("knots", "Boundary.knots", 
##                 "intercept")])
##   do.call("nsDeriv", a)
## }
Shat <- function(obj)
  {
    ## predicted survival for individuals (adjusted for covariates)
    newobj = survfit(obj,se.fit=FALSE)
    surv = newobj$surv
    rr = try(predict(obj,type="risk"),silent=T)
    ## case: only an intercept in the main formula with strata (it would be better to recognise this using attributes for newobj)
    if (inherits(rr,"try-error")) {
      if (try %in% c("Error in colSums(x[j, ] * weights[j]) : \n  'x' must be an array of at least two dimensions\n",
                  "Error in rowsum.default(x * weights, indx) : incorrect length for 'group'\n")) rr <- 1 else stop(rr)
    }
    surv2 = surv[match(obj$y[,ncol(obj$y)-1],newobj$time)]
    return(surv2^rr)
  }
replaceCall=function(obj,old,new) {
  if (is.atomic(obj) && length(obj)>1)
    return(as.call(c(quote(c),lapply(as.list(obj),replaceCall,old,new))))
  if (is.name(obj) || is.symbol(obj) || (is.atomic(obj) && length(obj)==1)) {
    if (obj==old) return(new)
    else return(obj)
  }
##   if (length(obj)==1 && length(obj[[1]])==1) {
##     if (obj==old) return(new)
##     else return(obj)
##   }
  as.call(lapply(obj,replaceCall,old,new))
}
replaceFormula=function(...) as.formula(replaceCall(...))
## replaceFormula(~f(a+b),quote(f),quote(g))
allCall=function(obj) {
  if (is.atomic(obj) && length(obj)==1) return(obj)
  if (is.atomic(obj) && length(obj)>1) return(as.call(c(quote(c),as.list(obj))))
  if (is.name(obj) || is.symbol(obj)) return(obj)
  as.call(lapply(obj,allCall))
}
## allCall(as.call(c(quote(ns),list(df=3,knots=c(1,2)))))[[2]]
vector2call=function(obj) {
  if (is.atomic(obj) && length(obj)==1) return(obj)
  if (is.atomic(obj) && length(obj)>1) return(as.call(c(quote(c),as.list(obj))))
  if (is.name(obj) || is.symbol(obj)) return(obj)
  lapply(obj,allCall) # is this correct?
}
## vector2call(list(df=3,knots=c(1,2)))
rhs=function(formula) 
  if (length(formula)==3) formula[[3]] else formula[[2]]
lhs <- function(formula) 
  if (length(formula)==3) formula[[2]] else NULL
"rhs<-" = function(formula,value) {
  newformula <- formula
  newformula[[length(formula)]] <- value
  newformula
}
"lhs<-" <- function(formula,value) {
  if (length(formula)==2)
    as.formula(as.call(c(formula[[1]],value,formula[[2]])))
  else {
    newformula <- formula
    newformula[[2]] <- value
    newformula
  }
}
## numerically calculate the gradient (func may return a vector)
grad <- function(func,x,...) # would shadow numDeriv::grad()
  {
    h <- .Machine$double.eps^(1/3)*ifelse(abs(x)>1,abs(x),1)
    temp <- x+h
    h.hi <- temp-x
    temp <- x-h
    h.lo <- x-temp
    twoeps <- h.hi+h.lo
    nx <- length(x)
    ny <- length(func(x,...))
    if (ny==0L) stop("Length of function equals 0")
    df <- if(ny==1L) rep(NA, nx) else matrix(NA, nrow=nx,ncol=ny)
    for (i in 1L:nx) {
      hi <- lo <- x
      hi[i] <- x[i] + h.hi[i]
      lo[i] <- x[i] - h.lo[i]
      if (ny==1L)
        df[i] <- (func(hi, ...) - func(lo, ...))/twoeps[i]
      else df[i,] <- (func(hi, ...) - func(lo, ...))/twoeps[i]
      }
    return(df)
  }
## fun: takes coef as its first argument
## requires: coef() and vcov() on the object
numDeltaMethod <- function(object,fun,...) {
  coef <- coef(object)
  est <- fun(coef,...)
  Sigma <- vcov(object)
  gd <- grad(fun,coef,...)
  se.est <- as.vector(sqrt(diag(t(gd) %*% Sigma %*% gd)))
  data.frame(Estimate = est, SE = se.est)
}
predictnl <- function (object, ...) 
  UseMethod("predictnl")
predictnl.default <- function(object,fun,newdata=NULL,...)
  { ## link=c(I,log,sqrt),invlink=NULL
##     link <- match.arg(link)
##     if (is.null(invlink))
##       invlink <- switch(deparse(substitute(link)),I=I,log=exp,sqrt=function(x) x^2)
    if (is.null(newdata) && !is.null(object$data))
      newdata <- object$data
    localf <- function(coef,...)
      {
        object$coefficients = coef
        fun(object,...)
      }
    numDeltaMethod(object,localf,newdata=newdata,...)
  }
setMethod("predictnl", "mle2", function(object,fun,newdata=NULL,...)
  {
    if (is.null(newdata) && !is.null(object@data))
      newdata <- object@data
    localf <- function(coef,...)
      {
        object@fullcoef = coef # changed from predictnl.default()
        fun(object,...)
      }
    numDeltaMethod(object,localf,newdata=newdata,...)
  })
## setMethod("predictnl", "mle", function(object,fun,...)
##   {
##     localf <- function(coef,...)
##       {
##         object@fullcoef = coef # changed from predictnl.default()
##         fun(object,...)
##       }
##     numDeltaMethod(object,localf,...)
##   })
predict.formula <- function(formula,data,newdata,na.action,type="model.matrix") 
{
  mf <- match.call(expand.dots = FALSE)
  type <- match.arg(type)
  m <- match(c("formula", "data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  xlevels <-.getXlevels(mt, mf)
  mfnew <- model.frame(mt, newdata, na.action=na.action, xlev=xlevels)
  if (!is.null(cl <- attr(mt, "dataClasses"))) .checkMFClasses(cl, mfnew)
  model.matrix(mt, mfnew, contrasts=contrasts)
}
`%call+%` <- function(left,right) call("+",left,right)
##
bread.stpm2 <- function (x, ...) {
  rval <- vcov(x) * nrow(x@y)
  dimnames(rval) <- list(names(coef(x)), names(coef(x)))
  return(rval)
}
estfun.stpm2 <- function(obj, weighted=FALSE, ...) {
  rr <- t(rstpm2:::grad(obj@logli,coef(obj)))
  colnames(rr) <- names(coef(obj))
  if (weighted)
    rr <- rr * obj@weights
  rr
}
meat.stpm2 <- 
function (x, adjust = FALSE, ...) 
{
    psi <- estfun.stpm2(x, ...)
    k <- NCOL(psi)
    n <- NROW(psi)
    rval <- crossprod(as.matrix(psi))/n
    if (adjust) 
        rval <- n/(n - k) * rval
    rownames(rval) <- colnames(rval) <- colnames(psi)
    return(rval)
}
sandwich.stpm2 <- 
function (x, bread. = bread.stpm2, meat. = meat.stpm2, ...) 
{
    if (is.function(bread.)) 
        bread. <- bread.(x)
    if (is.function(meat.)) 
        meat. <- meat.(x, ...)
    n <- NROW(estfun.stpm2(x))
    return(1/n * (bread. %*% meat. %*% bread.))
}
setOldClass("terms")
setClassUnion("listOrNULL",c("list","NULL"))
##setClassUnion("numericOrNULL",c("numeric","NULL"))
setOldClass("Surv")
setClass("stpm2", representation(xlevels="list",
                                 contrasts="listOrNULL",
                                 terms="terms",
                                 logli="function",
                                 ## weights="numericOrNULL",
                                 model.frame="list",
                                 x="matrix",
                                 xd="matrix",
                                 termsd="terms",
                                 Call="call",
                                 y="Surv"
                                 ),
         contains="mle2")

stpm2 <- function(formula, data,
                  df=3, logH.args=NULL, logH.formula=NULL,
                  tvc=NULL, tvc.formula=NULL,
                  control=list(parscale=0.1,maxit=300),
                  coxph.strata=NULL, weights=NULL, robust=FALSE,
                  bhazard=NULL, contrasts=NULL, subset=NULL, ...)
  {
    Call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "contrasts", "weights"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    ## mf$drop.unused.levels <- TRUE # include?
    mf[[1L]] <- as.name("model.frame")
    eventExpression <- lhs(formula)[[length(lhs(formula))]]
    delayed <- length(lhs(formula))==4
    timevar <- lhs(formula)[[if (delayed) 3 else 2]]
    ## set up the formulae
    if (is.null(logH.formula) && is.null(logH.args))
      logH.args$df <- df
    if (!is.null(logH.args) && is.null(logH.args$log))
      logH.args$log <- TRUE
    if (is.null(logH.formula))
      logH.formula <- as.formula(call("~",as.call(c(quote(nsx),call("log",timevar),
                                                    vector2call(logH.args)))))
    full.formula <- formula
    rhs(full.formula) <- rhs(formula) %call+% rhs(logH.formula)
    if (is.null(tvc.formula) && !is.null(tvc)) {
      tvc.formulas <-
        lapply(names(tvc), function(name)
               call(":",
                    as.name(name),
                    as.call(c(quote(nsx),
                              call("log",timevar),
                              vector2call(list(log=TRUE,df=tvc[[name]]))))))
      if (length(tvc.formulas)>1)
        tvc.formulas <- list(Reduce(`%call+%`, tvc.formulas))
      tvc.formula <- as.formula(call("~",tvc.formulas[[1]]))
    }
    if (!is.null(tvc.formula))
      rhs(full.formula) <- rhs(full.formula) %call+% rhs(tvc.formula)
    logHD.formula <- replaceFormula(logH.formula,quote(nsx),quote(nsxDeriv))
    if (!is.null(tvc.formula)) {
      tvcD.formula <- replaceFormula(tvc.formula,quote(nsx),quote(nsxDeriv))
      rhs(logHD.formula) <- logHD.formula[[2]] %call+% tvcD.formula[[2]]
    }
    ## set up primary terms objects (mt and mtd)
    mf$formula = full.formula
    ## mf$subset <- eventExpression # call("&",call("(",eventExpression),call("(",substitute(subset)))
    datae <- eval(call("subset",substitute(data),eventExpression),parent.frame()) # data required?
    mf$data <- quote(datae) # restricted to event times
    mfX <- mfd <- mf # copy
    ## mf <- eval(mf, parent.frame())
    mf <- eval(mf)
    mt <- attr(mf, "terms") # primary!
    xlev <- .getXlevels(mt, mf)
    mfd[[2]] <- logHD.formula
    ## mfd <- eval(mfd, parent.frame())
    mfd <- eval(mfd)
    mtd <- attr(mfd, "terms") # primary!
    ## design matrices
    ## mfX <- model.frame(mt, data, xlev = xlev, weights = weights)
    mfX$formula <- quote(mt)
    mfX$data <- quote(data)
    mfX2 <- mfXD <- mfX # copies
    mfX <- eval(mfX)
    if (!is.null(cl <- attr(mt, "dataClasses"))) 
      .checkMFClasses(cl, mfX)
    X <- model.matrix(mt, mfX, contrasts)
    wt <- model.weights(mfX)
    if (is.null(wt)) wt <- rep(1,nrow(X))
    ## mfXD <- model.frame(mtd, data, xlev=xlev, weights = weights)
    mfXD$formula <- quote(mtd)
    mfXD <- eval(mfXD)
    if (!is.null(cl <- attr(mtd, "dataClasses"))) 
      .checkMFClasses(cl, mfXD)
    XD <- model.matrix(mtd, mfXD, contrasts)[,-1,drop=FALSE]
    ##
    y <- model.extract(mfX,"response")
    if (!inherits(y, "Surv")) 
      stop("Response must be a survival object")
    type <- attr(y, "type")
    if (type != "right" && type != "counting") 
        stop(paste("stpm2 model doesn't support \"", type, "\" survival data", 
            sep = ""))
    event <- y[,ncol(y)]==1
    time <- y[,ncol(y)-1]
    ## initial values
    coxph.strata <- substitute(coxph.strata)
    coxph.formula <- formula
    if (!is.null(coxph.strata)) 
      rhs(coxph.formula) <- rhs(formula) %call+% call("strata",coxph.strata)
    coxph.obj <- coxph(coxph.formula,data=data,model=TRUE)
    ## coxph.obj <- eval.parent(substitute(coxph(formula,data),
    ##                                     list(formula=formula,data=data)))
    data$logHhat <- pmax(-18,log(-log(Shat(coxph.obj))))
    ##
    lm.formula <- full.formula
    lhs(lm.formula) <- quote(logHhat) # new response
    init <- coef(lm(lm.formula,data[event,],contrasts=contrasts))
    indexXD <- (length(coef(coxph.obj))+2):ncol(X)
    bhazard <- if (is.null(bhazard)) 0 else bhazard[event] # crude
    if (delayed && any(y[,1]>0)) {
      data2 <- data[y[,1]>0,,drop=FALSE] # data for delayed entry times
      mt2 <- delete.response(mt)
      ## hack: copy over times - so we can use the same term object
      y.names <- sapply(lhs(full.formula),deparse)
      data2[[y.names[3]]] <- data2[[y.names[2]]] 
      ## data2[[y.names[2]]] <- 0
      mfX2$formula <- quote(mt2)
      mfX2$data <- quote(data2)
      mfX2 <- eval(mfX2)
      ##mfX2 <- model.frame(mt2, data2, xlev=xlev, weights = weights)
      if (!is.null(cl <- attr(mt2, "dataClasses"))) 
        .checkMFClasses(cl, mfX2)
      X2 <- model.matrix(mt2, mfX2, contrasts)
      ## delayed.formula <- full.formula
      ## lhs(delayed.formula) <- NULL
      ## X2 = predict(delayed.formula, data, data2) ## delayed entry
      negll <- function(beta) {
        eta <- X %*% beta
        eta2 <- X2 %*% beta
        h <- (XD[event,] %*% beta[indexXD])*exp(eta[event]) + bhazard
        ## h <- (XD[event,] %*% beta[indexXD])*exp(eta[event])/time[event] + bhazard
        h[h<0] <- 1e-100
        ll <- sum(wt[event]*log(h)) +  sum(wt*exp(eta2)) -
          sum(wt*exp(eta))
        return(-ll)
      }
      logli <- function(beta) {
        eta <- X %*% beta
        eta2 <- X2 %*% beta
        h <- (XD[event,] %*% beta[indexXD])*exp(eta[event]) + bhazard
        h[h<0] <- 1e-100
        out <- exp(eta2) - exp(eta)
        out[event] <- out[event]+log(h)
        out <- out*wt
        return(out)
      }
    }
    else { # right censoring only
      negll <- function(beta) {
        eta <- X %*% beta
        h <- (XD[event,] %*% beta[indexXD])*exp(eta[event]) + bhazard
        ## h <- (XD[event,] %*% beta[indexXD])*exp(eta[event])/time[event] + bhazard
        h[h<0] <- 1e-100
        ll <- sum(wt[event]*log(h)) - sum(wt*exp(eta))
        return(-ll)
      }
      logli <- function(beta) {
        eta <- X %*% beta
        h <- (XD[event,] %*% beta[indexXD])*exp(eta[event]) + bhazard
        h[h<0] <- 1e-100
        out <- -exp(eta)
        out[event] <- out[event]+log(h)
        out <- out*wt
        return(out)
      }
    }
    ## MLE
    if (!is.null(control) && "parscale" %in% names(control)) {
      if (length(control$parscale)==1)
        control$parscale <- rep(control$parscale,length(init))
      if (is.null(names(control$parscale)))
        names(control$parscale) <- names(init)
    }
    parnames(negll) <- names(init)
    mle2 <- mle2(negll,init,vecpar=TRUE, control=control, ...)
    out <- new("stpm2",
               call = mle2@call,
               call.orig = mle2@call,
               coef = mle2@coef,
               fullcoef = mle2@fullcoef,
               vcov = mle2@vcov,
               min = mle2@min,
               details = mle2@details,
               minuslogl = mle2@minuslogl,
               method = mle2@method,
               data = data,
               formula = mle2@formula,
               optimizer = "optim",
               xlevels = .getXlevels(mt, mf),
               contrasts = attr(X, "contrasts"),
               logli = logli,
               ##weights = weights,
               Call = Call,
               terms = mt,
               model.frame = mf,
               x = X,
               xd = XD,
               termsd = mtd,
               y = y)
    if (robust) # kludge
      out@vcov <- sandwich.stpm2(out)
    return(out)
  }

##
## predict survival
cloglog <- function(x) log(-log(x))
cexpexp <- function(x) exp(-exp(x))
setMethod("predictnl", "stpm2",
          function(object,fun,newdata=NULL,link=c("I","log","cloglog"),...)
  {
    link <- match.arg(link)
    invlinkf <- switch(link,I=I,log=exp,cloglog=cexpexp)
    linkf <- eval(parse(text=link))
    if (is.null(newdata) && !is.null(object@data))
      newdata <- object@data
    localf <- function(coef,...)
      {
        object@fullcoef = coef
        linkf(fun(object,...))
      }
    dm <- numDeltaMethod(object,localf,newdata=newdata,...)
    out <- invlinkf(data.frame(Estimate=dm$Estimate,
                               lower=dm$Estimate-1.96*dm$SE,
                               upper=dm$Estimate+1.96*dm$SE))
    ## cloglog switches the bounds
    if (link=="cloglog") 
      out <- data.frame(Estimate=out$Estimate,lower=out$upper,upper=out$lower)
    return(out)
  })
##
setMethod("predict", "stpm2",
          function(object,newdata=NULL,type=c("surv","cumhaz","hazard","hr","sdiff","hdiff","loghazard","link"),
                   grid=FALSE,seqLength=300,
                   se.fit=FALSE,link=NULL,exposed=incrVar(var),var,...)
  {
    ## exposed is a function that takes newdata and returns the revised newdata
    ## var is a string for a variable that defines a unit change in exposure
    local <-  function (object, newdata=NULL, type="surv", exposed)
      {
        tt <- object@terms 
        if (is.null(newdata)) {
          ##mm <- X <- model.matrix(object) # fails (missing timevar)
          X <- object@x
          XD <- object@xd
          y <- model.response(object@model.frame) 
          time <- as.vector(y[,ncol(y)-1])
        }
        else {
          Terms <- delete.response(tt)
          m <- model.frame(Terms, newdata, xlev = object@xlevels)
          if (!is.null(cl <- attr(Terms, "dataClasses"))) 
            .checkMFClasses(cl, m)
          X <- model.matrix(Terms, m, contrasts.arg = object@contrasts)
          resp <- attr(Terms, "variables")[attr(Terms, "response")]
          ## similarly for the derivatives
          if (type %in% c("hazard","hr","sdiff","hdiff","loghazard")) {
            ttd <- object@termsd
            TermsD <- delete.response(ttd)
            md <- model.frame(TermsD, newdata, xlev = object@xlevels)
            if (!is.null(cld <- attr(TermsD, "dataClasses"))) 
              .checkMFClasses(cld, md)
            XD <- model.matrix(TermsD, md, contrasts.arg = object@contrasts)[,-1,drop=FALSE]
            ## how to elegantly extract the time variable?
            timevar <- if (length(tt[[2]])==3) tt[[2]][[2]] else tt[[2]][[3]]
            time <- model.matrix(as.formula(call("~",timevar)),newdata)[,-1,drop=TRUE]
            ##
          }
          if (type %in% c("hr","sdiff","hdiff")) {
            if (missing(exposed))
              stop("exposed needs to be specified for type in ('hr','sdiff','hdiff')")
            newdata2 <- exposed(newdata)
            m2 <- model.frame(Terms, newdata2, xlev = object@xlevels)
            if (!is.null(cl <- attr(Terms, "dataClasses"))) 
              .checkMFClasses(cl, m2)
            X2 <- model.matrix(Terms, m2, contrasts.arg = object@contrasts)
            md2 <- model.frame(TermsD, newdata2, xlev = object@xlevels)
            if (!is.null(cld <- attr(TermsD, "dataClasses"))) 
              .checkMFClasses(cld, md2)
            XD2 <- model.matrix(TermsD, md2, contrasts.arg = object@contrasts)[,-1,drop=FALSE]
          }
        }
        beta <- coef(object)
        cumHaz = exp(X %*% beta)
        Sigma = vcov(object)
        if (type=="link") { # delayed entry?
          return(X %*% beta)
        }
        if (type=="cumhaz") { # delayed entry?
          return(cumHaz)
        }
        if (type=="surv") { # delayed entry?
          return(exp(-cumHaz))
        }
        if (type=="sdiff")
          return(exp(-exp(X2 %*% beta)) - exp(-cumHaz))
        if (type=="hazard") {
          betaXD <- beta[(ncol(X)-ncol(XD)+1):ncol(X)]
          ##return((XD %*% betaXD)/time*cumHaz)
          return((XD %*% betaXD)*cumHaz)
        }
        if (type=="loghazard") {
          betaXD <- beta[(ncol(X)-ncol(XD)+1):ncol(X)]
          ##return((XD %*% betaXD)/time*cumHaz)
          return(log(XD %*% betaXD)+log(cumHaz))
        }
        if (type=="hdiff") {
          betaXD <- beta[(ncol(X)-ncol(XD)+1):ncol(X)]
          ##return((XD2 %*% betaXD)/time*exp(X2 %*% beta) - (XD %*% betaXD)/time*cumHaz)
          return((XD2 %*% betaXD)*exp(X2 %*% beta) - (XD %*% betaXD)/time*cumHaz)
        }
        if (type=="hr") {
          cumHazRatio = exp((X2 - X) %*% beta)
          betaXD <- beta[(ncol(X)-ncol(XD)+1):ncol(X)]
          return((XD2 %*% betaXD)/(XD %*% betaXD)*cumHazRatio)
        }
      }
    ##debug(local)
    type <- match.arg(type)
    if (is.null(newdata) && type %in% c("hr","sdiff","hdiff"))
      stop("Prediction using type in ('hr','sdiff','hdiff') requires newdata to be specified.")
    if (grid) {
      Terms <- object@terms
      timevar <- lhs(Terms)[[length(lhs(Terms))-1]]
      Y <- object@y
      event <- Y[,ncol(Y)]==1
      time <- Y[,ncol(Y)-1]
      eventTimes <- time[event]
      X <- seq(min(eventTimes),max(eventTimes),length=seqLength)[-1]
      data.x <- data.frame(X)
      names(data.x) <- deparse(timevar)
      newdata <- merge(newdata,data.x)
    }
    pred <- if (!se.fit) {
      local(object,newdata,type=type,exposed=exposed,
            ...)
    }
    else {
      if (is.null(link))
        link <- switch(type,surv="cloglog",cumhaz="log",hazard="log",hr="log",sdiff="I",
                       hdiff="I",loghazard="I",link="I")
      predictnl(object,local,link=link,newdata=newdata,type=type,
                exposed=exposed,...) 
    }
    attr(pred,"newdata") <- newdata
    ##if (grid) cbind(newdata,as.data.frame(pred)) else pred
    return(pred)
  })
##`%c%` <- function(f,g) function(...) g(f(...)) # function composition
## to do:
## (*) Stata-compatible knots
incrVar <- function(var,increment=1) {
  ##var <- deparse(substitute(var))
  ##function(data) "$<-"(data,var,"$"(data,var)+increment) # FAILS
  n <- length(var)
  if (n>1)
    increment <- rep(increment,n)
  function(data) {
    for (i in 1:n) {
      data[[var[i]]] <- data[[var[i]]] + increment[i]
    }
    data
  }
}
setMethod("plot", signature(x="stpm2", y="missing"),
          function(x,y,newdata,type="surv",
                      xlab="Time",line.col=1,ci.col="grey",
                      add=FALSE,ci=TRUE,rug=TRUE,
                      var=NULL,...) {
  y <- predict(x,newdata,type=type,var=var,grid=T,se.fit=T)
  ylab <- switch(type,hr="Hazard ratio",hazard="Hazard",surv="Survival",
                 sdiff="Survival difference",hdiff="Hazard difference")
  xx <- attr(y,"newdata")
  xx <- xx[,ncol(xx)]
  if (!add) matplot(xx, y, type="n", xlab=xlab, ylab=ylab, ...)
  if (ci) polygon(c(xx,rev(xx)), c(y[,2],rev(y[,3])), col=ci.col, border=ci.col)
  lines(xx,y[,1],col=line.col)
  if (rug) {
      Y <- x@y
      eventTimes <- Y[Y[,ncol(Y)]==1,ncol(Y)-1]
      rug(eventTimes,col=line.col)
    }
  return(invisible(y))
})

if (F) { ##### examples #####
require(foreign)
if (F) { # testing in open code
  install.packages("bbmle", repos="http://R-Forge.R-project.org")
  require(bbmle)
  brcancer=read.dta("brcancer.dta")
  brcancer=transform(brcancer,rate0=10^(-5+x1/100))
}
try(detach("package:bbmle",unload=TRUE),silent=TRUE)

try(detach("package:rstpm2",unload=TRUE),silent=TRUE)
require(rstpm2)
data(brcancer)
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(rectime),df=3,stata=TRUE)))

brcancer <- transform(brcancer,w=10)
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                     weights=w,robust=T,
                     logH.formula=~nsx(log(rectime),df=3,stata=TRUE)))


## sandwich variance estimator (from the sandwich package)

coeftest.stpm2 <- 
function (x, vcov. = NULL, df = NULL, ...) 
{
    est <- coef(x)
    if (is.null(vcov.)) 
        se <- vcov(x)
    else {
        if (is.function(vcov.)) 
            se <- vcov.(x)
        else se <- vcov.
    }
    se <- sqrt(diag(se))
    if (!is.null(names(est)) && !is.null(names(se))) {
        anames <- names(est)[names(est) %in% names(se)]
        est <- est[anames]
        se <- se[anames]
    }
    tval <- as.vector(est)/se
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    mthd <- "z"
    rval <- cbind(est, se, tval, pval)
    colnames(rval) <- cnames
    class(rval) <- "coeftest"
    attr(rval, "method") <- paste(mthd, "test of coefficients")
    return(rval)
}
## weights.stpm2 <- 
## function (object, ...) 
## {
##     wts <- object@weights
##     if (is.null(wts)) 
##         wts
##     else napredict(object@na.action, wts)
## }

require(sandwich)
coxph1 <- coxph(Surv(rectime,censrec==1)~hormon,data=brcancer)
update(coxph1,robust=TRUE)
sandwich(coxph1)
sandwich.stpm2(fit) # hurrah!


## require(lmtest)
## coeftest(coxph1)
## coeftest(coxph1,vcov.=sandwich(coxph1))
## coeftest(fit,sandwich(fit))


sandwich(fit)
sandwich(fit,bread.=bread.stpm2,meat.=meat.stpm2)


## some predictions
head(predict(fit,se.fit=T,type="surv"))
head(predict(fit,se.fit=T,type="hazard"))

## some plots
plot(fit,newdata=data.frame(hormon=0),type="hazard")

## time-varying coefficient
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                         tvc=list(hormon=3)))
anova(fit,fit.tvc) # compare with and without tvc

plot(fit.tvc,newdata=data.frame(hormon=0),type="hr",var="hormon")
                                        # no lines method: use add=TRUE
plot(fit.tvc,newdata=data.frame(hormon=1),type="hr",var="hormon",
     add=TRUE,ci=FALSE,line.col=2)

plot(fit.tvc,newdata=data.frame(hormon=0),type="sdiff",var="hormon")

plot(fit.tvc,newdata=data.frame(hormon=0),type="hdiff",var="hormon")

plot(fit.tvc,newdata=data.frame(hormon=0),type="hazard")
plot(fit.tvc,newdata=data.frame(hormon=1),type="hazard",line.col=2,ci=FALSE,add=TRUE)
## trace("predict", browser, exit=browser, signature = "stpm2")

set.seed(10101)
brcancer <- transform(brcancer, x=rlnorm(nrow(brcancer)))
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                     tvc.formula=~hormon:nsx(log(rectime),df=3)))



## cure model
require(foreign)
colon <- read.dta("http://www.pauldickman.com/survival/colon.dta")
popmort <- read.dta("http://www.pauldickman.com/survival/popmort.dta")
names(popmort)[match(c("_age","_year"),names(popmort))] <- c("X_age","X_year")
colon2 <- transform(transform(colon,
                              status=ifelse(`surv_mm`>120.5,1,status),
                              tm=pmin(`surv_mm`,120.5)/12,
                              sex=as.numeric(sex)),
                    X_age=pmin(floor(age+tm),99),
                    X_year=floor(yydx+tm))
colon2 <- merge(colon2,popmort)
summary(fit <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     bhazard=colon2$rate,
                     logH.formula=~nsx(log(tm),df=5,cure=T,stata=T))) # oops
head(predict(fit))
plot(fit,newdata=data.frame(hormon=1))
plot(fit,newdata=data.frame(hormon=1),type="hazard")
plot(fit,newdata=data.frame(hormon=1),type="cumhaz")

oldx <- 0:100
oldy <- (oldx-50)^2
oldy[c(20,30)] <- 0
old <- data.frame(x=oldx,y=oldy)
predict(lm(y~nsx(x,knots=c(25,50,75,95)),old)) # as per Stata
newx <- seq(min(oldx)/1.05,max(oldx)*1.05,length=101)
new <- data.frame(x=newx)
plot(oldx,oldy)
predict(lm(y~nsx(x,df=5,cure=T),old))
sum(oldy)
terms(lm(y~nsx(x,df=5,cure=T),old))
lm(y~nsx(x,df=5),old)


lines(newx,
      predict(lm(y~nsx(x,df=4,cure=F),old),newdata=new),
      type="l") # oops
lines(newx,
      predict(lm(y~nsx(x,df=3),old),newdata=new),
      lty=2)


summary(fit <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     bhazard=colon2$rate,
                     logH.formula=~nsx(log(tm),df=6,stata=T))) # okay
summary(fit <- stpm2(Surv(tm,status %in% 2:3)~I(year8594=="Diagnosed 85-94"),
                     data=colon2,
                     logH.formula=~nsx(log(tm),df=6,stata=T))) # okay

## Stata
## stata.knots=c(4.276666164398193, 6.214608192443848, 6.7833251953125, 7.806289196014404)
stataKnots <- function(x,df) {
  intKnots <- round((1:(df-1))/df,2) # yes, Paul implicitly rounded to 2 dp
  logx <- log(x)
  c(min(logx),quantile(logx,intKnots,type=2),max(logx))
}
stata.knots <- stataKnots(subset(brcancer,censrec==1)$rectime,3)
## sapply(1:9,function(type) log(quantile(subset(brcancer,censrec==1)$rectime,c(0.33,0.67),type=type)))
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
                     logH.args=list(knots=stata.knots[2:3],
                       Boundary.knots=stata.knots[c(1,4)])))
## formula specification for logH
summary(stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,
              logH.formula=~ns(log(rectime),df=3)))

pred <- predict(fit.tvc,newdata=data.frame(hormon=0:3),grid=T,se.fit=T,type="cumhaz")
pred.all <- cbind(pred,attr(pred,"newdata"))
require(lattice)
xyplot(Estimate ~ rectime, data=pred.all, group=hormon,type="l",xlab="Time")


## relative survival
brcancer <- transform(brcancer,rate0=10^(-5+x1/100))
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,bhazard=brcancer$rate0,df=3))
head(predict(fit,se.fit=T))

## delayed entry
brcancer2 <- transform(brcancer,startTime=ifelse(hormon==0,rectime*0.5,0))
## debug(stpm2)
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     logH.formula=~nsx(log(rectime),df=3,stata=T)))
head(predict(fit,se.fit=T))
## delayed entry and tvc
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=brcancer2,
                     tvc.formula=~hormon:nsx(log(rectime),df=3,stata=T)))
head(predict(fit,se.fit=T))



## multiple time scales
brcancer <- transform(brcancer,recyr=rectime/365.25)
## predictions from a simple model
summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon+x1,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50))))
head(predict(fit))
grid.x1 <- with(brcancer, seq(40,70,length=300))
newdata0 <- with(brcancer, data.frame(recyr=5,x1=grid.x1,hormon=0))
matplot(grid.x1,
        predict(fit,type="hr",newdata=newdata0,var="hormon",se.fit=T), type="l")
## predictions with multiple time scales
summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50)),
                     tvc.formula=~hormon:nsx(log(recyr+x1),df=2)))
matplot(grid.x1,
        predict(fit,type="hr",newdata=newdata0,var="hormon",se.fit=T), type="l")





brcancer <- transform(brcancer,recyr=rectime/365.25,entry=recyr/2)
summary(fit <- stpm2(Surv(entry,recyr,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50)),
                     tvc.formula=~hormon:nsx(log(recyr+x1),df=2)))


summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon+x1,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50))))


plot(grid.x1,
     predict(fit,type="hr",newdata=newdata0,var="hormon",se.fit=T)$fit, type="l")

plot(fit,newdata=data.frame(hormon=0,x1=50),var="hormon",type="hr")

head(predict(fit,type="hazard",newdata=newdata0))
head(predict(fit,type="hazard",newdata=transform(newdata0,hormon=1)))



newdata0 <- with(brcancer, data.frame(recyr=5+1,x1=grid.x1-1,hormon=0))
predict(fit,type="hr",newdata=newdata0,var="hormon")

summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon+x1,data=brcancer,
                     logH.formula=~nsx(log(recyr),df=3,centre=log(50)),tvc=list(hormon=3)))

brcancer <- transform(brcancer, startAge=x1, endAge=x1+rectime/365)
summary(fit <- stpm2(Surv(startAge,endAge,censrec==1)~hormon,data=brcancer,
                     logH.formula=~nsx(log(endAge),df=3,centre=log(50)),tvc=list(hormon=3)))


## some simulated data: H_{weibull}(t)=(t/b)^a
n <- 1000
sim1 <- data.frame(age=seq(20,70,length=n),x=rep(0:1,each=n/2))
y <- rweibull(1000,shape=1,scale=1)



with(brcancer, plot(density(x1[censrec==1])))

summary(fit <- stpm2(Surv(recyr,censrec==1)~hormon,data=brcancer,logH.formula=~nsx(log(recyr),df=3,stata=T)))




brcancer <- transform(brcancer,ageStart=rnorm(length(rectime),50,5))
brcancer <- transform(brcancer,ageStop=ageStart+rectime)
summary(fit <- stpm2(Surv(ageStart,ageStop,censrec==1)~hormon,data=brcancer,df=3))

brcancer3 <- transform(brcancer,startTime=ifelse(censrec==1,0,10))
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=subset(brcancer,rectime>10),df=3))
summary(fit <- stpm2(Surv(startTime,rectime,censrec==1)~hormon,data=subset(brcancer3,rectime>10),df=3))

## check the performance time
brcancer10 = do.call("rbind",lapply(1:100,function(i) brcancer))
system.time(summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,df=3,data=brcancer10)))


nsx(1:10,df=3) - ns(1:10,df=3)
nsx(1:10,df=3,centre=3)
nsx(1:10,df=3,centre=3,Boundary.knots=c(2,8),derivs=c(1,1))
nsx(1:10,df=3,cure=T)
nsxDeriv(1:10,df=3) - nsDeriv(1:10,df=3)
nsxDeriv(1:10,df=3,centre=5,derivs=c(1,1))
nsxDeriv(1:10,df=3,centre=5,cure=T)

nsDeriv(1:10,df=3) - nsDeriv2(1:10,df=3)

## bug with calling mle2
require(bbmle)
mle2a <- function(...)
  mle2(...)
## some data
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)
## some fits
(fit0 <- mle2(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d)) # okay
(fit0.2 <- mle2(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d,
              control=list(parscale=2))) # okay
(fit1 <- mle2a(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d)) # okay
(fit1.2 <- mle2a(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d,
              control=list(parscale=2))) # FAILS

} ## end of examples ##




## ## * stata
## cd c:\Users\marcle\Documents\Home\
## clear
## *webuse brcancer
## use brcancer
## stset rectime, f(censrec==1)
## cap program drop dopredictions
## program define dopredictions
##   preserve
##   predict hr, hrnumerator(hormon 1) ci
##   predict haz, hazard ci
##   predict surv, surv ci
##   predict sdiff, sdiff1(hormon 1) ci
##   list hr* in 1/5
##   list haz* surv* in 1/5
##   list sdiff* in 1/5
##   restore
## end

## * basic model
## stpm2 hormon, df(3) scale(h)
## dopredictions

## * tvc
## stpm2 hormon, df(3) tvc(hormon) dftvc(3) scale(h)
## dopredictions

## * delayed entry
## preserve
##   replace _t0 = rectime*0.5 if hormon==0
##   stpm2 hormon, df(3) scale(h)
##   dopredictions
## restore

## * relative survival
## preserve  
##   gen rate0=10^(-5+x1/100)
##   stpm2 hormon, df(3) scale(h) bhazard(rate0)
##   dopredictions
## restore

## * test speed
## clear all
## set mem 100m
## use brcancer
## stset rectime, f(censrec==1)
## expand 100
## timer clear
## timer on 1
## stpm2 hormon, df(3) scale(h)
## timer off 1
## timer list
 

## hazard.pm = function(obj,tm,X,XD) # obj$par
## {
##   Xlocal=predict(X,newx=log(tm))
##   XDlocal=predict(XD,newx=log(tm))
##   with(obj,
##        c((XDlocal %*% par)/tm*exp(Xlocal %*% par)))
## }
## with(list(df=df,x=seq(0,3,length=100)[-1]),
##      {
##        plot(x,hazard.pm(fit,x,X,XD),type="l",ylim=c(0,2))
##        lines(x,dweibull(x,shape=1)/pweibull(x,shape=1,lower=F),lty=2)
##      })
## ##
## require(deSolve)
## temp <- as.data.frame(ode(y=0,times=seq(0,10,length=100)[-1],
##                           func=function(t,state,parameters=NULL) list(exp(sin(2*pi*log(t))))))
## plot(temp,type="l")
## temp <- transform(temp, cum=`1`,logcum=log(`1`))
## with(temp,plot(log(time),logcum))
## temp1 <- temp[-1,]
## fit <- glm(log(cum)~log(time)+sin(2*pi*log(time))+cos(2*pi*log(time)),data=temp1)
## lines(log(temp1$time),predict(fit))
## ## In summary:
## ## we can model using sine and cosine terms for the log-cumulative hazard - for log(time).

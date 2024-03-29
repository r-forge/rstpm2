\name{stpm2}
\Rdversion{1.1}
\alias{stpm2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Flexible parametric survival model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{

stpm2(formula, data, df=3, logH.args=NULL, logH.formula=NULL, tvc=NULL,
                  tvc.formula=NULL,
                  control=list(parscale=0.1,maxit=300), init=FALSE,
                  coxph.strata=NULL, weights=NULL, robust=FALSE,
                  bhazard=NULL, contrasts=NULL, subset=NULL, ...)  }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula object, with the response on the left of a \code{~} operator, and
the terms on the right.  The response must be a survival object as
returned by the \code{\link{Surv}} function. The terms should include linear
terms for any time-varying coefficients. [required]
}
  \item{data}{
a data.frame in which to interpret the variables named in
the \code{formula} argument. [at present: required]
}
  \item{df}{
an integer that describes the degrees of freedom for the \code{ns}
function for modelling the baseline log-cumulative hazard (default=3). 
}
  \item{logH.args}{
a list describing the arguments for the \code{ns} function for modelling
the baseline log-cumulative hazard (default=NULL). Use this for
changing the knot placement.
}
  \item{logH.formula}{
a formula for describing the baseline log-cumulative hazard function
(default=NULL). Only one of \code{df}, \code{logH.args} or
\code{logH.formula} is required.
}
  \item{tvc}{
a list with the names of the time-varying coefficients and the degrees
of freedom (e.g. \code{tvc=list(x=3)} specifies \code{x} as a time-varying
coefficient with 3 degrees of freedom).
}
  \item{tvc.formula}{
a formula for describing the time-varying coefficients. If a
time-varying coefficient is being model, then only one of
\code{tvc} and \code{tvc.formula} is required.
}
  \item{bhazard}{
a vector for the background hazard for relative survival estimation. At
present, this does not use \code{data}. At present, this is required for
all individuals - although it is only used at the event times.
}
\item{control}{\code{control} argument passed to \code{optim}.}
\item{init}{\code{init} should be \code{FALSE}, such that initial
                  values will be determined internally, or a numeric
                  vector of initial values.}
\item{coxph.strata}{variable in the \code{data} argument for
  stratification of the \code{coxph} model fit for estimating initial values.}
\item{weights}{an optional vector of 'prior weights' to be used in the
  fitting process. Should be \code{NULL} or a numeric vector.}
\item{robust}{Boolean used to determine whether to use a robust variance
  estimator.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{\link{model.matrix.default}}.
}
\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
\item{\dots}{
  additional arguments to be passed to the \code{\link{mle2}} .
}
}
\details{
  
  
  The implementation extends the \code{mle2} object from the
  \code{bbmle} package.
}
\value{
An \code{stpm2-class} object. 
  %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(brcancer)
summary(fit <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3))

## some predictions
head(predict(fit,se.fit=TRUE,type="surv"))
head(predict(fit,se.fit=TRUE,type="hazard"))

## some plots
plot(fit,newdata=data.frame(hormon=0),type="hazard")
plot(fit,newdata=data.frame(hormon=0),type="surv")

## time-varying coefficient
summary(fit.tvc <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=3,
                     tvc=list(hormon=3)))
anova(fit,fit.tvc) # compare with and without tvc

## some more plots
plot(fit.tvc,newdata=data.frame(hormon=0),type="hr",var="hormon")
                                        # no lines method: use add=TRUE
plot(fit.tvc,newdata=data.frame(hormon=1),type="hr",var="hormon",
     add=TRUE,ci=FALSE,line.col=2)

plot(fit.tvc,newdata=data.frame(hormon=0),type="sdiff",var="hormon")

plot(fit.tvc,newdata=data.frame(hormon=0),type="hdiff",var="hormon")

plot(fit.tvc,newdata=data.frame(hormon=0),type="hazard")
plot(fit.tvc,newdata=data.frame(hormon=1),type="hazard",line.col=2,ci=FALSE,add=TRUE)

## compare number of knots
hormon0 <- data.frame(hormon=0)
plot(fit,type="hazard",newdata=hormon0)
AIC(fit)
for (df in 4:6) {
    fit.new <- stpm2(Surv(rectime,censrec==1)~hormon,data=brcancer,df=df)
    plot(fit.new,type="hazard",newdata=hormon0,add=TRUE,ci=FALSE,line.col=df)
    print(AIC(fit.new))
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{stpm2-class}
\Rdversion{1.1}
\docType{class}
\alias{stpm2-class}
\alias{plot,stpm2,missing-method}
\alias{predict,stpm2-method}
\alias{predictnl,stpm2-method}

\title{Class "stpm2" ~~~ }
\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("stpm2", ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
	 \describe{
    \item{\code{xlevels}:}{Object of class \code{"list"} ~~ }
    \item{\code{contrasts}:}{Object of class \code{"listOrNULL"} ~~ }
    \item{\code{terms}:}{Object of class \code{"terms"} ~~ }
    \item{\code{model.frame}:}{Object of class \code{"list"} ~~ }
    \item{\code{x}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{xd}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{termsd}:}{Object of class \code{"terms"} ~~ }
    \item{\code{Call}:}{Object of class \code{"character"} ~~ }
    \item{\code{y}:}{Object of class \code{"Surv"} ~~ }
    \item{\code{call}:}{Object of class \code{"language"} ~~ }
    \item{\code{call.orig}:}{Object of class \code{"language"} ~~ }
    \item{\code{coef}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{fullcoef}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{vcov}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{min}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{details}:}{Object of class \code{"list"} ~~ }
    \item{\code{minuslogl}:}{Object of class \code{"function"} ~~ }
    \item{\code{method}:}{Object of class \code{"character"} ~~ }
    \item{\code{data}:}{Object of class \code{"list"} ~~ }
    \item{\code{formula}:}{Object of class \code{"character"} ~~ }
    \item{\code{optimizer}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{mle2}"}, directly.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "stpm2", y = "missing")}: ... }
    \item{predict}{\code{signature(object = "stpm2")}: ... }
    \item{predictnl}{\code{signature(object = "stpm2")}: ... }
	 }
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
	or \code{\linkS4class{CLASSNAME}} for links to other classes
}
\examples{
showClass("stpm2")
}
\keyword{classes}

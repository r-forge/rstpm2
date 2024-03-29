\name{predict-methods}
\docType{methods}
\alias{predict-methods}
\alias{predict,stpm2-method}
\title{ Predicted values for an stpm2 fit}
\description{
  Given an \code{stpm2} fit and an optional list of new data, return predictions
}
\section{Methods}{
  \describe{
    \item{object= "stpm2"}{an \code{stpm2} fit}
}}

\usage{

    \S4method{predict}{stpm2}(object, newdata=NULL,
               type=c("surv","cumhaz","hazard","hr","sdiff","hdiff","loghazard","link"),
               grid=FALSE,seqLength=300,
               se.fit=FALSE,link=NULL,exposed=incrVar(var),var,...)
}
\arguments{
  \item{object}{an \code{stpm2} object}
  \item{newdata}{optional list of new data (required if type in
               ("hr","sdiff","hdiff")). For type in
               ("hr","sdiff","hdiff"), this defines the unexposed
               newdata. This can be combined with \code{grid} to get a
               regular set of event times (i.e. newdata would \emph{not}
               include the event times). }
  \item{type}{specify the type of prediction:
    \itemize{
      \item{"surv"}{survival probabilities}
      \item{"cumhaz"}{cumulative hazard}
      \item{"hazard"}{hazard}
      \item{"hr"}{hazard ratio}
      \item{"sdiff"}{survival difference}
      \item{"hdiff"}{hazard difference}
      \item{"loghazard"}{log hazards}
      \item{"link"}{X beta}
    }
  }
  \item{grid}{whether to merge newdata with a regular sequence of event
    times (default=FALSE)}
  \item{seqLength}{length of the sequence used when \code{grid=TRUE}}
  \item{se.fit}{whether to calculate confidence intervals (default=FALSE)}
  \item{link}{allows a different link for the confidence interval
    calculation (default=NULL, such that
    switch(type,surv="cloglog",cumhaz="log",hazard="log",hr="log",sdiff="I",
    hdiff="I",loghazard="I",link="I"))}
  \item{exposed}{a function that takes newdata and returns a transformed
               data-frame for those exposed (defaults to incrementing ``var'')}
  \item{var}{specify the variable name or names for the
               exposed/unexposed (names are given as characters)}
  \item{\dots}{additional arguments (for generic compatibility)}
}

\value{
  A data-frame with components \code{Estimate}, \code{lower} and
  \code{upper}, with an attribute "newdata" for the \code{newdata}
  data-frame.  
}

\details{
  The confidence interval estimation is based on the delta method using
  numerical differentiation.
}

\seealso{\code{\link{stpm2}}}
\keyword{methods}
%%\keyword{ ~~ other possible keyword(s)}

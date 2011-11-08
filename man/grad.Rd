\name{grad}
\Rdversion{1.1}
\alias{grad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gradient function
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
grad(func, x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
%%     ~~Describe \code{func} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (func, x, ...) 
{
    h <- .Machine$double.eps^(1/3) * ifelse(abs(x) > 1, abs(x), 
        1)
    temp <- x + h
    h.hi <- temp - x
    temp <- x - h
    h.lo <- x - temp
    twoeps <- h.hi + h.lo
    nx <- length(x)
    ny <- length(func(x, ...))
    if (ny == 0L) 
        stop("Length of function equals 0")
    df <- if (ny == 1L) 
        rep(NA, nx)
    else matrix(NA, nrow = nx, ncol = ny)
    for (i in 1L:nx) {
        hi <- lo <- x
        hi[i] <- x[i] + h.hi[i]
        lo[i] <- x[i] - h.lo[i]
        if (ny == 1L) 
            df[i] <- (func(hi, ...) - func(lo, ...))/twoeps[i]
        else df[i, ] <- (func(hi, ...) - func(lo, ...))/twoeps[i]
    }
    return(df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

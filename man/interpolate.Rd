\name{interpolate}
\alias{interpolate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
interpolate(stand.bp, stand.time, sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stand.bp}{
%%     ~~Describe \code{stand.bp} here~~
}
  \item{stand.time}{
%%     ~~Describe \code{stand.time} here~~
}
  \item{sample}{
%%     ~~Describe \code{sample} here~~
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
function (stand.bp, stand.time, sample) 
{
    if (length(sample) > 0) {
        res <- NULL
        for (i in 1:length(sample)) {
            a <- sample[i] - stand.time
            if (length(unique(a > 0)) == 2) {
                pick <- which(a == min(a[a >= 0]))
                bp <- stand.bp[pick:(pick + 1)]
                time <- stand.time[pick:(pick + 1)]
                stand.lm <- lm(bp ~ time)
                res.i <- predict(stand.lm, newdata = data.frame(time = sample[i]))
                res[i] <- res.i
            }
            else {
                res[i] <- NA
            }
        }
        return(res)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

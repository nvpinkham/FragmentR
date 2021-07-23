\name{pairwise.dist}
\alias{pairwise.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pairwise.dist(str.1 = database[[1]], str.2, window = 2.5, method = "bc")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{str.1}{
%%     ~~Describe \code{str.1} here~~
}
  \item{str.2}{
%%     ~~Describe \code{str.2} here~~
}
  \item{window}{
%%     ~~Describe \code{window} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (str.1 = database[[1]], str.2, window = 2.5, method = "bc") 
{
    same.1 <- same.2 <- NULL
    for (i in 1:nrow(str.1)) {
        a <- abs(str.1$bp[i] - str.2$bp)
        if (min(a) <= window) {
            same.1 <- c(same.1, i)
            same.2 <- c(same.2, which(a == min(a))[1])
        }
    }
    str.1$bp[same.1] <- str.2$bp[same.2]
    h.1 <- str.1$height
    names(h.1) <- str.1$bp
    h.2 <- str.2$height
    names(h.2) <- str.2$bp
    h.11 <- aggregate(h.1, list(names(h.1)), sum)
    h.1 <- h.11$x
    names(h.1) <- h.11$Group.1
    h.21 <- aggregate(h.2, list(names(h.2)), sum)
    h.2 <- h.21$x
    names(h.2) <- h.21$Group.1
    aa <- names(h.2)[!names(h.2) \%in\% names(h.1)]
    bb <- vector(length = length(aa))
    names(bb) <- aa
    h.1 <- c(h.1, bb)
    h.1 <- h.1[order(names(h.1))]
    aa <- names(h.1)[!names(h.1) \%in\% names(h.2)]
    bb <- vector(length = length(aa))
    names(bb) <- aa
    h.2 <- c(h.2, bb)
    h.2 <- h.2[order(names(h.2))]
    hs <- rbind(h.1, h.2)
    hs <- hs * 25
    hs <- round(hs)
    set.seed(42)
    hs <- vegan::rrarefy(hs, min(rowSums(hs)))
    jd <- vegan::vegdist(hs, method = "jaccard")
    bc <- vegan::vegdist(hs)
    if (method == "bc") {
        return(bc)
    }
    if (method == "jd") {
        return(jd)
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

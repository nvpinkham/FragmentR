\name{find.match2}
\alias{find.match2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find.match2(file.query, ladder = rox.ladder, database = cd.db, meth = "bc", log.trans = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.query}{
%%     ~~Describe \code{file.query} here~~
}
  \item{ladder}{
%%     ~~Describe \code{ladder} here~~
}
  \item{database}{
%%     ~~Describe \code{database} here~~
}
  \item{meth}{
%%     ~~Describe \code{meth} here~~
}
  \item{log.trans}{
%%     ~~Describe \code{log.trans} here~~
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
function (file.query, ladder = rox.ladder, database = cd.db, 
    meth = "bc", log.trans = F) 
{
    bp.query <- sum.fsa(file.query, ladder, cutoff = 100, sd = 3)
    bp.query <- filter.bp(bp.query)
    if (log.trans) {
        bp.query$height <- log(bp.query$height)
        bp.query$area <- log(abs(bp.query$area))
    }
    if (nrow(bp.query) > 0) {
        res <- sapply(database, function(database) pairwise.dist(str.1 = database, 
            str.2 = bp.query, method = meth))
        names(res) <- names(database)
        return(sort(res)[1:10])
    }
    else {
        return("no match found")
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

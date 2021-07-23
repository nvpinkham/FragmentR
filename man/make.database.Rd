\name{make.database}
\alias{make.database}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.database(folder = "F-RibotypingFiles", ladder, cutoff = 100, sd = 3, bp.min = 200, bp.max = 600)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{folder}{
%%     ~~Describe \code{folder} here~~
}
  \item{ladder}{
%%     ~~Describe \code{ladder} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{sd}{
%%     ~~Describe \code{sd} here~~
}
  \item{bp.min}{
%%     ~~Describe \code{bp.min} here~~
}
  \item{bp.max}{
%%     ~~Describe \code{bp.max} here~~
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
function (folder = "F-RibotypingFiles", ladder, cutoff = 100, 
    sd = 3, bp.min = 200, bp.max = 600) 
{
    print("this may take a while please be patient")
    db.fsa <- list.files(folder, full.names = T, recursive = T, 
        pattern = ".fsa")
    str.bp <- sum.fsa(files[1], ladder, cutoff, sd)
    bp.query <- filter.bp(str.bp, bp.min, bp.max)
    summary.list <- lapply(db.fsa, sum.fsa, ladder, cutoff, sd)
    db <- lapply(db.fsa, function(db.fsa) sum.fsa(file_path = db.fsa, 
        cutoff = cutoff, ladder = rox.ladder, sd = sd))
    db.filtered <- lapply(db, filter.bp)
    names(db.filtered) <- db.fsa
    return(db.filtered)
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

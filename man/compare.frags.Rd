\name{compare.frags}
\alias{compare.frags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compare.frags(query.file = "Promega_001_96-well-plate_ready-to-go/001-G09.fsa", query.ladder = rox.ladder, query.channel = 1, query.channel.ladder = 4, hit.file = "F-RibotypingFiles//FP313/7440.F09.FP313.fsa", hit.ladder = rox.ladder, hit.channel = 1, hit.channel.ladder = 4, meth = "bc", cutoff = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{query.file}{
%%     ~~Describe \code{query.file} here~~
}
  \item{query.ladder}{
%%     ~~Describe \code{query.ladder} here~~
}
  \item{query.channel}{
%%     ~~Describe \code{query.channel} here~~
}
  \item{query.channel.ladder}{
%%     ~~Describe \code{query.channel.ladder} here~~
}
  \item{hit.file}{
%%     ~~Describe \code{hit.file} here~~
}
  \item{hit.ladder}{
%%     ~~Describe \code{hit.ladder} here~~
}
  \item{hit.channel}{
%%     ~~Describe \code{hit.channel} here~~
}
  \item{hit.channel.ladder}{
%%     ~~Describe \code{hit.channel.ladder} here~~
}
  \item{meth}{
%%     ~~Describe \code{meth} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
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
function (query.file = "Promega_001_96-well-plate_ready-to-go/001-G09.fsa", 
    query.ladder = rox.ladder, query.channel = 1, query.channel.ladder = 4, 
    hit.file = "F-RibotypingFiles//FP313/7440.F09.FP313.fsa", 
    hit.ladder = rox.ladder, hit.channel = 1, hit.channel.ladder = 4, 
    meth = "bc", cutoff = 100) 
{
    str.1 <- sum.fsa(file_path = query.file, ladder = query.ladder, 
        channel.query = query.channel, channel.ladder = query.channel.ladder, 
        cutoff = cutoff)
    str.1 <- filter.bp(str.1)
    str.2 <- sum.fsa(file_path = hit.file, ladder = query.ladder, 
        channel.query = hit.channel, channel.ladder = hit.channel.ladder, 
        cutoff = cutoff)
    str.2 <- filter.bp(str.2)
    cdiff1 <- MLPA::read.fsa(query.file)
    cdiff2 <- MLPA::read.fsa(hit.file)
    cdiff1 <- cdiff1[, query.channel]
    cdiff2 <- cdiff2[, hit.channel]
    cdiff1[cdiff1 < 0] <- 0
    cdiff2[cdiff2 < 0] <- 0
    bp.1 <- str.1$bp
    time.1 <- str.1$time
    lm.1 <- lm(bp.1 ~ time.1)
    str1.x <- predict(lm.1, newdata = data.frame(time.1 = 1:length(cdiff1)))
    cdiff1 <- (cdiff1/max(cdiff1[str1.x > 150 & str1.x < 650])) * 
        100
    plot(str1.x, cdiff1, type = "l", xlim = c(150, 650), ylim = c(0, 
        100), xlab = "bp", ylab = "peak intensity", lwd = 2)
    points(str1.x[str.1$time], cdiff1[str.1$time])
    text(str1.x[str.1$time], cdiff1[str.1$time], round(str.1$bp, 
        0), pos = 2)
    bp.2 <- str.2$bp
    time.2 <- str.2$time
    lm.2 <- lm(bp.2 ~ time.2)
    str2.x <- predict(lm.2, newdata = data.frame(time.2 = 1:length(cdiff2)))
    cdiff2 <- (cdiff2/max(cdiff2[str.2 > 150 & str.2 < 650])) * 
        100
    lines(str2.x, cdiff2, type = "l", col = 2, lwd = 2, lty = 1)
    points(str2.x[str.2$time], cdiff2[str.2$time], col = 2)
    text(str2.x[str.2$time], cdiff2[str.2$time], round(str.2$bp, 
        0), pos = 4, col = 2)
    dist <- pairwise.dist(str.1, str.2, method = meth)
    a <- strsplit(query.file, "/")[[1]]
    b <- strsplit(hit.file, "/")[[1]]
    ribo <- b[length(b) - 1]
    title(paste0(round(dist, 3), "BC distance from\n", a[length(a)], 
        "\nto refference ", b[length(b)], " (", ribo, ")"))
    legend("topright", lty = 1, col = c(1, 2), legend = c(a[length(a)], 
        b[length(b)]))
    res <- c(query.file, hit.file, dist, ribo)
    names(res) <- c("query.file", "hit.file", "BCdist", "ribotype")
    res <- as.matrix(res)
    res <- as.data.frame(t(res))
    res$BCdist <- as.numeric(res$BCdist)
    return(res)
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

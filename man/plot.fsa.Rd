\name{plot.fsa}
\alias{plot.fsa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.fsa(file_path, ladder, cutoff = 250, revise.cutoff = T, sd = 3, channel.ladder = 4, channel.query = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_path}{
%%     ~~Describe \code{file_path} here~~
}
  \item{ladder}{
%%     ~~Describe \code{ladder} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{revise.cutoff}{
%%     ~~Describe \code{revise.cutoff} here~~
}
  \item{sd}{
%%     ~~Describe \code{sd} here~~
}
  \item{channel.ladder}{
%%     ~~Describe \code{channel.ladder} here~~
}
  \item{channel.query}{
%%     ~~Describe \code{channel.query} here~~
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
function (file_path, ladder, cutoff = 250, revise.cutoff = T, 
    sd = 3, channel.ladder = 4, channel.query = 1) 
{
    ladder.id <- deparse(substitute(ladder))
    cdiff <- MLPA::read.fsa(file_path)
    plot(cdiff[, 1], type = "n", ylim = range(cdiff) * 1.05, 
        ylab = "", xlab = "")
    abline(h = cutoff, lty = 3)
    num.channels <- length(attr(cdiff, "colors"))
    res <- list()
    for (i in 1:num.channels) {
        channel.i <- cdiff[, i]
        channel.i.2 <- cdiff[, i]
        lines(channel.i, col = attr(cdiff, "colors")[i], type = "l")
        channel.i[channel.i < cutoff] <- cutoff
        time <- find_peaks(channel.i, 50, 50)
        height <- channel.i.2[time]
        if (revise.cutoff & i == channel.ladder) {
            tallest <- rev(sort(height))[1:length(ladder)]
            new.cutoff <- median(tallest) - (sd(tallest) * sd)
            channel.i[channel.i < new.cutoff] <- new.cutoff
            time <- find_peaks(channel.i, 150, 150)
            int <- diff(time)
            time <- time[c(T, int < median(int) * 2)]
            height <- channel.i.2[time]
            abline(h = new.cutoff, lty = 3)
            text(1000, new.cutoff, paste("revised cutoff", round(new.cutoff)), 
                pos = 3)
        }
        points(time, height, bg = attr(cdiff, "colors")[i], pch = 21, 
            col = i)
        area <- find_area(x = channel.i.2, peaks = time)
        res[[i]] <- as.data.frame(cbind(time, height, area))
    }
    names(res) <- names(attr(cdiff, "colors"))
    legend("topright", pch = 21, col = 1:num.channels, pt.bg = attr(cdiff, 
        "colors"), legend = paste(names(attr(cdiff, "colors")), 
        "  wavelength =", attr(cdiff, "wavelengths"), "NM"))
    legend("bottomright", paste("file path =", file_path))
    l.r <- rev(ladder)
    h.r <- rev(res[[channel.ladder]]$height)
    t.r <- rev(res[[channel.ladder]]$time)
    height <- rev(h.r[1:min(c(length(l.r), length(h.r)))])
    time <- rev(t.r[1:min(c(length(l.r), length(t.r)))])
    ladder <- rev(l.r[1:min(c(length(l.r), length(h.r)))])
    text(time, height, paste("        ", ladder)[1:length(ladder)], 
        srt = 45, cex = 0.8, col = "grey")
    text(time[length(time)], height[length(height)], ladder.id, 
        col = attr(cdiff, "colors")[channel.ladder], pos = 4)
    bp <- interpolate(stand.bp = ladder, stand.time = time, sample = res[[channel.query]]$time)
    res.bp <- cbind(res[[channel.query]], bp)
    text(res.bp$time, res.bp$height, round(res.bp$bp), font = 4, 
        pos = 3)
    return(res.bp)
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

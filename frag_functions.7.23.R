

packages.needed <- c("MLPA", "bayestestR", "vegan", "zip")

packages2install <- !packages.needed %in% rownames(installed.packages())
packages.needed <- packages.needed[packages2install]

packages.needed <-c(packages.needed)

install.packages(packages.needed)

require("MLPA")
require("bayestestR")
require("vegan")
require("zip")# needed for windows users

ladder <- c(50,  75, 100, 125, 150, 200, 250, 300, 350, 400, 450, 475, 500,
            550, 600, 650, 700, 750, 800, 850, 900, 950, 1000)

#' Add together two numbers
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
add <- function(x, y) {
  return(x + y)
}

scalar1 <- function(x) {x / sqrt(sum(x^2))}

Closest <- function (x, a, which = FALSE, na.rm = FALSE){
  if (na.rm)
    x <- x[!is.na(x)]
  mdist <- min(abs(x - a))
  if (is.na(mdist))
    res <- NA
  else {
    idx <- DescTools::IsZero(abs(x - a) - mdist)
    if (which == TRUE)
      res <- which(idx)
    else res <- x[idx]
  }
  return(res)
}

local.max <- function(a, b){
  return(a[b - 21 + find_peaks(a[(b - 20) : (b + 20)], m = 20, n = 20)])
}

find_peaks <- function (x, m=50, n=50){ # default m is set to %; look at the points +/- % in the vector;
  # I addin "n" so you can chnage the range you look both behind and in front of the point
  shape <- diff(sign(diff(x, na.pad = FALSE))) # use sign to show the dirrection of change
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1 #
    z <- ifelse(z > 0, z, 1)
    w <- i + n + 1 #
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

find_area <- function(x, peaks,  m=50, n=50){

  if(length(peaks) > 0){

    #incase oeaks are at begining or end of vwctor
    peaks <- peaks + 51
    x <- c(rep(0, 50), x, rep(0, 50))

    area <- NULL


    for(j in 1:length(peaks)){

      # plot(x[(peaks[j] - m) : (peaks[j] + n)])

      area[j] <- bayestestR::auc(1 : 101 , x[(peaks[j] - m) : (peaks[j] + n)], "trapezoid")

    }
    return(area)
  }
}

extrapolate <- function(stand.bp, stand.time, sample){

  stand.lm <- lm(stand.bp ~ stand.time)

  res <- predict(stand.lm,
                 newdata = data.frame(stand.time = sample))


  return(res)

}

interpolate <- function(stand.bp, stand.time, sample){

  if(length(sample) > 0){

    res <- NULL

    for(i in 1:length(sample)){

      a <- sample[i] - stand.time

      if(length(unique(a > 0)) == 2){

        pick <- which(a == min(a[a >= 0]))

        bp <- stand.bp[pick : (pick + 1)]
        time <- stand.time[pick : (pick + 1)]

        stand.lm <- lm(bp ~ time)

        res.i <- predict(stand.lm,
                         newdata = data.frame(time = sample[i]))

        #res[i] <- round(res.i)
        res[i] <- res.i

      }else{

        #print("can't interpolate")
        res[i] <- NA
      }

    }

    return(res)
  }
}

filter.bp <- function(res.bp, bp.min = 200,  bp.max = 600){

  res <- res.bp[res.bp$bp >= bp.min & res.bp$bp <= bp.max & !is.na(res.bp$bp) , ]
  return(res)
}

plot.fsa <- function(file_path,
                          ladder,
                          cutoff = 250,
                          revise.cutoff = T,
                     sd = 3,
                          channel.ladder = 4,
                          channel.query = 1){

  ladder.id <- deparse(substitute(ladder))
  cdiff <- MLPA::read.fsa(file_path)

  plot(cdiff[, 1], type = "n",
       ylim = range(cdiff) * 1.05,
       ylab = "", xlab = "")

  abline(h = cutoff, lty = 3)

  num.channels <- length(  attr(cdiff, "colors"))

  res <- list()# save output

  for(i in 1 : num.channels){

    channel.i <- cdiff[, i]
    channel.i.2 <- cdiff[, i]
    lines(channel.i,
          col = attr(cdiff, "colors")[i],
          type = "l")

    channel.i[channel.i <  cutoff] <-  cutoff
    time <- find_peaks(channel.i, 50, 50)
    height <- channel.i.2[time]

    if(revise.cutoff & i == channel.ladder){

     tallest <- rev(sort(height))[1:length(ladder)]
     new.cutoff <-  median(  tallest ) - (sd(  tallest ) * sd)

     channel.i[channel.i <  new.cutoff] <-  new.cutoff
     time <- find_peaks(channel.i, 150, 150)

     int <- diff(time)
     time <- time[c(T, int < median(int) * 2)]

     height <- channel.i.2[time]

     abline(h = new.cutoff, lty = 3)
     text(1000, new.cutoff,
          paste("revised cutoff", round(new.cutoff)), pos = 3)


    }

    points(time,   height,
           bg = attr(cdiff, "colors")[i],
           pch = 21,
           col = i)

    area <- find_area(x = channel.i.2, peaks = time)
    res[[i]] <- as.data.frame(cbind(time, height, area))

  }

  names(res) <- names(attr(cdiff, "colors"))

  legend("topright",
         pch = 21,
         col = 1:   num.channels,
         pt.bg = attr(cdiff, "colors"),
         legend = paste(names(attr(cdiff, "colors")),
                        "  wavelength =",
                        attr(cdiff, "wavelengths"),
                        "NM"))
  legend("bottomright", paste("file path =", file_path))

  l.r <- rev(ladder)
  h.r <- rev(res[[channel.ladder]]$height)
  t.r <- rev(res[[channel.ladder]]$time)

  height <- rev(h.r[1 : min(c(length(l.r), length(h.r)))])
  time <- rev(t.r[1 : min(c(length(l.r), length(t.r)))])
  ladder <- rev(l.r[1 : min(c(length(l.r), length(h.r)))])

  text(time, height,
       paste("        ", ladder)[1:length(ladder)],
       srt = 45, cex = .8, col = "grey")

  text(time[length(time)],
       height[length(height)],
       ladder.id,
       col = attr(cdiff, "colors")[channel.ladder],
       pos = 4)

  bp <- interpolate(stand.bp = ladder, stand.time = time, sample = res[[channel.query]]$time)

  # extrapolated bp:
  #  ext <- extrapolate(stand.bp = ladder, stand.time = time, sample = res[[channel.query]]$time)
  #  ext[!is.na(int)] <- NA

  res.bp <- cbind(res[[channel.query]], bp)
 # res.bp <- res.bp[res.bp$bp >= 200 & !is.na(res.bp$bp),]

  text(res.bp$time, res.bp$height,
       round(res.bp$bp), font = 4, pos = 3)

  return(res.bp)
}

sum.fsa <- function(file_path,
                    ladder,
                    cutoff = 250,
                    sd = 3,
                    revise.cutoff = T,
                    channel.ladder = 4,
                    channel.query = 1){

  cdiff <- MLPA::read.fsa(file_path)

  num.channels <- length(  attr(cdiff, "colors"))

  res <- list()# save output

  for(i in 1 : num.channels){

    channel.i <- cdiff[, i]
    channel.i.2 <- cdiff[, i]

    channel.i[channel.i <  cutoff] <-  cutoff
    time <- find_peaks(channel.i, 50, 50)
    height <- channel.i.2[time]

    if(revise.cutoff & i == channel.ladder){

      tallest <- rev(sort(height))[1:length(ladder)]
      new.cutoff <-  median(  tallest ) - (sd(  tallest ) * sd)


      channel.i[channel.i <  new.cutoff] <-  new.cutoff

      time <- find_peaks(channel.i, 150, 150)# larger window for ladder

      int <- diff(time)
      time <- time[c(T, int < median(int) * 2)]

      height <- channel.i.2[time]


    }

    area <- find_area(x = channel.i.2, peaks = time)
    res[[i]] <- as.data.frame(cbind(time, height, area))

  }

  names(res) <- names(attr(cdiff, "colors"))

  l.r <- rev(ladder)
  h.r <- rev(res[[channel.ladder]]$height)
  t.r <- rev(res[[channel.ladder]]$time)

  height <- rev(h.r[1 : min(c(length(l.r), length(h.r)))])
  time <- rev(t.r[1 : min(c(length(l.r), length(t.r)))])
  ladder <- rev(l.r[1 : min(c(length(l.r), length(h.r)))])

  bp <- interpolate(stand.bp = ladder, stand.time = time, sample = res[[channel.query]]$time)

  # extrapolated bp:
  #  ext <- extrapolate(stand.bp = ladder, stand.time = time, sample = res[[channel.query]]$time)
  #  ext[!is.na(int)] <- NA

  res.bp <- cbind(res[[channel.query]], bp)
 # res.bp <- res.bp[res.bp$bp >= 200 & !is.na(res.bp$bp),]

  return(res.bp)
}


pairwise.dist <- function(str.1 = database[[1]],  str.2, window = 2.5, method = "bc"){

  same.1 <- same.2 <- NULL

  for(i in 1:nrow(str.1)){

    a <- abs(str.1$bp[i] - str.2$bp)

    if(min(a) <= window){

      same.1 <- c(same.1, i)
      same.2 <- c(same.2, which(a == min(a))[1])
      # same.2 <- c(same.2, which(a <= 5))
    }
  }

  str.1$bp[same.1] <- str.2$bp[same.2]

  h.1 <- str.1$height
  names(h.1) <- str.1$bp

  h.2 <- str.2$height
  names(h.2) <- str.2$bp

  ######################## deduplicate names

  h.11 <- aggregate(h.1, list(names(h.1)), sum)
  h.1 <- h.11$x
  names(h.1) <- h.11$Group.1


  h.21 <- aggregate(h.2, list(names(h.2)), sum)
  h.2 <- h.21$x
  names(h.2) <- h.21$Group.1

  #######################

  aa <- names(h.2)[!names(h.2) %in% names(h.1) ]
  bb <- vector(length = length(aa))
  names(bb) <- aa
  h.1 <- c(h.1, bb)
  h.1 <- h.1[order(names(h.1))]

  aa <- names(h.1)[!names(h.1) %in% names(h.2) ]
  bb <- vector(length = length(aa))
  names(bb) <- aa
  h.2 <- c(h.2, bb)
  h.2 <- h.2[order(names(h.2))]

  hs <- rbind(h.1, h.2)
  hs <- hs * 25# needed for the log transformation, rarefy needs integers
  hs <- round(hs)

  set.seed(42)
  hs <- vegan::rrarefy(hs, min(rowSums(hs)))

  jd <- vegan::vegdist(hs, method = "jaccard")
  bc <- vegan::vegdist(hs)

  if(method == "bc"){

    return(bc)
  }
  if(method == "jd"){

    return(jd)
  }
}

compare.frags <- function(query.file =  "Promega_001_96-well-plate_ready-to-go/001-G09.fsa",
                          query.ladder = rox.ladder,
                          query.channel= 1,
                          query.channel.ladder = 4,

                          hit.file = "F-RibotypingFiles//FP313/7440.F09.FP313.fsa" ,
                          hit.ladder = rox.ladder,
                          hit.channel= 1,
                          hit.channel.ladder = 4,

                          meth = "bc", cutoff = 100){



  str.1 <- sum.fsa(file_path = query.file,
                   ladder = query.ladder,
                   channel.query = query.channel,
                   channel.ladder = query.channel.ladder,
                   cutoff = cutoff)
  str.1 <- filter.bp(str.1)

  str.2 <- sum.fsa(file_path = hit.file,
                   ladder = query.ladder,
                   channel.query = hit.channel,
                   channel.ladder = hit.channel.ladder,
                   cutoff = cutoff)
  str.2 <- filter.bp(str.2)

  cdiff1 <- MLPA::read.fsa(query.file)
  cdiff2 <- MLPA::read.fsa(hit.file)

  cdiff1 <-   cdiff1[, query.channel]
  cdiff2 <-   cdiff2[, hit.channel]

  cdiff1[cdiff1 < 0] <- 0
  cdiff2[cdiff2 < 0] <- 0

  bp.1 <- str.1$bp
  time.1 <- str.1$time

  lm.1 <- lm(bp.1 ~ time.1)
  str1.x <- predict(lm.1, newdata = data.frame(time.1 = 1 : length(cdiff1)))

  cdiff1 <- (cdiff1 / max(cdiff1[  str1.x > 150 &   str1.x < 650])) * 100
  # normaliz to peak value in the 150 to 650 bp range

  plot(str1.x, cdiff1,  type = "l", xlim = c(150, 650), ylim = c(0,100),
       xlab = "bp", ylab = "peak intensity", lwd = 2)

  points( str1.x[str.1$time],
          cdiff1[str.1$time])

  text( str1.x[str.1$time],
        cdiff1[str.1$time],
        round(str.1$bp, 0), pos = 2)

  bp.2 <- str.2$bp
  time.2 <- str.2$time

  lm.2 <- lm(bp.2 ~ time.2)
  str2.x <- predict(lm.2, newdata = data.frame(time.2 = 1 : length(cdiff2)))

  cdiff2 <- (cdiff2 / max(cdiff2[  str.2 > 150 &   str.2 < 650])) * 100

  lines(str2.x, cdiff2, type = "l", col = 2, lwd = 2, lty = 1)

  points( str2.x[str.2$time],
          cdiff2[str.2$time],
          col = 2)

  text( str2.x[str.2$time],
        cdiff2[str.2$time],
        round(str.2$bp, 0), pos = 4, col = 2)

  ################ make hit plot ################
  dist <- pairwise.dist(str.1, str.2, method = meth)

  a <- strsplit(query.file, "/")[[1]]
  b <- strsplit(hit.file, "/")[[1]]

  ribo <- b[length(b) - 1]# Check to make sure each folder is a ribotype

  title(paste0(round(dist, 3), "BC distance from\n",
               a[length(a)], "\nto refference ",
               b[length(b)], " (",
               ribo, ")"))

  legend("topright", lty =1, col = c(1,2),
         legend = c( a[length(a)],
                     b[length(b)]))

  res <- c(query.file, hit.file, dist, ribo)
  names(res) <- c("query.file", "hit.file", "BCdist", "ribotype")
  res <- as.matrix(res)
  res <- as.data.frame(t(res))
  res$BCdist <- as.numeric(res$BCdist)

  return(res)
}

find.match2 <- function(file.query, ladder = rox.ladder, database = cd.db, meth= "bc", log.trans = F){

  bp.query <- sum.fsa(file.query, ladder, cutoff = 100, sd = 3)
  bp.query <- filter.bp(bp.query)

  if(log.trans){
    bp.query$height <- log(bp.query$height)
    bp.query$area <- log(abs(bp.query$area))
  }

  if( nrow(bp.query) > 0 ){

    res <- sapply(database,  function(database) pairwise.dist(str.1 = database, str.2 = bp.query, method = meth))
    names(res) <- names(database)
    return(sort(res)[1:10])
  }else{
    return("no match found")
  }
}



make.database <- function(folder = "F-RibotypingFiles", ladder,
                          cutoff = 100, sd = 3, bp.min = 200, bp.max = 600){

  print("this may take a while please be patient")

  db.fsa <- list.files(folder, full.names = T, recursive = T, pattern = ".fsa")
  str.bp <- sum.fsa(files[1], ladder, cutoff, sd)
  bp.query <- filter.bp(str.bp, bp.min, bp.max)

  summary.list <- lapply(db.fsa, sum.fsa, ladder, cutoff, sd)

  db <- lapply(db.fsa, function(db.fsa) sum.fsa(file_path = db.fsa,
                                                cutoff = cutoff,
                                                ladder = rox.ladder,
                                                sd = sd))
  db.filtered <- lapply(db, filter.bp)
  names(db.filtered) <- db.fsa

  return(db.filtered)

}


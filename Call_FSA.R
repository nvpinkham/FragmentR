

source("frag_functions.6.7.R")

rox.ladder <- c(50,  75, 100, 125, 150, 200, 250, 300, 350, 400, 450, 475, 500,
                550, 600, 650, 700, 750, 800, 850, 900, 950, 1000)

# Restore the object
cd.db <- readRDS(file = "Cdiff_DB_list.5.4.rds")
cd.log <- cd.db

for(i in 1 : length(cd.db)){
  
  cd.log[[i]]$area <- log(abs(cd.db[[i]]$area) + 0.0000001 )
  cd.log[[i]]$height <- log(cd.db[[i]]$height + 0.0000001 )
}



dir.create("Files_to_analyzed")
dir.create("Files_analyzed")

p1 <- list.files("Files_to_analyzed/", 
                 recursive = T, full.names = T, pattern = ".fsa")

ids <- strsplit(p1, split = "[/]")
IDs <- NULL
#ids <- sapply(ids, "[[", 2)

ids <-  gsub("Files_to_analyzed", "",  p1)
ids <-  gsub("/", "_",  ids)

jpegs <- gsub("fsa", "jpeg", IDs)

results.dir <- paste("Results", Sys.Date())
dir.create(results.dir)

results <-  as.data.frame(matrix(ncol = 9, nrow = length(p1)))
colnames(results) <- c("query_file", "Confidence", 
                       "Hit_1",  "ribo_1",  "Dist_1", "Hit_log" , "Dist_log", 
                       "hit_1_jpeg", "chormatogram_jpeg")

results$query_file <- IDs
querys.db <- list()

for(i in 1 : nrow(results)){
  
  print(i)
  
  i.res <- find.match2(file.query = p1[i], ladder = rox.ladder, database = cd.log, log.trans = F)
  i.res.log <- find.match2(file.query = p1[i], ladder = rox.ladder, database = cd.log, log.trans = T)
  
  
  if(i.res[1] !=  "no match found"){
   
    hit.file <- paste0(results.dir, "/hit_",  jpegs[i])
    
    jpeg(    hit.file)
    
    res1 <- compare.frags(query.file = p1[i], hit.file = names(i.res[1]), log.trans = F)
    results$Hit_1[i] <- res1$hit.file
    results$ribo_1[i] <- res1$ribotype
    results$Dist_1[i] <- res1$BCdist
    results$hit_1_jpeg[i] <- hit.file
    
    dev.off()
  
    res.log <- compare.frags(query.file = p1[i], hit.file = names(i.res.log[1]), log.trans = T)
    results$Hit_log[i] <- res.log$hit.file
    results$Dist_log[i] <- res.log$BCdist
    
  }else{
    
    results$Hit_1[i] <- i.res 
 
  }
  
  chrom.file <- paste0(results.dir, "/chrom_",  jpegs[i])
  jpeg(  chrom.file)
  plot.fsa(file_path = p1[i], ladder = rox.ladder)
  dev.off()
  
  results$chormatogram_jpeg[i] <- chrom.file
  
  file.copy(from = p1[i], to = paste("Files_analyzed/", ids[i]), overwrite = T)
#  file.remove(from = p1[i])
  
}


results$Confidence[results$Dist_1 < 0.1] <- "good match"
results$Confidence[results$Dist_1 >=  0.1] <- "questionable match"
results$Confidence[results$Dist_1 >=  0.2] <- "poor match"


results <- results[order(results$Dist_1) , ]
write.csv(results, paste0(results.dir, "/SUMMARY.csv"))

zip::zip(zipfile = "Results.zip", files = results.dir)



packages.needed <- c("MLPA", "bayestestR", "vegan", "zip")

packages2install <- !packages.needed %in% rownames(installed.packages())
packages.needed <- packages.needed[packages2install]

packages.needed <-c(packages.needed)

install.packages(packages.needed)

require("MLPA")
require("bayestestR")
require("vegan")
require("zip")# needed for windows users

if(!file.exists("Files_to_analyze")){
  dir.create("Files_to_analyze")
}





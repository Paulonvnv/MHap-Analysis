
# if(!require(adegenet)){
#   options(
#     repos = c(
#       zkamvar = "https://zkamvar.r-universe.dev",
#       CRAN = "https://cloud.r-project.org"
#     )
#   )
#   install.packages("adegenet")
#   library(adegenet)
# }else{
#   library(adegenet)
# }

if(!require(ade4)){
  install.packages("ade4")
  library(ade4)
}else{
  library(ade4)
}

# if(!require(poppr)){
#   install.packages("poppr")
#   library(poppr)
# }else{
#   library(poppr)
# }

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}else{
  library(dplyr)
}

if(!require(magrittr)){
  install.packages("magrittr")
  library(magrittr)
}else{
  library(magrittr)
}

if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)
}else{
  library(tidyr)
}

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}else{
  library(ggplot2)
}

if(!require(cowplot)){
  install.packages("cowplot")
  library(cowplot)
}else{
  library(cowplot)
}


if(!require(vegan)){
  install.packages("vegan")
  library(vegan)
}else{
  library(vegan)
}

if(!require(parallel)){
  install.packages("parallel")
  library(parallel)
}else{
  library(parallel)
}

if(!require(ape)){
  install.packages("ape")
  library(ape)
}else{
  library(ape)
}

# if(!require(pegas)){
#   install.packages("pegas")
#   library(pegas)
# }else{
#   library(pegas)
# }

if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}else{
  library(RColorBrewer)
}

if(!require(Hmisc)){
  install.packages('Hmisc')
  library(Hmisc)
}else{
  library(Hmisc)
}

if(!require(ggpubr)){
  install.packages('ggpubr')
  library(ggpubr)
}else{
  library(ggpubr)
}

if(!require(doMC)){
  install.packages('doMC')
  library(doMC)
}else{
  library(doMC)
}

if(!require(svMisc)){
  install.packages('svMisc')
  library(svMisc)
}else{
  library(svMisc)
}

if(!require(Biostrings)){
  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("Biostrings")
  library(Biostrings)
}else{
  library(Biostrings)
}

if(!require(Rcpp)){
  install.packages('Rcpp')
  library(Rcpp)
}else{
  library(Rcpp)
}

if(!require(DT)){
  install.packages('DT')
  library(DT)
}else{
  library(DT)
}

if(!require(stringr)){
  install.packages('stringr')
  library(stringr)
}else{
  library(stringr)
}

# if(!require(XLConnect)){
#   install.packages('XLConnect')
#   library(XLConnect)
# }


if(!require(tmap)){
  install.packages('tmap')
  library(tmap)
}


if(!require(sp)){
  install.packages('sp')
  library(sp)
}



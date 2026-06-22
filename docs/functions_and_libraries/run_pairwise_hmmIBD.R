#!/bin/r env

library(argparse)
library(stringr)
# Parse arguments ----

parser = ArgumentParser()

# PATHs to references and functions

parser$add_argument("-wd", "--wd", 
                    help="Path to input and output files or folders")

parser$add_argument("-fd", "--fd", default = 'null',
                    help="Path to function files")


# Upload of raw data

# parser$add_argument("-ampseqt", "--ampseq_tsvfolder",  default = 'null',
#                     help="Name of the folder where all the slots of the ampseq object are stored in tsv files")

parser$add_argument("-locit", "--loci_tsvfolder",  default = 'null',
                    help="Name of the folder where all the slots of the loci object are stored in tsv files")

parser$add_argument("-rGenomet", "--rGenome_tsvfolder",  default = 'null',
                    help="Name of the folder where all the slots of the rGenome object are stored in tsv files")

# Paramters to export data and results

parser$add_argument("-o", "--out",
                    help="string that define the prefix to be used for naming all output files and headers in the reports.")


# Parameters for IBD and Connectivity report

## UGER parameters for IBD calculation using task arrays
### Only "nTasks must be defined in the json file, the other two parameters are defined automatically

parser$add_argument("-t", "--nTasks", default = 1,
                    help="Number of Tasks arrays to split the estimation of IBD using tasks arrays in UGER")

parser$add_argument("-tid", "--Task_id", default = 'null',
                    help="Tasks array ID defined by UGER")

parser$add_argument("-ibd_step", "--ibd_step", default = 'null',
                    help="Step of the estimation of IBD, ")

## Other parameters for IBD

parser$add_argument("-pairs", "--pairs", default = 'null',
                    help="string with the file name of the pairwise comparisons")

parser$add_argument("-nChunks", "--nChunks", default = 500,
                    help="Number of chunks to subdivide the pairwise comparisons to reduce the consumption of RAM memory")

parser$add_argument("-parallel", "--parallel", default = FALSE,
                    help="Boolean to allow parallelization of the estimation of IBD")


parser$add_argument("-max_k", "--max_k", default = 200,
                    help="maximum value for switching paramente k of the estimation of IBD")

parser$add_argument("-by", "--by", default = NULL,
                    help="Boolean to allow parallelization of the estimation of IBD")

parser$add_argument("-max_alleles", "--max_alleles", default = NULL,
                    help="Boolean to allow parallelization of the estimation of IBD")


# Defining and checking variables ----

print("starting to parse variables")
args = parser$parse_args()

## PATHs to references and functions----

# Working directory
wd = gsub("/$", "", args$wd)

print(paste0('wd: ', wd))
setwd(wd) 

# Tools or functions directory
fd = gsub("/$", "", args$fd)
print(paste0('fd: ', fd))


## Upload of raw data----

# Starting ampseq file in tsv format
# ampseq_tsvfolder = args$ampseq_tsvfolder
# ampseq_tsvfolder = if(ampseq_tsvfolder == 'null'){NULL}else{file.path(wd, args$ampseq_tsvfolder)}
# print(paste0('ampseq_tsvfolder: ', ampseq_tsvfolder))

# Starting loci file in tsv format
loci_tsvfolder = args$loci_tsvfolder
loci_tsvfolder = if(loci_tsvfolder == 'null'){NULL}else{file.path(wd, args$loci_tsvfolder)}
print(paste0('loci_tsvfolder: ', loci_tsvfolder))

# Starting loci file in tsv format
rGenome_tsvfolder = args$rGenome_tsvfolder
rGenome_tsvfolder = if(rGenome_tsvfolder == 'null'){NULL}else{file.path(wd, args$rGenome_tsvfolder)}
print(paste0('rGenome_tsvfolder: ', rGenome_tsvfolder))



## Parameters to export data and results----

# output pattern
output = args$out
print(paste0('output: ', output))

### Parameters for IBD and Connectivity report----

#### Parameters for UGER----
# nTasks
nTasks = as.integer(args$nTasks)
print(paste0('nTasks: ', nTasks))

# Task_id
Task_id = as.integer(args$Task_id)
print(paste0('Task_id: ', Task_id))

# ibd_step
ibd_step = as.character(args$ibd_step)
print(paste0('ibd_step: ', ibd_step))

#### Other IBD parameters----

# pairs
pairs = as.character(args$pairs)
if(pairs == 'null'){
  pairs_df = NULL
}else{
    pairs_df = read.table(file.path(wd, pairs), sep = '\t', header = T)}
print(paste0('pairs: ', pairs))

# parallel
parallel = as.logical(args$parallel)
print(paste0('parallel: ', parallel))

# nChunks
nChunks = as.integer(args$nChunks)
print(paste0('nChunks: ', nChunks))

# max_k
max_k = as.integer(args$max_k)
print(paste0('max_k: ', max_k))

# by
by = as.character(args$by)
by = if(by == 'null'){NULL}else{file.path(wd, by)}
print(paste0('by: ', by))

# max_alleles
max_alleles = as.integer(args$max_alleles)

if(is.na(max_alleles)){
  max_alleles = NULL
}

print(paste0('max_alleles: ', max_alleles))


print("All variables checked")

# Check packages and functions----

print("Loading libraies and functions")
source(file.path(fd,'amplseq_required_libraries.R'))
source(file.path(fd,'amplseq_functions.R'))

# Upload object and convert to loci object----

if(ibd_step == 'pairwise'){print("starting to upload genotypic data")
  # if(!is.null(ampseq_tsvfolder)){
  #   print("Uploading genotipic data in ampseq format from tsv file")
  #   ampseq_object = read_ampseq(file = ampseq_tsvfolder, format = 'tsv')
  #   
  #   loci_object = ampseq2loci(ampseq_object)
  #   
  # }else 
  if(!is.null(loci_tsvfolder)){
    print("Uploading genotipic data in loci format from tsv file")
    loci_object = read_loci(file = loci_tsvfolder, format = 'tsv')
    
  }else if(!is.null(rGenome_tsvfolder)){
    print("Uploading genotipic data in rGenome format from tsv file")
    rGenome_object = read_rGenome(file = rGenome_tsvfolder, format = 'tsv')
    
  }}


# Genetic Relatedness ----

if(nTasks > 1 & !is.null(nTasks)){
  
  if(ibd_step == 'pairwise'){
    
    # call hmmIBD and PCA functions from Rcpp
    sourceCpp(file.path(fd,'hmmloglikelihood.cpp'))
    
    if(!file.exists(file.path(wd, 'pairwise_ibd_chunks'))){
      system(paste0('mkdir ', file.path(wd, 'pairwise_ibd_chunks')))
    }
    
    chunks = round(seq(1, nChunks + 1, length.out = nTasks + 1))
    low = chunks[Task_id]
    high = chunks[Task_id + 1] - 1
    
    pairwise_relatedness = NULL
    
    for(w in low:high){
      start = Sys.time()
      pairwise_relatedness = rbind(pairwise_relatedness,
                                   pairwise_hmmIBD(loci_object, parallel = FALSE, w = w, n = nChunks, pairs = NULL, max_k = max_k, by = NULL, freq_table = NULL, max_alleles = NULL))
      time_diff = Sys.time() - start
      
      print(paste0('step ', w, ' done in ', time_diff, ' secs'))
      
    }
    
    write.csv(pairwise_relatedness,
              file.path(wd, 'pairwise_ibd_chunks', paste0('pairwise_ibd_chunk_', Task_id, '.csv')),
              quote = FALSE,
              row.names = FALSE)
    
  }else if(ibd_step == 'merge'){
    
    pairwise_relatedness = NULL
    
    for(file in list.files(file.path(wd, 'pairwise_ibd_chunks'))){
      
      pairwise_relatedness = rbind(pairwise_relatedness,
                                   read.csv(file.path(wd, 'pairwise_ibd_chunks', file)))
      
    }
    
    write.csv(pairwise_relatedness,
              file.path(wd, paste0(output, '_pairwise_ibd', '.csv')),
              quote = FALSE,
              row.names = FALSE)
    
  }
  
}else if(nTasks == 1 | is.null(nTasks)){
  
  # call hmmIBD and PCA functions from Rcpp
  sourceCpp(file.path(fd,'hmmloglikelihood.cpp'))
  
  pairwise_relatedness = NULL
  
  for(w in nChunks){
    start = Sys.time()
    pairwise_relatedness = rbind(pairwise_relatedness,
                                 pairwise_hmmIBD(loci_object, parallel = FALSE, w = w, n = nChunks, pairs = NULL, max_k = max_k, by = NULL, freq_table = NULL, max_alleles = NULL))
    time_diff = Sys.time() - start
    
    print(paste0('step ', w, ' done in ', time_diff, ' secs'))
    
  }
  
  write.csv(pairwise_relatedness,
            file.path(wd, paste0(output, '_pairwise_ibd', '.csv')),
            quote = FALSE,
            row.names = FALSE)
  
}



#!/bin/r env

library(argparse)
library(stringr)

parser = ArgumentParser()

parser$add_argument("-i", "--path_to_fq", 
                    help="Path to input to files and output files")

parser$add_argument("-o", "--output_folder", 
                    help="Path to function files and reference files")



args = parser$parse_args()

path_to_fq = args$path_to_fq
path_to_fq = gsub('/$', '', path_to_fq)

output_folder = args$output_folder

for(folder in list.files(path_to_fq)){
  
  for(file in list.files(file.path(path_to_fq, folder))){
    
    
    if(!file.exists(output_folder)){
      system(paste0('mkdir ', output_folder))
    }
    
    system(paste0('cp ', file.path(path_to_fq, folder, file), ' ', output_folder))
    
  }
  
}
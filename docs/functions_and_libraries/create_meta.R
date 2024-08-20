#!/bin/r env

library(argparse)
library(stringr)

parser = ArgumentParser()

parser$add_argument("-i", "--path_to_fq", 
                    help="Path to input to files and output files")

parser$add_argument("-o", "--outpout_file", 
                    help="Path to function files and reference files")

parser$add_argument("-fw", "--pattern_fw", 
                    help="name of the gzvcf file")

parser$add_argument("-rv", "--pattern_rv", 
                    help="Prefix of output files")

args = parser$parse_args()

path_to_fq = args$path_to_fq
outpout_file = args$outpout_file

pattern_fw = args$pattern_fw
pattern_rv = args$pattern_rv

meta_df = NULL
for(file in list.files(path_to_fq)){
  if(grepl(pattern_fw, file)){
    meta_df = rbind(meta_df, data.frame(id = gsub(pattern_fw, '', file),
                                        ip1 = file.path(path_to_fq,file),
                                        ip2 = file.path(path_to_fq,gsub(pattern_fw, pattern_rv, file))))
  }
}

write.table(meta_df, file = outpout_file, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)



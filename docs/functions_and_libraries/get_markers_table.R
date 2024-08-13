
#!/bin/r env

library(argparse)
library(stringr)

parser = ArgumentParser()

parser$add_argument("-i", "--input_file", 
                    help="Path to input to files and output files")

parser$add_argument("-o", "--outpout_file", 
                    help="Path to function files and reference files")

parser$add_argument("-h", "--header", 
                    help="boolean for header")

parser$add_argument("-s", "--separator", default = ',',
                    help="separator pattern")

parser$add_argument("-gff", "--ref_gff", default = 'null',
                    help="name of .gff file containing coordinates of genomic regions")

args = parser$parse_args()

input_file = args$path_to_fq
outpout_file = args$outpout_file

header = as.logical(args$header)
separator = args$separator

marker_table = read.table(input_file, header = header, sep = separator)

marker_table = marker_table[,1:4]

names(marker_table) = c("amplicon",
                        "chromosome",
                        "start",
                        "end")

marker_table[['pos']] = (marker_table[['start']] + marker_table[['end']])/2
marker_table[['length']] == marker_table[['end']] - marker_table[['start']] + 1
  
marker_table[["distance"]] = Inf

for(chromosome in levels(as.factor(marker_table[["chromosome"]]))){
  for(amplicon in 1:(nrow(marker_table[marker_table[["chromosome"]] == chromosome,])-1)){
    marker_table[marker_table[["chromosome"]] == chromosome,][amplicon, "distance"] = 
      marker_table[marker_table[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - marker_table[marker_table[["chromosome"]] == chromosome,][amplicon, "pos"]
  }
}

write.csv(marker_table, outpout_file, quote = FALSE, row.names = FALSE)



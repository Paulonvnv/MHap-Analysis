#!/bin/r env

library(argparse)
library(stringr)

parser = ArgumentParser()

parser$add_argument("-i", "--path_to_fq", 
                    help="Path to input to files and output files")

parser$add_argument("-o", "--output_path", 
                    help="Path to function files and reference files")

parser$add_argument("-fw", "--pattern_fw", 
                    help="name of the gzvcf file")

parser$add_argument("-rv", "--pattern_rv", 
                    help="Prefix of output files")

parser$add_argument("-nseqs", "--max_nseqs", 
                    help="Maximum number of sequences per fastq file")

parser$add_argument("-fsize", "--max_file_size", 
                    help="Maximum size of the file to be downsize. Files above the threshold will be omitted.")

parser$add_argument("-ow", "--overwrite", default = FALSE,
                    help="Boolean to overwrite output_path if it exist")

args = parser$parse_args()

path_to_fq = args$path_to_fq
output_path = args$output_path

pattern_fw = args$pattern_fw
pattern_rv = args$pattern_rv

max_nseqs = as.integer(args$max_nseqs)

max_file_size = as.integer(args$max_file_size)

overwrite = as.logical(args$overwrite)

# path_to_fq = '~/Documents/Github/MHap-Analysis/docs/data/Pviv_example/Run2/Fastq'
# output_path = '~/Documents/Github/MHap-Analysis/docs/data/Pviv_example/Run2/Fastq_downsized3'
# 
# pattern_fw = '_L001_R1_001.fastq.gz'
# pattern_rv = '_L001_R2_001.fastq.gz'
# 
# max_nseqs = 200000
# 
# overwrite = FALSE

if(!file.exists(output_path)){
  
  print(paste0('Creating output directory ', output_path))
  system(paste0('mkdir ', output_path))
  files = list.files(path_to_fq)
}else{
  
  if(overwrite){
    
    print(paste0('Output path already exist, existing files will be overwritten. To complement already exisisting files and avoid overwriting, use the following argument: -ow FALSE'))
    system(paste0('rm -r ', output_path))
    system(paste0('mkdir ', output_path))
    files = list.files(path_to_fq)
    
  }else{
    
    print(paste0('Output path already exist, only pending files will be downsize. To overwrite all exisisting files use the following argument: -ow TRUE'))
    files = list.files(path_to_fq)
    processed_files = list.files(output_path)
    incomplete_files = processed_files[grepl('fastq$', processed_files)]
    processed_files = processed_files[!grepl('fastq$', processed_files)]
    
    processed_fws = processed_files[grepl(pattern_fw, processed_files)]
    processed_rvs = processed_files[grepl(pattern_rv, processed_files)]
    
    incomplete_fws = processed_fws[!(gsub(pattern_fw, '', processed_fws) %in% gsub(pattern_rv, '', processed_rvs))]
    incomplete_rvs = processed_rvs[!(gsub(pattern_rv, '', processed_rvs) %in% gsub(pattern_fw, '', processed_fws))]
    
    if(length(incomplete_fws) > 0){
      incomplete_files = c(incomplete_files, incomplete_fws)
      processed_files = processed_files[processed_files != incomplete_fws]
    }
    
    if(length(incomplete_rvs) > 0){
      incomplete_files = c(incomplete_files, incomplete_rvs)
      processed_files = processed_files[processed_files != incomplete_rvs]
    }
    
    files = files[!(files %in% processed_files)]
    
    for(incomplete_file in incomplete_files){
      system(paste0('rm ', file.path(output_path, incomplete_file)))
    }
    
  }
  
}


for(file in files){
  if(grepl(pattern_fw, file)){
    
    print(paste0('Start working on file number ', which(files == file), ' out of ', length(files)))
    
    fastq_path_1 = file.path(path_to_fq,file)
    fastq_path_2 = file.path(path_to_fq,gsub(pattern_fw, pattern_rv, file))
    
    fw_size = file.size(fastq_path_1)
    
    if(fw_size < max_file_size){
      
      n_sequences  = as.numeric(system(
        paste0('zcat -f ', 
               fastq_path_1,
               ' | wc -l'), 
        intern = TRUE))/4
      
      if(n_sequences > 0 & !is.null(n_sequences) & !is.na(n_sequences)){
        
        print(paste0('There are ', n_sequences, ' sequences in the file ', file))
        
        fastqs1 = system(
          paste0('zcat -f ', 
                 fastq_path_1), 
          intern = TRUE)
        
        print('Forward sequences imported into R environment.')
        
        fastqs2 = system(
          paste0('zcat -f ', 
                 fastq_path_2), 
          intern = TRUE)
        
        print('Reverse sequences imported into R environment.')
        
        seq_names_pos = seq(1, (4*n_sequences), 4)
        
        seqs_pos = seq_names_pos + 1
        
        plus_pos = seq_names_pos + 2
        
        seq_qc_pos = seq_names_pos + 3
        
        seqs1 = fastqs1[seqs_pos]
        seqs2 = fastqs2[seqs_pos]
        
        no_Ns_pos = which(!grepl('N{4,}', seqs1) & !grepl('N{4,}', seqs2))
        
        print(paste0(length(no_Ns_pos), " sequences are no N's."))
        
        if(length(no_Ns_pos) > max_nseqs){
          
          print('The number of sequences is above the max limit, sequences will be randomly sampled.')
          positions = sort(sample(no_Ns_pos, max_nseqs, replace = F))
          
          filtered_pos = sort(c(seq_names_pos[positions],
                                seqs_pos[positions],
                                plus_pos[positions],
                                seq_qc_pos[positions]))
          
          downsized_fastqs1 = fastqs1[filtered_pos]
          downsized_fastqs2 = fastqs2[filtered_pos]
          
        }else{
          
          positions = no_Ns_pos
          
          filtered_pos = sort(c(seq_names_pos[positions],
                                seqs_pos[positions],
                                plus_pos[positions],
                                seq_qc_pos[positions]))
          
          downsized_fastqs1 = fastqs1[filtered_pos]
          downsized_fastqs2 = fastqs2[filtered_pos]
        }
        
        write(downsized_fastqs1, file = file.path(output_path, gsub('.gz$', '', file)))
        system(paste0('gzip ', file.path(output_path, gsub('.gz$', '', file))))
        
        write(downsized_fastqs2, file = file.path(output_path, gsub('.gz$', '', gsub(pattern_fw, pattern_rv, file))))
        system(paste0('gzip ', file.path(output_path, gsub('.gz$', '', gsub(pattern_fw, pattern_rv, file)))))
        
        rm(list = c(
          'file',
          'seq_names_pos',
          'seqs_pos',
          'plus_pos',
          'seq_qc_pos',
          'seqs1',
          'seqs2',
          'no_Ns_pos',
          'filtered_pos',
          'downsized_fastqs1',
          'downsized_fastqs2',
          'n_sequences'
        ))  
      }else{
        
        print(paste0('The file ', file, ' is empty'))
        system(paste0('cp ', fastq_path_1, ' ', file.path(output_path, file)))
        
      }
      
    }else{
     print(paste0('The size of the file is ', fw_size, ', and it is above the max_file_size threshold', max_file_size, 'This sample will be skiped')) 
    }
  }
}



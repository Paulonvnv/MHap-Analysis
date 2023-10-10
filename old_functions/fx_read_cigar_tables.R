
fx_read_cigar_tables = function(paths, sample_id_pattern){
  
  # Uploading data sets and merge CIGARs tables
  # and generate a metadata table for sample names, runs and sample types (controls vs samples of interest)
  cigar_tables = NULL
  metadata = NULL
  
  for(run in list.files(paths)){
    cigar_run = read.table(file.path(paths, run, "dada2/run_dada2/CIGARVariants.out.tsv"))
    names(cigar_run) = cigar_run[1,]
    samples = sapply(strsplit(cigar_run[-1,1], "_"), function(x) x[1])
    samples[!grepl(sample_id_pattern,samples)] = paste(samples[!grepl(sample_id_pattern,samples)], run, sep = "_")
    cigar_run = cigar_run[-1, -1]
    cigar_run = apply(cigar_run, 2, function(x) as.integer(x))
    rownames(cigar_run) = paste(1:length(samples), run, samples, sep = '/')
    cigar_tables[[run]] = cigar_run

    rm(list = c("run", "cigar_run", "samples"))
  }
  
  # Combining cigar tables
  
  cigar_table = NULL
  
  for(run in names(cigar_tables)){
    temp_cigar_table = data.frame(samples = rownames(cigar_tables[[run]]), cigar_tables[[run]])
    names(temp_cigar_table) = c("samples", colnames(cigar_tables[[run]]))
    temp_cigar_table %<>% pivot_longer(cols = all_of(names(temp_cigar_table)[-1]), names_to = "alleles", values_to = "counts")
    cigar_table = rbind(cigar_table, temp_cigar_table)
    rm(run)
  }
  
  cigar_table %<>% pivot_wider(names_from = "alleles", values_from = "counts")
  sample_names = gsub("^.+/","",cigar_table[['samples']])
  
  duplicated_samples = sample_names[duplicated(sample_names)]
  
  cigar_table_replicates = cigar_table[duplicated(sample_names),]
  cigar_table = cigar_table[!duplicated(sample_names),]
  
  for(sample in duplicated_samples){
    
    temp_data = rbind(cigar_table[grepl(sample,
                                        cigar_table$samples),],
                      cigar_table_replicates[grepl(sample,
                                                   cigar_table_replicates$samples),])
    
    temp_test = rowSums(temp_data[,-1], na.rm = T)
    
    cigar_table[grepl(sample,
                      cigar_table$samples),] = temp_data[which.max(temp_test),]
    
  }
  
  metadata = data.frame(samples = gsub("^.+/","",cigar_table[['samples']]),
                        run = gsub('/.+$', '',gsub("^\\d+/","",cigar_table[['samples']])),
                        order_in_plate = gsub("/.+$","",cigar_table[['samples']]), typeofSamp = NA)
  
  
  metadata[grepl(sample_id_pattern, metadata[["samples"]]),][["typeofSamp"]] = "Samples"
  metadata[!grepl(sample_id_pattern, metadata[["samples"]]),][["typeofSamp"]] = "Controls"
  
  

  cigar_table = cigar_table[, -1]
  cigar_table = as.matrix(cigar_table)
  
  rownames(cigar_table) = metadata$samples
  
  cigar_table[is.na(cigar_table)] = 0
  
  return(list(cigar_table = cigar_table,
              metadata = metadata))
  
}

# Functions to create, upload, joint, write or convert formats that store genetic data ----

## cigar S4class and create_cigar----

## cigar S4 class

setClass('cigar', slots = c(
  cigar_table = "ANY",
  metadata = "ANY",
  asv_table = 'ANY',
  asv_seqs = 'ANY'
))

## cigar constructor
create_cigar = function(cigar_table = NULL,
                   metadata = NULL,
                   asv_table = NULL,
                   asv_seqs = NULL){
  
  obj = new('cigar')
  obj@cigar_table = cigar_table
  obj@metadata = metadata
  obj@asv_table = asv_table
  obj@asv_seqs = asv_seqs
  
  return(obj)
}


## read_cigar_tables----

read_cigar_tables = function(paths = NULL,
                             cigar_files = NULL,
                             asv_table_files = NULL,
                             asv2cigar_files = NULL,
                             asv_seq_files = NULL,
                             zero_read_sample_list = NULL,
                             sample_id_pattern = '.',
                             mask_Ns = TRUE){
  
  # Uploading data sets and merge CIGARs tables
  # and generate a metadata table for sample names, runs and sample types (controls vs samples of interest)
  cigar_tables = NULL
  metadata = NULL
  
  if(!is.null(paths)&is.null(cigar_files)){
    for(Run in list.files(paths)){
      
      cigar_run = NULL
      asv_table_run = NULL
      asv_seqs = NULL
      
      if(file.exists(file.path(paths, Run, "dada2/run_dada2/CIGARVariants.out.tsv"))){
        cigar_run = read.table(file.path(paths, Run, "dada2/run_dada2/CIGARVariants.out.tsv"), header = T, check.names = FALSE)
        
      }else if(file.exists(file.path(paths, Run, "dada2/run_dada2/CIGARVariants_Bfilter.out.tsv"))){
        cigar_run = read.table(file.path(paths, Run, "dada2/run_dada2/CIGARVariants_Bfilter.out.tsv"), header = T, check.names = FALSE)
      }else{
        print(paste0('Cigar file ', file.path(paths, Run, "dada2/run_dada2/CIGARVariants_Bfilter.out.tsv"), ' not found'))
      }
      
      if(file.exists(file.path(paths, Run, "dada2/run_dada2/ASV_to_CIGAR.out.txt"))){
        asv2cigar_run = read.table(file.path(paths, Run, "dada2/run_dada2/ASV_to_CIGAR.out.txt"),
                                   header = T,
                                   check.names = FALSE)
        
      }else{
        print(paste0('asv2cigar file ', file.path(paths, Run, "dada2/run_dada2/ASVTable.txt"), ' not found'))
      }
      
      if(file.exists(file.path(paths, Run, "dada2/run_dada2/ASVTable.txt"))){
        asv_table_run = read.table(file.path(paths, Run, "dada2/run_dada2/ASVTable.txt"),
                                   header = T,
                                   check.names = FALSE)
        asv_table_run = left_join(asv_table_run, asv2cigar_run, by = join_by(hapid == ASV))
        
      }else{
        print(paste0('asv_table file ', file.path(paths, Run, "dada2/run_dada2/ASVTable.txt"), ' not found'))
      }
      
      if(file.exists(file.path(paths, Run, "dada2/run_dada2/ASVSeqs.fasta"))){
        asv_seqs = Biostrings::readDNAStringSet(file.path(paths, Run, "dada2/run_dada2/ASVSeqs.fasta"))
      }else{
        print(paste0('asv_seqs file ', file.path(paths, Run, "dada2/run_dada2/ASVSeqs.fasta"), ' not found'))
      }
      
      if(file.exists(file.path(paths, Run, "dada2/run_dada2/zeroReadSamples.txt"))){
        ZeroReadSamples = read.table(file.path(paths, Run, "dada2/run_dada2/zeroReadSamples.txt") , header = T)
      }else{
        print(paste0('asv_seqs file ', file.path(paths, Run, "dada2/run_dada2/zeroReadSamples.txt"), ' not found'))
      }
      
      # Formating CIGAR table
      
      if(!is.null(cigar_run)){
        samples = cigar_run[,1]
        samples[!grepl(sample_id_pattern,samples)] = paste(samples[!grepl(sample_id_pattern,samples)], Run, sep = "_")
        cigar_run = cigar_run[, -1]
        cigar_run = apply(cigar_run, 2, function(x) as.integer(x))
        rownames(cigar_run) = paste(1:length(samples), Run, samples, sep = '/')
        cigar_tables[[Run]][['cigar_table']] = cigar_run
        cigar_tables[[Run]][['asv_table']] = asv_table_run
        cigar_tables[[Run]][['asv_seqs']] = asv_seqs
        
        if(length(ZeroReadSamples[[1]]) > 0 ){
          
          ZeroReadSamples_names = gsub('_prim_1.fq.gz$', '', ZeroReadSamples[[1]])
          ZeroReadSamples_names[!grepl(sample_id_pattern,ZeroReadSamples_names)] = paste(ZeroReadSamples_names[!grepl(sample_id_pattern,ZeroReadSamples_names)], Run, sep = "_")
          
          ZeroReadSamples = matrix(0, nrow = length(ZeroReadSamples[[1]]),
                                   ncol = ncol(cigar_tables[[Run]][['cigar_table']]),
                                   dimnames = list(ZeroReadSamples_names,
                                                   colnames(cigar_tables[[Run]][['cigar_table']]))
          )
          
          rownames(ZeroReadSamples) = paste(1:length(ZeroReadSamples_names), Run, ZeroReadSamples_names, sep = '/')
          
          cigar_tables[[Run]][['cigar_table']] = rbind(cigar_tables[[Run]][['cigar_table']],
                                                       ZeroReadSamples)
          
        }
        
        
        #rm(list = c("asv2cigar_run", "samples"))
        
      }
      
      rm(list = c("Run", "cigar_run", "asv_table_run", "asv_seqs"))
      
    }
  }else if(is.null(paths)&!is.null(cigar_files)){
    
    
    for(file in 1:length(cigar_files)){
      cigar_run = read.table(cigar_files[file], header = T, check.names = FALSE
      )
      #samples = gsub('_S\\d+$','', cigar_run[,1])
      samples = cigar_run[,1]
      samples[!grepl(sample_id_pattern,samples)] = paste0(samples[!grepl(sample_id_pattern,samples)], '_file', which(cigar_files == cigar_files[file]))
      cigar_run = cigar_run[, -1]
      cigar_run = apply(cigar_run, 2, function(x) as.integer(x))
      rownames(cigar_run) = paste(1:length(samples), cigar_files[file], samples, sep = '/')
      cigar_tables[[cigar_files[file]]][['cigar_table']] = cigar_run
      
      if(!is.null(asv_table_files) | !is.null(paths) & !is.null(asv2cigar_files) & !is.null(asv_seq_files) & !is.null(zero_read_sample_list)){
        
        if(length(asv_table_files) == length(cigar_files) &
           length(asv2cigar_files) == length(cigar_files) & 
           length(asv_seq_files) == length(cigar_files) &
           length(zero_read_sample_list == length(cigar_files))){
          
          # Uploading asv2cigar
          asv2cigar_run = read.table(asv2cigar_files[file],
                                     header = T,
                                     check.names = FALSE)
          
          # Uploading asv_table
          asv_table_run = read.table(asv_table_files[file],
                                     header = T,
                                     check.names = FALSE)
          
          asv_table_run = left_join(asv_table_run, asv2cigar_run, by = join_by(hapid == ASV))
          
          # Saving asv_table in cigar_tables object
          cigar_tables[[cigar_files[file]]][['asv_table']] = asv_table_run
          
          # Uploading asv_seqs
          asv_seqs = Biostrings::readDNAStringSet(asv_seq_files[file])
          
          # Saving asv_seqs in cigar_tables object
          cigar_tables[[cigar_files[file]]][['asv_seqs']] = asv_seqs
          
          # adding ZeroReadSamples
          
          ZeroReadSamples = read.table(zero_read_sample_list[file], header = T)
          
          ZeroReadSamples = matrix(0, nrow = length(ZeroReadSamples[[1]]),
                                   ncol = ncol(cigar_tables[[cigar_files[file]]][['cigar_table']]),
                                   dimnames = list(paste((nrow(cigar_tables[[cigar_files[file]]][['cigar_table']]) + 
                                                            1):(nrow(cigar_tables[[cigar_files[file]]][['cigar_table']]) + 
                                                                  length(ZeroReadSamples[[1]])), 
                                                         cigar_files[file], 
                                                         gsub('_prim_1.fq.gz$', '', ZeroReadSamples[[1]]), sep = '/'),
                                                   
                                                   colnames(cigar_tables[[cigar_files[file]]][['cigar_table']]))
          )
          
          cigar_tables[[cigar_files[file]]][['cigar_table']] = rbind(cigar_tables[[cigar_files[file]]][['cigar_table']],
                                                                     ZeroReadSamples
          )
          
          
          #rm(list = c("asv2cigar_run", "asv_table_run", "asv_seqs"))
          
        }else{
          
          print('The length of the following list of files asv_table_files, asv2cigar_files, and asv_seq_files, do not concide with the length of the cigar_files. The final cigar_object will not contain these slots')
          
        }
        
        
      }else{
        
        print('One of the following arguments: asv_table_files, asv2cigar_files, asv_seqs_files, has not been provided. The final cigar_object will not contain these slots')
      }
      
      rm(list = c("file", "cigar_run", "samples"))
    }
  }
  
  # Combining cigar tables
  
  cigar_table = NULL
  asv_table = NULL
  asv_seqs = NULL
  
  for(Run in names(cigar_tables)){
    temp_cigar_table = data.frame(Sample_id = rownames(cigar_tables[[Run]][['cigar_table']]), cigar_tables[[Run]][['cigar_table']])
    names(temp_cigar_table) = c("Sample_id", colnames(cigar_tables[[Run]][['cigar_table']]))
    temp_cigar_table %<>% pivot_longer(cols = all_of(names(temp_cigar_table)[-1]), names_to = "alleles", values_to = "counts")
    cigar_table = rbind(cigar_table, temp_cigar_table)
    
    if(!is.null(asv_table_files) | !is.null(paths)){
      temp_asv_seqs1 = cigar_tables[[Run]][['asv_seqs']]
      temp_asv_table1 = cigar_tables[[Run]][['asv_table']]
      
      # if the are shared ASVs between the previous and new data set, 
      if(sum(as.character(temp_asv_seqs1) %in% as.character(asv_seqs)) > 0){
        
        # Detect common ASVs from the previous and the new cigar table and sum up read counts, sample counts and bimera detection
        
        # asvs in the new data set present in the previous data set
        temp_asv_table2 = temp_asv_table1[as.character(temp_asv_seqs1) %in% as.character(asv_seqs),]
        temp_asv_seqs2 = temp_asv_seqs1[as.character(temp_asv_seqs1) %in% as.character(asv_seqs)]
        rownames(temp_asv_table2) = as.character(temp_asv_seqs2)
        
        # asvs in the previous data set present in the new data set
        temp_asv_table3 = asv_table[as.character(asv_seqs) %in% as.character(temp_asv_seqs1),]
        temp_asv_seqs3 = asv_seqs[as.character(asv_seqs) %in% as.character(temp_asv_seqs1)]
        rownames(temp_asv_table3) = as.character(temp_asv_seqs3)
        
        # Update the unique identifier of the asv (hapid)
        temp_asv_table2 = temp_asv_table2[rownames(temp_asv_table3),]
        temp_asv_table2[['hapid']] = temp_asv_table3[['hapid']] 
        
        # Update the read counts, the sample counts, and the bimera status
        temp_asv_table3[['total_reads']] = temp_asv_table3[['total_reads']] + temp_asv_table2[['total_reads']]
        temp_asv_table3[['total_samples']] = temp_asv_table3[['total_samples']] + temp_asv_table2[['total_samples']]
        temp_asv_table3[['bimera']] = temp_asv_table3[['bimera']] | temp_asv_table2[['bimera']]
        
        # Impute the inconsistent cigar_strings in the cigar_table
        unconsitent_cigar_strings = data.frame(temp_asv_table2[temp_asv_table2$CIGAR != temp_asv_table3$CIGAR,], CIGAR2 = temp_asv_table3[temp_asv_table3$CIGAR != temp_asv_table2$CIGAR,][['CIGAR']])
        
        for(allele in unconsitent_cigar_strings$CIGAR){
          
          cigar_table[cigar_table$alleles == paste(unconsitent_cigar_strings[unconsitent_cigar_strings$CIGAR == allele, c('Amplicon', 'CIGAR')], collapse = ','),][['alleles']] = 
            paste(unconsitent_cigar_strings[unconsitent_cigar_strings$CIGAR == allele, c('Amplicon', 'CIGAR2')], collapse = ',')
          
        }
        
        # ASVs in the previous data set not present in the new one
        temp_asv_table4 = asv_table[!(as.character(asv_seqs) %in% as.character(temp_asv_seqs1)),]
        
        asv_table = rbind(temp_asv_table4, temp_asv_table3)
        rownames(asv_table) = gsub('ASV','',asv_table$hapid)
        
      }
      
      temp_asv_table1 = temp_asv_table1[!(as.character(temp_asv_seqs1) %in% as.character(asv_seqs)),]
      temp_asv_seqs1 = temp_asv_seqs1[!(as.character(temp_asv_seqs1) %in% as.character(asv_seqs))]
      
      names(temp_asv_seqs1) = paste0('ASV',(1 + length(asv_seqs)):(length(temp_asv_seqs1) + length(asv_seqs)))
      temp_asv_table1$hapid = names(temp_asv_seqs1)
      
      asv_seqs = DNAStringSet(c(as.character(asv_seqs), as.character(temp_asv_seqs1)))
      
      asv_table = rbind(asv_table, temp_asv_table1)
      
      rm(list = c('Run', 'temp_cigar_table', 'temp_asv_seqs1', 'temp_asv_table1'))
    }else{
      rm(list = c('Run', 'temp_cigar_table'))
    }
    
    
  }
  
  if(!is.null(asv_table_files) | !is.null(paths)){
    asv_table[['CIGAR_masked']] = asv_table[['CIGAR']]
  }
  
  if(mask_Ns){
    
    if(!is.null(asv_table_files) | !is.null(paths)){
      
      asv_table[['CIGAR_masked']] = gsub('(\\d+N|\\d+(I|D)=[ATGC]N+)', '', asv_table[['CIGAR_masked']])
      
      asv_table %<>% mutate(CIGAR_masked = case_when(
        CIGAR_masked == '' & !is.na(CIGAR_masked) ~ '.',
        CIGAR_masked != '' & !is.na(CIGAR_masked) ~ CIGAR_masked
      ))
      
      cigar_strings = apply(asv_table[!is.na(asv_table$CIGAR_masked) &
                                        grepl('N', asv_table$CIGAR),], 1, function(cigar_string){
                                          paste(cigar_string['Amplicon'], cigar_string['CIGAR'], sep = ',')
                                        })
      
      if(length(cigar_strings) > 0){
        
        cigar_strings_replacements = apply(asv_table[!is.na(asv_table$CIGAR_masked) &
                                                       grepl('N', asv_table$CIGAR),], 1, function(cigar_string){
                                                         paste(cigar_string['Amplicon'], cigar_string['CIGAR_masked'], sep = ',')
                                                       })
        
        
        for(cigar_string in 1:length(cigar_strings)){
          
          cigar_table$alleles = gsub(paste0('^',cigar_strings[cigar_string], '$'),
                                     cigar_strings_replacements[cigar_string],
                                     cigar_table$alleles
          )
        }
        
      }
      
    }
    
    if(sum(grepl('(\\d+N{1}|\\d+(I|D)=[ATGC]N+)', cigar_table$alleles)) > 0){
      cigar_strings2 = unique(cigar_table$alleles[grepl('(\\d+N{1}|\\d+(I|D)=[ATGC]N+)', cigar_table$alleles)])
      
      cigar_strings_replacements2 = gsub('^.+,', '', cigar_strings2)
      cigar_strings_replacements2 = gsub('(\\d+N|\\d+(I|D)=[ATGC]N+)', '', cigar_strings_replacements2)
      
      cigar_strings_replacements2[cigar_strings_replacements2 == ''] = '.'
      
      cigar_strings_replacements2 = paste(gsub(',.+$', '', cigar_strings2),
                                          cigar_strings_replacements2,
                                          sep = ','
      )
      
      for(cigar_string in 1:length(cigar_strings2)){
        
        cigar_table$alleles = gsub(paste0('^',cigar_strings2[cigar_string], '$'),
                                   cigar_strings_replacements2[cigar_string],
                                   cigar_table$alleles
        )
      }
      
    }
    
    cigar_table %<>% summarise(counts = sum(counts), .by = c(Sample_id, alleles))
    
  }
  
  
  cigar_table %<>% pivot_wider(names_from = "alleles", values_from = "counts")
  sample_names = gsub("^.+/|_S\\d+$","",cigar_table[['Sample_id']])
  
  duplicated_samples = sample_names[duplicated(sample_names)]
  
  cigar_table_replicates = cigar_table[duplicated(sample_names),]
  cigar_table = cigar_table[!duplicated(sample_names),]
  
  for(sample in duplicated_samples){
    
    temp_data = rbind(cigar_table[grepl(sample,
                                        cigar_table$Sample_id),],
                      cigar_table_replicates[grepl(sample,
                                                   cigar_table_replicates$Sample_id),])
    
    temp_test = rowSums(temp_data[,-1], na.rm = T)
    
    cigar_table[grepl(sample,
                      cigar_table$Sample_id),] = temp_data[which.max(temp_test),]
    
  }
  
  metadata = data.frame(Sample_id = gsub("^.+/|_S\\d+$","",cigar_table[['Sample_id']]),
                        order_in_plate = as.integer(gsub("^S|_$","",str_extract(cigar_table[['Sample_id']], "(S\\d+$|S\\d+_)"))), typeofSamp = NA)
  
  if(!is.null(paths)){
    
    metadata[['Run']] = gsub('/.+$', '',
                             gsub("^\\d+/",
                                  "",
                                  cigar_table[['Sample_id']]))
    
  }else if(!is.null(cigar_files) & length(cigar_files) > 1){
    
    metadata[['Run']] = gsub('(^/|_?CIGARVariants|_Bfilter_)',
                             '',
                             str_extract(cigar_table[['Sample_id']], 
                                         pattern = '/(\\w+_)?CIGARVariants(_Bfilter_)?(w+|\\d+)?'))
    
  }else if(!is.null(cigar_files) & length(cigar_files) == 1){
    
    metadata[['Run']] = 'Run1'
    
  }
  
  
  if(sample_id_pattern != '.'){
    metadata[grepl(sample_id_pattern, metadata[["Sample_id"]]),][["typeofSamp"]] = "Samples"
    metadata[!grepl(sample_id_pattern, metadata[["Sample_id"]]),][["typeofSamp"]] = "Controls"
  }else{
    metadata[["typeofSamp"]] = "Samples"
  }
  
  cigar_table = cigar_table[, -1]
  cigar_table = as.matrix(cigar_table)
  
  rownames(cigar_table) = metadata$Sample_id
  
  cigar_table[is.na(cigar_table)] = 0
  
  cigar_object = create_cigar(cigar_table = cigar_table, 
                              asv_table = asv_table,
                              asv_seqs = asv_seqs,
                              metadata = metadata)
  
  return(cigar_object)
  
}


## ampseq S4class and create_ampseq----

## cigar S4 class

setClass('ampseq', slots = c(
  gt = "ANY",
  asv_table = 'ANY',
  asv_seqs = 'ANY',
  metadata = "ANY",
  markers = "ANY",
  loci_performance = "ANY",
  pop_summary = "ANY",
  controls = "ANY",
  discarded_loci = "ANY",
  discarded_samples = "ANY",
  plots = "ANY"
))

## ampseq constructor
create_ampseq = function(gt = NULL,
                         asv_table = NULL,
                         asv_seqs = NULL,
                         metadata = NULL,
                         markers = NULL,
                         loci_performance = NULL,
                         pop_summary = NULL,
                         controls = NULL,
                         discarded_loci = NULL,
                         discarded_samples = NULL,
                         plots = NULL){
  obj = new('ampseq')
  obj@gt = gt
  obj@asv_table = asv_table
  obj@asv_seqs = asv_seqs
  obj@metadata = metadata
  obj@markers = markers
  obj@loci_performance = loci_performance
  obj@pop_summary = pop_summary
  obj@controls = controls
  obj@discarded_loci = discarded_loci
  obj@discarded_samples = discarded_samples
  obj@plots = plots
  
  return(obj)
}

## cigar2ampseq----

cigar2ampseq = function(cigar_object, min_abd = 1, min_ratio = .1, markers = NULL, markers_pattern = NULL, remove_controls = F){
  
  cigar_table = cigar_object@cigar_table
  asv_table = cigar_object@asv_table
  asv_seqs = cigar_object@asv_seqs
  
  if(!is.null(markers)){
    ampseq_loci_vector = markers[["amplicon"]]
    
  }else if(!is.null(markers_pattern)){
    all_alleles = gsub(markers_pattern, '', colnames(cigar_table))
    all_alleles = paste0(paste(c(all_alleles[all_alleles != '.'], '\\.'), collapse = '|'), '$')
    ampseq_loci_vector = unique(gsub(all_alleles, '',colnames(cigar_table)))
    
  }else if(is.null(markers) & is.null(markers_pattern)){
    ampseq_loci_vector = unique(sapply(strsplit(colnames(cigar_table), ","), function(x) x[1]))
  }
  
  ampseq_loci_abd_table = matrix(NA, nrow = nrow(cigar_table), ncol = length(ampseq_loci_vector), dimnames = list(rownames(cigar_table), ampseq_loci_vector))
  
  for(sample in rownames(ampseq_loci_abd_table)){
    for(locus in colnames(ampseq_loci_abd_table)){
      alleles = cigar_table[sample, grepl(paste0("^",locus,',|\\.)'), colnames(cigar_table))]
      
      if(length(alleles) == 1){
        names(alleles) <- colnames(cigar_table)[grepl(paste0("^",locus,',|\\.)'), colnames(cigar_table))]
      }
      
      if(length(alleles[which(alleles > min_abd)]) == 1){
        
        if(!is.null(markers_pattern)){
          ampseq_loci_abd_table[sample, locus] = paste(gsub(locus,'',names(alleles[which(alleles > min_abd)]), ","), alleles[which(alleles > min_abd)], sep = ":")
        }else{
          ampseq_loci_abd_table[sample, locus] = paste(gsub(paste0(locus, '.'), '', names(alleles[which(alleles > min_abd)]), ","), alleles[which(alleles > min_abd)], sep = ":")
        }
        
      }else if(length(alleles[which(alleles > min_abd)]) > 1){
        allele_names = names(alleles)
        alleles = as.integer(alleles)
        names(alleles) = allele_names
        alleles = sort(alleles, decreasing = T)
        
        if(!is.null(markers_pattern)){
          ampseq_loci_abd_table[sample, locus] = gsub(locus, "", paste(paste(names(alleles[alleles/max(alleles) > min_ratio]), alleles[alleles/max(alleles) > min_ratio], sep = ":"), collapse = "_"))  
        }else{
          ampseq_loci_abd_table[sample, locus] = gsub(paste0(locus, "."), "", paste(paste(names(alleles[alleles/max(alleles) > min_ratio]), alleles[alleles/max(alleles) > min_ratio], sep = ":"), collapse = "_"))  
        }
        
      }
    }
    rm(list = c("sample", "locus", "alleles"))
  }
  
  rm(ampseq_loci_vector)
  
  # Removing controls
  if(remove_controls){
    controls_ampseq_loci_abd_table = ampseq_loci_abd_table[cigar_object@metadata[cigar_object@metadata[["typeofSamp"]] == "Controls",][["Sample_id"]],]
    controls_metadata = cigar_object@metadata[cigar_object@metadata[["typeofSamp"]] == "Controls",]
    
    ampseq_loci_abd_table = ampseq_loci_abd_table[cigar_object@metadata[cigar_object@metadata[["typeofSamp"]] == "Samples",][["Sample_id"]],]
    metadata = cigar_object@metadata[cigar_object@metadata[["typeofSamp"]] == "Samples",]
    
    ampseq_object = create_ampseq(gt = ampseq_loci_abd_table,
                                  asv_table = asv_table,
                                  asv_seqs = asv_seqs,
                                  metadata = metadata,
                                  controls = list(gt = controls_ampseq_loci_abd_table,
                                                  metadata = controls_metadata),
                                  markers = markers,
                                  loci_performance = NULL,
                                  pop_summary = NULL)
    
    return(ampseq_object)
    
  }else{
    
    metadata = cigar_object@metadata
    
    ampseq_object = create_ampseq(gt = ampseq_loci_abd_table,
                                  asv_table = asv_table,
                                  asv_seqs = asv_seqs,
                                  metadata = metadata,
                                  markers = markers,
                                  loci_performance = NULL,
                                  pop_summary = NULL)
    
    return(ampseq_object)
    
  }
  
}

## join_ampseq----

join_ampseq = function(ampseq_obj_list = NULL){
  
  gt = NULL
  asv_table = NULL
  asv_seqs = NULL
  metadata = NULL
  markers = NULL
  
  for(obj in 1:length(ampseq_obj_list)){
    
    obj = ampseq_obj_list[[obj]]
    
    if(is.null(gt)){ # for the first object
      
      gt = obj@gt
      asv_table = obj@asv_table
      asv_seqs = obj@asv_seqs
      metadata = obj@metadata
      markers = obj@markers
      
    }else{# for subsequent objects
      
      temp_gt1 = obj@gt
      temp_asv_table1 = obj@asv_table
      temp_asv_seqs1 = obj@asv_seqs
      temp_metadata1 = obj@metadata
      temp_markers1 = obj@markers
      
      # if there are shared asvs
      if(sum(as.character(temp_asv_seqs1) %in% as.character(asv_seqs)) > 0){
        
        # Detect common ASVs from the previous and the new cigar table and sum up read counts, sample counts and bimera detection
        
        # asvs in the new data set present in the previous data set
        temp_asv_table2 = temp_asv_table1[as.character(temp_asv_seqs1) %in% as.character(asv_seqs),]
        temp_asv_seqs2 = temp_asv_seqs1[as.character(temp_asv_seqs1) %in% as.character(asv_seqs)]
        rownames(temp_asv_table2) = as.character(temp_asv_seqs2)
        
        # asvs in the previous data set present in the new data set
        temp_asv_table3 = asv_table[as.character(asv_seqs) %in% as.character(temp_asv_seqs1),]
        temp_asv_seqs3 = asv_seqs[as.character(asv_seqs) %in% as.character(temp_asv_seqs1)]
        rownames(temp_asv_table3) = as.character(temp_asv_seqs3)
        
        # Update the unique identifier of the asv (hapid)
        temp_asv_table2 = temp_asv_table2[rownames(temp_asv_table3),]
        temp_asv_table2[['hapid']] = temp_asv_table3[['hapid']] 
        
        # Update the read counts, the sample counts, and the bimera status
        temp_asv_table3[['total_reads']] = temp_asv_table3[['total_reads']] + temp_asv_table2[['total_reads']]
        temp_asv_table3[['total_samples']] = temp_asv_table3[['total_samples']] + temp_asv_table2[['total_samples']]
        temp_asv_table3[['bimera']] = temp_asv_table3[['bimera']] | temp_asv_table2[['bimera']]
        
        # Impute the inconsistent cigar_strings in the cigar_table
        unconsitent_cigar_strings = data.frame(temp_asv_table2[temp_asv_table2$CIGAR != temp_asv_table3$CIGAR,], 
                                               CIGAR2 = temp_asv_table3[temp_asv_table3$CIGAR != temp_asv_table2$CIGAR,][['CIGAR']],
                                               CIGAR_masked2 = temp_asv_table3[temp_asv_table3$CIGAR != temp_asv_table2$CIGAR,][['CIGAR']])
        
        if(nrow(unconsitent_cigar_strings) > 0){
          for(pos in 1:nrow(unconsitent_cigar_strings)){
            
            # cigar_string_replacment = unconsitent_cigar_strings[pos, ][['CIGAR2']]
            # cigar_string_pattern = unconsitent_cigar_strings[pos, ][['CIGAR']]
            cigar_string_masked_replacment = unconsitent_cigar_strings[pos, ][['CIGAR_masked2']]
            cigar_string_masked_pattern = unconsitent_cigar_strings[pos, ][['CIGAR_masked']]
            
            mhap = unconsitent_cigar_strings[pos, ][['Amplicon']]
            
            for(samp in 1:nrow(temp_gt1)){
              
              cigar_string_masked_pattern_in_samp = str_extract(temp_gt1[samp, mhap], paste0('(^|_)',cigar_string_masked_pattern, ':'))
              
              if(!is.na(cigar_string_masked_pattern_in_samp)){
                
                cigar_string_masked_replacment_in_samp = gsub(cigar_string_masked_pattern,
                                                              cigar_string_masked_replacment,
                                                              cigar_string_masked_pattern_in_samp)
                
                temp_gt1[samp, mhap] = gsub(cigar_string_masked_pattern_in_samp,
                                            cigar_string_masked_replacment_in_samp, 
                                            temp_gt1[samp, mhap])
                
              }
              rm(samp)
            }
            
            rm(list = c('pos', 
                        'mhap', 
                        'cigar_string_masked_replacment', 
                        'cigar_string_masked_pattern', 
                        'cigar_string_masked_pattern_in_samp',
                        'cigar_string_masked_replacment_in_samp'))
            
          }
        }
        
        rm(unconsitent_cigar_strings)
        
        # ASVs in the previous data set not present in the new one
        temp_asv_table4 = asv_table[!(as.character(asv_seqs) %in% as.character(temp_asv_seqs1)),]
        
        asv_table = rbind(temp_asv_table4, temp_asv_table3)
        rownames(asv_table) = gsub('ASV','',asv_table$hapid)
        
        rm(list = c('temp_asv_seqs2', 'temp_asv_seqs3', 'temp_asv_table2', 'temp_asv_table3', 'temp_asv_table4'))
        
      }
      
      temp_asv_table1 = temp_asv_table1[!(as.character(temp_asv_seqs1) %in% as.character(asv_seqs)),]
      temp_asv_seqs1 = temp_asv_seqs1[!(as.character(temp_asv_seqs1) %in% as.character(asv_seqs))]
      
      names(temp_asv_seqs1) = paste0('ASV',(1 + length(asv_seqs)):(length(temp_asv_seqs1) + length(asv_seqs)))
      temp_asv_table1$hapid = names(temp_asv_seqs1)
      
      # merging asv_seqs
      asv_seqs = DNAStringSet(c(as.character(asv_seqs), as.character(temp_asv_seqs1)))
      
      # merging asv_table
      asv_table = rbind(asv_table, temp_asv_table1)
      
      # merging marker table
      unshared_attributes = names(temp_markers1)[!(names(temp_markers1) %in% names(markers))]
      unshared_attributes = c(unshared_attributes, names(markers)[!(names(markers) %in% names(temp_markers1))])
      
      unshared_loci = temp_markers1$amplicon[!(temp_markers1$amplicon %in% markers$amplicon)]
      unshared_loci = c(unshared_loci, markers$amplicon[!(markers$amplicon %in% temp_markers1$amplicon)])
      
      if((length(unshared_loci) != 0 &
          length(unshared_attributes) != 0) |
         (length(unshared_loci) == 0 &
          length(unshared_attributes) != 0)){
        
        shared_fields = names(markers)[names(markers) %in% names(temp_markers1)]
        markers = full_join(markers, temp_markers1, by = shared_fields)
        
        rm(shared_fields)
        
        markers %<>% arrange(chromosome, start)
        
      }else if(length(unshared_loci) != 0 &
               length(unshared_attributes) == 0){
        
        unshared_loci_in_tablen = markers$amplicon[!(markers$amplicon %in% temp_markers1$amplicon)]
        temp_markers2 = temp_markers1[temp_markers1$amplicon %in% unshared_loci_in_tablen, names(markers)]
        
        markers = rbind(markers, temp_markers2)
        markers %<>% arrange(chromosome, start)
        
        rm(temp_markers2)
        
      }
      
      # merging gt genotype tables
      
      temp_gt2 = matrix(NA,
                        ncol = nrow(markers),
                        nrow = nrow(gt),
                        dimnames = list(rownames(gt),
                                        markers$amplicon)
      )
      
      for(amplicon in colnames(gt)){
        
        temp_gt2[,amplicon] = gt[, amplicon]
        rm(amplicon)
      }
      
      temp_gt3 = matrix(NA,
                        ncol = nrow(markers),
                        nrow = nrow(temp_gt1),
                        dimnames = list(rownames(temp_gt1),
                                        markers$amplicon)
      )
      
      for(amplicon in colnames(temp_gt1)){
        
        temp_gt3[,amplicon] = temp_gt1[, amplicon]
        rm(amplicon)
        
      }
      
      gt = rbind(temp_gt2, temp_gt3)
      
      # merging_metadata
      
      shared_metadata_fields = names(temp_metadata1)[(names(temp_metadata1) %in% names(metadata))]
      metadata =full_join(metadata, temp_metadata1, by = shared_metadata_fields)
      
      rm(list = c('temp_gt1', 
                  'temp_gt2', 
                  'temp_gt3',
                  'temp_markers1',
                  'unshared_attributes',
                  'unshared_loci',
                  'shared_metadata_fields', 'temp_asv_seqs1', 'temp_asv_table1'))
      
    }
    
    rm(list = c('obj'))
    
  }
  
  ampseq_object = create_ampseq(gt = gt,
                                asv_table = asv_table,
                                asv_seqs = asv_seqs,
                                metadata = metadata,
                                markers = markers)
  
  return(ampseq_object)
  
}

## write_ampseq----

write_ampseq = function(ampseq_object, format = c('excel', 'csv', 'json'), name = 'wb.xlsx'){
  
  if(format == 'excel'){
    
    if(file.exists(name)){
      system(paste0('rm ', name))
    }
    
    excel_wb = loadWorkbook(name, create = T)
    
    for(temp_slot in c('gt', 
                       'metadata', 
                       'markers', 
                       'loci_performance', 
                       'asv_table', 
                       'asv_seqs', 
                       'discarded_loci',
                       'discarded_samples',
                       'controls')){
      
      if(temp_slot == 'gt'){
        
        temp_sheet = data.frame(Sample_id = rownames(slot(ampseq_object, temp_slot)),
                                as.data.frame(slot(ampseq_object, temp_slot)))
        
      }else if(temp_slot == 'asv_seqs'){
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_sheet = data.frame(asv_id = names(slot(ampseq_object, temp_slot)),
                                  asv_seq = as.character(slot(ampseq_object, temp_slot)))
          
        }else{
          temp_sheet = NULL
        }
        
      }else if(temp_slot == 'markers'){
        
        temp_sheet = as.data.frame(slot(ampseq_object, temp_slot))
        temp_sheet[is.infinite(temp_sheet[['distance']]),][['distance']] = NA
        
      }else if(temp_slot == 'discarded_loci'){
        
        temp_sheet = NULL
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_discarded_loci = slot(ampseq_object, temp_slot)
          
          # write gt
          
          temp_discarded_loci_gt = data.frame(Sample_id = rownames(temp_discarded_loci[['gt']]),
                                              as.data.frame(temp_discarded_loci[['gt']]))
          
          createSheet(excel_wb, name = 'discarded_loci_gt')
          
          writeWorksheet(excel_wb,
                         temp_discarded_loci_gt,
                         sheet = 'discarded_loci_gt',
                         header = T)
          
          # write markers
          
          temp_discarded_loci_markers = as.data.frame(temp_discarded_loci[['markers']])
          temp_discarded_loci_markers[is.infinite(temp_discarded_loci_markers[['distance']]),][['distance']] = NA
          
          createSheet(excel_wb, name = 'discarded_loci_markers')
          
          writeWorksheet(excel_wb,
                         temp_discarded_loci_markers,
                         sheet = 'discarded_loci_markers',
                         header = T)
          
          # write loci_performance
          
          temp_discarded_loci_loci_performance = as.data.frame(temp_discarded_loci[['loci_performance']])
          
          createSheet(excel_wb, name = 'discarded_loci_loci_performance')
          
          writeWorksheet(excel_wb,
                         temp_discarded_loci_loci_performance,
                         sheet = 'discarded_loci_loci_performance',
                         header = T)
          
        }
        
        
      }else if(temp_slot == 'discarded_samples'){
        
        temp_sheet = NULL
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_discarded_samples = slot(ampseq_object, temp_slot)
          
          # write gt
          
          temp_discarded_samples_gt = data.frame(Sample_id = rownames(temp_discarded_samples[['gt']]),
                                              as.data.frame(temp_discarded_samples[['gt']]))
          
          createSheet(excel_wb, name = 'discarded_samples_gt')
          
          writeWorksheet(excel_wb,
                         temp_discarded_samples_gt,
                         sheet = 'discarded_samples_gt',
                         header = T)
          
          # write metadata
          
          temp_discarded_samples_metadata = as.data.frame(temp_discarded_samples[['metadata']])
          
          createSheet(excel_wb, name = 'discarded_samples_metadata')
          
          writeWorksheet(excel_wb,
                         temp_discarded_samples_metadata,
                         sheet = 'discarded_samples_metadata',
                         header = T)
        }
        
      }else if(temp_slot == 'controls'){
        
        temp_sheet = NULL
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_controls = slot(ampseq_object, temp_slot)
          
          # write gt
          
          temp_controls_gt = data.frame(Sample_id = rownames(temp_controls[['gt']]),
                                                 as.data.frame(temp_controls[['gt']]))
          
          createSheet(excel_wb, name = 'controls_gt')
          
          writeWorksheet(excel_wb,
                         temp_controls_gt,
                         sheet = 'controls_gt',
                         header = T)
          
          # write metadata
          
          temp_controls_metadata = as.data.frame(temp_controls[['metadata']])
          
          createSheet(excel_wb, name = 'controls_metadata')
          
          writeWorksheet(excel_wb,
                         temp_controls_metadata,
                         sheet = 'controls_metadata',
                         header = T)
        }
        
      }else{
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          temp_sheet = as.data.frame(slot(ampseq_object, temp_slot))
        }else{
          temp_sheet = NULL
        }
        
        
      }
      
      if(!is.null(temp_sheet)){
        createSheet(excel_wb, name = temp_slot)
        
        writeWorksheet(excel_wb,
                       temp_sheet,
                       sheet = temp_slot,
                       header = T)
      }
      
      
    }
    
    saveWorkbook(excel_wb)
    
  }else if(format == 'csv'){
    
    if(file.exists(name)){
      system(paste0('rm -r ', name))
    }
    
    system(paste0('mkdir ', name))
    
    for(temp_slot in c('gt', 
                       'metadata', 
                       'markers', 
                       'loci_performance', 
                       'asv_table', 
                       'asv_seqs', 
                       'discarded_loci',
                       'discarded_samples',
                       'controls')){
      
      if(temp_slot == 'gt'){
        
        temp_sheet = data.frame(Sample_id = rownames(slot(ampseq_object, temp_slot)),
                                as.data.frame(slot(ampseq_object, temp_slot)))
        
      }else if(temp_slot == 'asv_seqs'){
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_sheet = data.frame(asv_id = names(slot(ampseq_object, temp_slot)),
                                  asv_seq = as.character(slot(ampseq_object, temp_slot)))
          
        }else{
          temp_sheet = NULL
        }
        
      }else if(temp_slot == 'markers'){
        
        temp_sheet = as.data.frame(slot(ampseq_object, temp_slot))
        temp_sheet[is.infinite(temp_sheet[['distance']]),][['distance']] = NA
        
      }else if(temp_slot == 'discarded_loci'){
        
        temp_sheet = NULL
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_discarded_loci = slot(ampseq_object, temp_slot)
          
          # write gt
          
          temp_discarded_loci_gt = data.frame(Sample_id = rownames(temp_discarded_loci[['gt']]),
                                              as.data.frame(temp_discarded_loci[['gt']]))
          
          write.csv(temp_discarded_loci_gt,
                    file.path(name, 'discarded_loci_gt.csv'), 
                    quote = F, row.names = F)
          
          # write markers
          
          temp_discarded_loci_markers = as.data.frame(temp_discarded_loci[['markers']])
          temp_discarded_loci_markers[is.infinite(temp_discarded_loci_markers[['distance']]),][['distance']] = NA
          
          write.csv(temp_discarded_loci_markers,
                    file.path(name, 'discarded_loci_markers.csv'), 
                    quote = F, row.names = F)
          
          # write loci_performance
          
          temp_discarded_loci_loci_performance = as.data.frame(temp_discarded_loci[['loci_performance']])
          
          write.csv(temp_discarded_loci_loci_performance,
                    file.path(name, 'discarded_loci_loci_performance.csv'), 
                    quote = F, row.names = F)
          
        }
        
        
      }else if(temp_slot == 'discarded_samples'){
        
        temp_sheet = NULL
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_discarded_samples = slot(ampseq_object, temp_slot)
          
          # write gt
          
          temp_discarded_samples_gt = data.frame(Sample_id = rownames(temp_discarded_samples[['gt']]),
                                                 as.data.frame(temp_discarded_samples[['gt']]))
          
          write.csv(temp_discarded_samples_gt,
                    file.path(name, 'discarded_samples_gt.csv'), 
                    quote = F, row.names = F)
          
          # write metadata
          
          temp_discarded_samples_metadata = as.data.frame(temp_discarded_samples[['metadata']])
          
          write.csv(temp_discarded_samples_metadata,
                    file.path(name, 'discarded_samples_metadata.csv'),
                    quote = F, row.names = F)
        }
        
      }else if(temp_slot == 'controls'){
        
        temp_sheet = NULL
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_controls = slot(ampseq_object, temp_slot)
          
          # write gt
          
          temp_controls_gt = data.frame(Sample_id = rownames(temp_controls[['gt']]),
                                                 as.data.frame(temp_controls[['gt']]))
          
          write.csv(temp_controls_gt,
                    file.path(name, 'controls_gt.csv'), 
                    quote = F, row.names = F)
          
          # write metadata
          
          temp_controls_metadata = as.data.frame(temp_controls[['metadata']])
          
          write.csv(temp_controls_metadata,
                    file.path(name, 'controls_metadata.csv'),
                    quote = F, row.names = F)
        }
        
      }else{
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          temp_sheet = as.data.frame(slot(ampseq_object, temp_slot))
        }else{
          temp_sheet = NULL
        }
        
        
      }
      
      if(!is.null(temp_sheet)){
        
        write.csv(temp_sheet, paste0(file.path(name, temp_slot), '.csv'), quote = F, row.names = F)
        
      }
      
      
    }
    
    
  }else if(format == 'json'){
    # In development
  }
  
}

## read_ampseq----

read_ampseq = function(file = NULL, format = 'excel'){
  
  ampseq_object = create_ampseq()
  
  if(format == 'excel'){
    
    temp_wb = loadWorkbook(file)
    
    for(sheet in getSheets(temp_wb)){
      if(sheet == 'gt'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = as.matrix(temp_sheet[,-1])
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(ampseq_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet == 'asv_seqs'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_names = temp_sheet[[1]]
        temp_sheet = DNAStringSet(temp_sheet[[2]])
        names(temp_sheet) = temp_sheet_names
        
        slot(ampseq_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet %in% c('metadata', 'loci_performance')){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(ampseq_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet == 'markers'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet[is.na(temp_sheet[['distance']]),][['distance']] = Inf
        slot(ampseq_object, sheet, check = TRUE) = temp_sheet
        
      }else if(grepl('discarded_loci', sheet)){
        
        if(sheet == 'discarded_loci_gt'){
          temp_sheet = readWorksheet(temp_wb, sheet = sheet)
          temp_sheet_rownames = temp_sheet[,1]
          temp_sheet = as.matrix(temp_sheet[,-1])
          rownames(temp_sheet) = temp_sheet_rownames
          
          ampseq_object@discarded_loci[['gt']] = temp_sheet
          
        }else if(sheet == 'discarded_loci_markers'){
          
          temp_sheet = readWorksheet(temp_wb, sheet = sheet)
          temp_sheet[is.na(temp_sheet[['distance']]),][['distance']] = Inf
          ampseq_object@discarded_loci[['markers']] = temp_sheet
          
        }else if(sheet == 'discarded_loci_loci_performance'){
          
          temp_sheet = readWorksheet(temp_wb, sheet = sheet)
          temp_sheet_rownames = temp_sheet[,1]
          rownames(temp_sheet) = temp_sheet_rownames
          ampseq_object@discarded_loci[['loci_performance']] = temp_sheet
          
        }
        
      }else if(grepl('discarded_samples', sheet)){
        
        if(sheet == 'discarded_samples_gt'){
          temp_sheet = readWorksheet(temp_wb, sheet = sheet)
          temp_sheet_rownames = temp_sheet[,1]
          temp_sheet = as.matrix(temp_sheet[,-1])
          rownames(temp_sheet) = temp_sheet_rownames
          
          ampseq_object@discarded_samples[['gt']] = temp_sheet
          
        }else if(sheet == 'discarded_samples_metadata'){
          
          temp_sheet = readWorksheet(temp_wb, sheet = sheet)
          temp_sheet_rownames = temp_sheet[,1]
          rownames(temp_sheet) = temp_sheet_rownames
          ampseq_object@discarded_samples[['metadata']] = temp_sheet
          
        }
        
      }else if(grepl('controls', sheet)){
        
        if(sheet == 'controls_gt'){
          temp_sheet = readWorksheet(temp_wb, sheet = sheet)
          temp_sheet_rownames = temp_sheet[,1]
          temp_sheet = as.matrix(temp_sheet[,-1])
          rownames(temp_sheet) = temp_sheet_rownames
          
          ampseq_object@controls[['gt']] = temp_sheet
          
        }else if(sheet == 'controls_metadata'){
          
          temp_sheet = readWorksheet(temp_wb, sheet = sheet)
          temp_sheet_rownames = temp_sheet[,1]
          rownames(temp_sheet) = temp_sheet_rownames
          ampseq_object@controls[['metadata']] = temp_sheet
          
        }
        
      }else{
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        slot(ampseq_object, sheet, check = TRUE) = temp_sheet
        
      }
    }
    
  }else if(format == 'csv'){
    
    for(sheet in list.files(file)){
      if(sheet == 'gt.csv'){
        
        temp_sheet = read.csv(file.path(file, sheet))
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = as.matrix(temp_sheet[,-1])
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(ampseq_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet == 'asv_seqs.csv'){
        
        temp_sheet = read.csv(file.path(file, sheet))
        temp_sheet_names = temp_sheet[[1]]
        temp_sheet = DNAStringSet(temp_sheet[[2]])
        names(temp_sheet) = temp_sheet_names
        
        slot(ampseq_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet %in% c('metadata.csv', 'loci_performance.csv')){
        
        temp_sheet = read.csv(file.path(file, sheet))
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(ampseq_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet == 'markers.csv'){
        
        temp_sheet = read.csv(file.path(file, sheet))
        temp_sheet[is.na(temp_sheet[['distance']]),][['distance']] = Inf
        slot(ampseq_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
        
      }else if(grepl('discarded_loci', sheet)){
        
        if(sheet == 'discarded_loci_gt.csv'){
          temp_sheet = read.csv(file.path(file, sheet))
          temp_sheet_rownames = temp_sheet[,1]
          temp_sheet = as.matrix(temp_sheet[,-1])
          rownames(temp_sheet) = temp_sheet_rownames
          
          ampseq_object@discarded_loci[['gt']] = temp_sheet
          
        }else if(sheet == 'discarded_loci_markers.csv'){
          
          temp_sheet = read.csv(file.path(file, sheet))
          temp_sheet[is.na(temp_sheet[['distance']]),][['distance']] = Inf
          ampseq_object@discarded_loci[['markers']] = temp_sheet
          
        }else if(sheet == 'discarded_loci_loci_performance.csv'){
          
          temp_sheet = read.csv(file.path(file, sheet))
          temp_sheet_rownames = temp_sheet[,1]
          rownames(temp_sheet) = temp_sheet_rownames
          ampseq_object@discarded_loci[['loci_performance']] = temp_sheet
          
        }
        
      }else if(grepl('discarded_samples', sheet)){
        
        if(sheet == 'discarded_samples_gt.csv'){
          temp_sheet = read.csv(file.path(file, sheet))
          temp_sheet_rownames = temp_sheet[,1]
          temp_sheet = as.matrix(temp_sheet[,-1])
          rownames(temp_sheet) = temp_sheet_rownames
          
          ampseq_object@discarded_samples[['gt']] = temp_sheet
          
        }else if(sheet == 'discarded_samples_metadata.csv'){
          
          temp_sheet = read.csv(file.path(file, sheet))
          temp_sheet_rownames = temp_sheet[,1]
          rownames(temp_sheet) = temp_sheet_rownames
          ampseq_object@discarded_samples[['metadata']] = temp_sheet
          
        }
        
      }else if(grepl('controls', sheet)){
        
        if(sheet == 'controls_gt.csv'){
          temp_sheet = read.csv(file.path(file, sheet))
          temp_sheet_rownames = temp_sheet[,1]
          temp_sheet = as.matrix(temp_sheet[,-1])
          rownames(temp_sheet) = temp_sheet_rownames
          
          ampseq_object@controls[['gt']] = temp_sheet
          
        }else if(sheet == 'controls_metadata.csv'){
          
          temp_sheet = read.csv(file.path(file, sheet))
          temp_sheet_rownames = temp_sheet[,1]
          rownames(temp_sheet) = temp_sheet_rownames
          ampseq_object@controls[['metadata']] = temp_sheet
          
        }
        
      }else{
        
        temp_sheet = read.csv(file.path(file, sheet))
        slot(ampseq_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
        
      }
    }
    
  }else if(format == 'json'){
    # In development
  }
  
  return(ampseq_object)
  
}




## ampseq2loci----

ampseq2loci = function(ampseq_object){
  
  ampseq_loci_abd_table = ampseq_object@gt
  
  ampseq_allele_freq= NULL
  
  for(locus in colnames(ampseq_loci_abd_table)){
    alleles = levels(as.factor(unlist(strsplit(gsub(":[0-9]+", "", ampseq_loci_abd_table[,locus]), "_"))))
    ampseq_allele_freq[[locus]] = sort(sapply(alleles, function(allele){
      sum(grepl(paste("(^|_)", allele, ":", sep = ""), ampseq_loci_abd_table[,locus]))/sum(!is.na(unlist(strsplit(ampseq_loci_abd_table[,locus], "_"))))
    }), decreasing = T)
  }
  
  freq_table = matrix(NA, nrow = length(ampseq_allele_freq), ncol = max(sapply(ampseq_allele_freq, length)),
                      dimnames = list(names(ampseq_allele_freq),
                                      paste("Allele",seq(1,max(sapply(ampseq_allele_freq, length)),1), sep = "_")))
  
  for (locus in names(ampseq_allele_freq)){
    for(allele in 1:length(ampseq_allele_freq[[locus]])){
      freq_table[locus, allele] = ampseq_allele_freq[[locus]][allele]
    }
  }
  
  freq_table[is.na(freq_table)] = 0
  
  loci_table = gsub(":[0-9]+","" ,ampseq_loci_abd_table)
  
  for(locus in colnames(loci_table)){
    for (sample in rownames(loci_table)) {
      loci_table[sample,locus] = paste(which(names(ampseq_allele_freq[[locus]]) %in% strsplit(loci_table[sample,locus], "_")[[1]]), collapse = "_")
    }
  }
  
  loci_table[loci_table == ""] = NA
  
  return(list(loci_table = loci_table,
              freq_table = freq_table,
              ampseq_alleles_per_locus  = ampseq_allele_freq,
              markers = ampseq_object@markers))
  
}


# Functions to filter or mask alleles, loci, or samples on the ampseq object----

# filter_samples----


filter_samples = function(ampseq_object, v){
  
  obj = ampseq_object
  
  obj@gt = ampseq_object@gt[v,]
  obj@metadata = ampseq_object@metadata[v,]
  
  return(obj)
}

# filter_loci----

filter_loci = function(ampseq_object, v){
  
  obj = ampseq_object
  
  if(is.logical(v)){
    obj@gt = ampseq_object@gt[, v]
    obj@markers = ampseq_object@markers[v, ]
    
    if(!is.null(obj@loci_performance)){
      obj@loci_performance = ampseq_object@loci_performance[v, ]
    }
    
  }else if(is.character(v)){
    obj@gt = ampseq_object@gt[, colnames(ampseq_object@gt) %in% v]
    obj@markers = ampseq_object@markers[ampseq_object@markers$amplicon %in% v, ]
    
    if(!is.null(obj@loci_performance)){
      obj@loci_performance = ampseq_object@loci_performance[ampseq_object@loci_performance$loci %in% v, ]
    }
    
  }
  
  obj@markers[["distance"]] = Inf
  
  for(chromosome in levels(as.factor(obj@markers[["chromosome"]]))){
    for(amplicon in 1:(nrow(obj@markers[obj@markers[["chromosome"]] == chromosome,])-1)){
      obj@markers[
        obj@markers[["chromosome"]] == chromosome,
      ][amplicon, "distance"] = 
        obj@markers[obj@markers[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - 
        obj@markers[obj@markers[["chromosome"]] == chromosome,][amplicon, "pos"]
    }
  }
  
  return(obj)
}

# Mask alternative alleles----
setGeneric("mask_alt_alleles", function(obj = NULL, ref_fasta = NULL, mask_formula = "dVSITES_ij > 0.3", homopolymer_length = 5) standardGeneric("mask_alt_alleles"))

setMethod("mask_alt_alleles", signature(obj = "ampseq"),
          
          function(obj = NULL,  ref_fasta = NULL, mask_formula = "dVSITES_ij > 0.3", homopolymer_length = 5){
            
            gt = obj@gt
            mhaps = obj@markers
            asv_table = obj@asv_table
            
            mhaps$homopolymer_regions = NA
            
            ref_sequences = readDNAStringSet(ref_fasta)
            
            homopolymer_pattern = '(A{length,}|T{length,}|G{length,}|C{length,})'
            
            homopolymer_pattern = gsub('length', homopolymer_length, homopolymer_pattern)
            
            for(mhap in mhaps$amplicon){
              
              homopolymers = unlist(str_extract_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern))
              homopolymers_location = str_locate_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern)
              
              mhaps[mhaps$amplicon == mhap, ][['homopolymer_regions']] = paste(paste(homopolymers,
                                                                                     paste(homopolymers_location[[1]][,'start'],
                                                                                           homopolymers_location[[1]][,'end'], sep = '-'), sep = ':'),
                                                                               collapse = ',')
              
            }
            
            mhaps %<>% mutate(homopolymer_regions = case_when(
              homopolymer_regions == '' ~ NA,
              homopolymer_regions != '' ~ homopolymer_regions
            ))
            
            alt = sapply(colnames(gt), function(mhap){
              alt = unique(unlist(strsplit(gsub(':\\d+', '',gt[,mhap]), '_')))
              
              alt = paste(alt[!is.na(alt) & alt != '.'], collapse = ',')
            })
            
            gt_masked = gt
            
            gt = gsub(':\\d+', '',gt)
            
            # Heterozygous positions
            HetPos = matrix(grepl('_', gt), ncol = ncol(gt), nrow = nrow(gt))
            
            # Check if formula is correct
            
            if(grepl("(h_ij|h_ijminor|p_ij|P_ij|H_ij|H_ijminor|nVSITES_ij|dVSITES_ij|nSNVs_ij|dSNVs_ij|nINDELs_ij|dINDELs_ij|flanking_INDEL|SNV_in_homopolymer|INDEL_in_homopolymer|bimera)(<|>|!|=)+", mask_formula)){
              stop("All mathematical and logical operators must be separated by blank spaces in mask_formula")
            }
            
            # modify mask_formula
            
            if(grepl("SNV_in_homopolymer ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "SNV_in_homopolymer (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "SNV_in_homopolymer (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("SNV_in_homopolymer ", "ASVs_attributes_table_temp[['SNV_in_homopolymer']] ", mask_formula)
              }else{
                stop("Filter SNV_in_homopolymer is been called but there are spelling issues in this part of the formula")
              }
            }
            
            
            if(grepl("INDEL_in_homopolymer ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "INDEL_in_homopolymer (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "INDEL_in_homopolymer (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("INDEL_in_homopolymer ", "ASVs_attributes_table_temp[['INDEL_in_homopolymer']] ", mask_formula)
              }else{
                stop("Filter INDEL_in_homopolymer is been called but there are spelling issues in this part of the formula")
              }
            }
            
            
            if(grepl("flanking_INDEL ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "flanking_INDEL (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "flanking_INDEL (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("flanking_INDEL ", "ASVs_attributes_table_temp[['flanking_INDEL']] ", mask_formula)
              }else{
                stop("Filter flanking_INDEL is been called but there are spelling issues in this part of the formula")
              }
            }
            
            
            if(grepl("bimera ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "bimera (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "bimera (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("bimera ", "ASVs_attributes_table_temp[['bimera']] ", mask_formula)
              }else{
                stop("Filter bimera is been called but there are spelling issues in this part of the formula")
              }
            }
            
            
            if(grepl("dINDELs_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "dINDELs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "dINDELs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("dINDELs_ij ", "ASVs_attributes_table_temp[['dINDELs_ij']] ", mask_formula)
              }else{
                stop("Filter dINDELs_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            
            if(grepl("nINDELs_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "nINDELs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "nINDELs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("nINDELs_ij ", "ASVs_attributes_table_temp[['nINDELs_ij']] ", mask_formula)
              }else{
                stop("Filter nINDELs_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            
            if(grepl("dSNVs_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "dSNVs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "dSNVs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("dSNVs_ij ", "ASVs_attributes_table_temp[['dSNVs_ij']] ", mask_formula)
              }else{
                stop("Filter dSNVs_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            
            if(grepl("nSNVs_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "nSNVs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "nSNVs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("nSNVs_ij ", "ASVs_attributes_table_temp[['nSNVs_ij']] ", mask_formula)
              }else{
                stop("Filter nSNVs_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            
            if(grepl("dVSITES_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "dVSITES_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "dVSITES_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("dVSITES_ij ", "ASVs_attributes_table_temp[['dVSITES_ij']] ", mask_formula)
              }else{
                stop("Filter dVSITES_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            
            if(grepl("nVSITES_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "nVSITES_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "nVSITES_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("nVSITES_ij ", "ASVs_attributes_table_temp[['nVSITES_ij']] ", mask_formula)
              }else{
                stop("Filter nVSITES_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            
            if(grepl("h_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "h_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "h_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("h_ij ", "ASVs_attributes_table_temp[['h_ij']] ", mask_formula)
              }else{
                stop("Filter h_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            if(grepl("h_ijminor ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "h_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "h_ijminor (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("h_ijminor ", "ASVs_attributes_table_temp[['h_ijminor']] ", mask_formula)
              }else{
                stop("Filter h_ijminor is been called but there are spelling issues in this part of the formula")
              }
            }
            
            if(grepl("p_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "p_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "p_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("p_ij ", "ASVs_attributes_table_temp[['p_ij']] ", mask_formula)
              }else{
                stop("Filter p_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            if(grepl("P_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "P_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "P_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("P_ij ", "ASVs_attributes_table_temp[['P_ij']] ", mask_formula)
              }else{
                stop("Filter P_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            if(grepl("H_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "H_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "H_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("H_ij ", "ASVs_attributes_table_temp[['H_ij']] ", mask_formula)
              }else{
                stop("Filter H_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            if(grepl("H_ijminor ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "H_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "H_ijminor (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("H_ijminor ", "ASVs_attributes_table_temp[['H_ijminor']] ", mask_formula)
              }else{
                stop("Filter H_ijminor is been called but there are spelling issues in this part of the formula")
              }
            }
            
            mask_formula_check = str_split(mask_formula, "&|\\|")[[1]]
            mask_formula_check  = mask_formula_check[!grepl("ASVs_attributes_table_temp", mask_formula_check)]
            
            
            if(length(mask_formula_check) > 0){
              for(wrong_filter in mask_formula_check){
                print(paste0("Spelling error with filter ", wrong_filter))
              }
              stop("Execution halted, revise mask_filter argument.\nPossible filters are:\nh_ij, h_ijminor, p_ij, P_ij, H_ij, H_ijminor, nVSITES_ij, dVSITES_ij, nSNVs_ij, dSNVs_ij, nINDELs_ij, dINDELs_ij, flanking_INDEL, SNV_in_homopolymer, INDEL_in_homopolymer, bimera")
            }
            
            for(mhap in 1:ncol(gt)){
              temp_gts = gt[,mhap] # genotypes observed in that site
              alleles = strsplit(alt[mhap], ',')[[1]] # alternative alleles observed in that site
              
              if(length(alleles) != 0){
                
                # Vector of presence or absence of each alternative allele
                h_ij = t(sapply(alleles,
                                function(allele){
                                  P_ij = grepl(allele, temp_gts)
                                  
                                  H_ijminor = grepl(paste0('_',allele), temp_gts)
                                  
                                  # Samples where alternative alleles are present and the site is heterozygous
                                  
                                  H_ij = (P_ij == 1 & HetPos[,mhap] == 1)    
                                  
                                  H_ijminor = (H_ijminor == 1 & HetPos[,mhap] == 1)
                                  
                                  # Oberved Heterozygosity in mhap j
                                  
                                  Hetob_j = sum(HetPos[,mhap])/nrow(HetPos)
                                  
                                  # Number and Density of variant sites in the allele i of the mhap j
                                  nVSITES_ij = length(str_extract_all(allele, '\\d+')[[1]])
                                  dVSITES_ij = nVSITES_ij/mhaps[mhap,][['length']]
                                  
                                  VSITES_ij = unlist(strsplit(allele, '\\d+'))[-1]
                                  
                                  # Number and Density of SNPs and INDELs in the allele i of the mhap j
                                  nSNVs_ij = sum(nchar(VSITES_ij) == 1)
                                  nINDELs_ij = sum(nchar(VSITES_ij) != 1)
                                  
                                  dSNVs_ij = nSNVs_ij/mhaps[mhap,][['length']]
                                  dINDELs_ij = nINDELs_ij/mhaps[mhap,][['length']]
                                  
                                  # SNV or INDELs in homopolymer regions
                                  
                                  positions = as.integer(unlist(strsplit(allele, '(\\.|[ATCG]+|(I|D)=[ATCG]+)')))
                                  polymorphisms = unlist(strsplit(gsub('^\\d+', '', allele), '\\d+'))
                                  type_polymorphism = sapply(polymorphisms, function(polymorphism){
                                    if(nchar(polymorphism) == 1){
                                      'SNV'}else if(nchar(polymorphism) > 1){
                                        'INDEL'
                                      }
                                  })
                                  
                                  # SNV_in_homopolymer
                                  
                                  if('SNV' %in% type_polymorphism & !is.na(mhaps[mhap,][['homopolymer_regions']])){
                                    
                                    snv_positions = positions[type_polymorphism == 'SNV']
                                    snv_polymorphisms = polymorphisms[type_polymorphism == 'SNV']
                                    
                                    start = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '\\d+-'))))
                                    end = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '-\\d+'))))
                                    
                                    masked_region = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos], end[pos], 1)
                                    })))
                                    
                                    SNV_in_homopolymer = sum(snv_positions %in% masked_region) > 0
                                    
                                  }else{
                                    SNV_in_homopolymer = FALSE
                                  }
                                  
                                  
                                  # SNV_in_homopolymer pattern and replacement
                                  
                                  if(SNV_in_homopolymer == 1){
                                    
                                    SNV_in_homopolymer_pattern = allele
                                    
                                    snv_positions_remined = snv_positions[!(snv_positions %in% masked_region)]
                                    snv_polymorphisms_remined = snv_polymorphisms[!(snv_positions %in% masked_region)]
                                    
                                    indel_postions_to_add = positions[type_polymorphism == 'INDEL']
                                    indel_polymorphisms_to_add = polymorphisms[type_polymorphism == 'INDEL']
                                    
                                    SNV_in_homopolymer_replacement = data.frame(positions = c(snv_positions_remined,
                                                                                              indel_postions_to_add),
                                                                                polymorphisms = c(snv_polymorphisms_remined,
                                                                                                  indel_polymorphisms_to_add)
                                    )
                                    
                                    SNV_in_homopolymer_replacement = SNV_in_homopolymer_replacement[order(SNV_in_homopolymer_replacement$positions),]
                                    
                                    
                                    SNV_in_homopolymer_replacement = gsub(' ',
                                                                          '',
                                                                          paste(apply(SNV_in_homopolymer_replacement, 
                                                                                      1, 
                                                                                      function(pos){paste(pos, collapse = '')}), 
                                                                                collapse = ''))
                                    
                                    if(SNV_in_homopolymer_replacement == ''){
                                      SNV_in_homopolymer_replacement = '.'
                                    }
                                    
                                  }else{
                                    
                                    SNV_in_homopolymer_pattern = NA
                                    SNV_in_homopolymer_replacement = NA
                                  }
                                  
                                  
                                  # INDEL_in_homopolymer
                                  
                                  if('INDEL' %in% type_polymorphism & !is.na(mhaps[mhap,][['homopolymer_regions']])){
                                    
                                    indel_positions = positions[type_polymorphism == 'INDEL']
                                    indel_polymorphisms = polymorphisms[type_polymorphism == 'INDEL']
                                    
                                    start = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '\\d+-'))))
                                    end = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '-\\d+'))))
                                    
                                    masked_region = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos], end[pos], 1)
                                    })))
                                    
                                    INDEL_in_homopolymer = sum(indel_positions %in% masked_region) > 0
                                    
                                  }else{
                                    INDEL_in_homopolymer = FALSE
                                  }
                                  
                                  
                                  # INDEL_in_homopolymer pattern and replacement
                                  
                                  if(INDEL_in_homopolymer == 1){
                                    
                                    INDEL_in_homopolymer_pattern = allele
                                    
                                    indel_positions_remined = indel_positions[!(indel_positions %in% masked_region)]
                                    indel_polymorphisms_remined = indel_polymorphisms[!(indel_positions %in% masked_region)]
                                    
                                    snv_postions_to_add = positions[type_polymorphism == 'SNV']
                                    snv_polymorphisms_to_add = polymorphisms[type_polymorphism == 'SNV']
                                    
                                    INDEL_in_homopolymer_replacement = data.frame(positions = c(indel_positions_remined,
                                                                                                snv_postions_to_add),
                                                                                  polymorphisms = c(indel_polymorphisms_remined,
                                                                                                    snv_polymorphisms_to_add)
                                    )
                                    
                                    INDEL_in_homopolymer_replacement = INDEL_in_homopolymer_replacement[order(INDEL_in_homopolymer_replacement$positions),]
                                    
                                    
                                    INDEL_in_homopolymer_replacement = gsub(' ',
                                                                            '',
                                                                            paste(apply(INDEL_in_homopolymer_replacement, 
                                                                                        1, 
                                                                                        function(pos){paste(pos, collapse = '')}), 
                                                                                  collapse = ''))
                                    
                                    if(INDEL_in_homopolymer_replacement == ''){
                                      INDEL_in_homopolymer_replacement = '.'
                                    }
                                    
                                  }else{
                                    
                                    INDEL_in_homopolymer_pattern = NA
                                    INDEL_in_homopolymer_replacement = NA
                                  }
                                  
                                  
                                  # INDELs in flanking regions
                                  
                                  flanking_INDEL = as.integer(grepl(paste0(paste('(^1', mhaps[mhap,][['length']] + 1, sep = '|'),')(I|D)(=|\\.)[ATGC]+') ,allele))
                                  
                                  if(flanking_INDEL == 1){
                                    flanking_INDEL_pattern = str_extract_all(allele, paste0(paste('(^1', paste0('[ATGC]?', mhaps[mhap,][['length']] + 1), sep = '|'),')(I|D)(=|\\.)[ATGC]+\\d?'))[[1]]
                                    
                                    flanking_INDEL_replacement = gsub(paste0(paste('(^1', mhaps[mhap,][['length']] + 1, sep = '|'),')(I|D)(=|\\.)[ATGC]+'), '', flanking_INDEL_pattern)
                                    
                                    flanking_INDEL_replacement = ifelse(flanking_INDEL_replacement == '', '.', flanking_INDEL_replacement)
                                  }else{
                                    
                                    flanking_INDEL_pattern = NA
                                    
                                    flanking_INDEL_replacement = NA
                                  }
                                  
                                  # Bimeras
                                  
                                  bimera = unique(asv_table[asv_table[['Amplicon']] == colnames(gt)[mhap] &
                                                              !is.na(asv_table[['Amplicon']]) &
                                                              asv_table[['CIGAR_masked']] == allele
                                                            ,][['bimera']])
                                  
                                  if(length(bimera) >= 2){
                                    bimera = FALSE
                                  }else if(length(bimera) == 0){
                                    bimera = FALSE
                                  }
                                  
                                  c(sum(P_ij, na.rm = T), 
                                    sum(H_ij, na.rm = T), 
                                    sum(H_ijminor, na.rm = T),
                                    sum(H_ij, na.rm = T)/sum(P_ij, na.rm = T),
                                    ifelse(is.na(sum(H_ijminor, na.rm = T)/sum(H_ij, na.rm = T)), 
                                           0,
                                           sum(H_ijminor, na.rm = T)/sum(H_ij, na.rm = T)),
                                    
                                    Hetob_j,
                                    
                                    nVSITES_ij,
                                    dVSITES_ij,
                                    
                                    nSNVs_ij,
                                    dSNVs_ij,
                                    
                                    nINDELs_ij,
                                    dINDELs_ij,
                                    
                                    mhaps[mhap,][['homopolymer_regions']],
                                    
                                    SNV_in_homopolymer,
                                    SNV_in_homopolymer_pattern,
                                    SNV_in_homopolymer_replacement,
                                    
                                    INDEL_in_homopolymer,
                                    INDEL_in_homopolymer_pattern,
                                    INDEL_in_homopolymer_replacement,
                                    
                                    flanking_INDEL,
                                    paste(flanking_INDEL_pattern, collapse = '||'),
                                    paste(flanking_INDEL_replacement, collapse = '||'),
                                    bimera
                                    
                                  )
                                  
                                }, simplify = T))
                
                
                
                ASVs_attributes_table_temp = as.data.frame(cbind(alleles, h_ij))
                
                names(ASVs_attributes_table_temp) = c('Allele',
                                                      'P_ij',
                                                      'H_ij',
                                                      'H_ijminor',
                                                      'h_ij',
                                                      'h_ijminor',
                                                      'Hetob_j',
                                                      'nVSITES_ij',
                                                      'dVSITES_ij',
                                                      'nSNVs_ij',
                                                      'dSNVs_ij',
                                                      'nINDELs_ij',
                                                      'dINDELs_ij',
                                                      'homopolymer_regions',
                                                      'SNV_in_homopolymer',
                                                      'SNV_in_homopolymer_pattern',
                                                      'SNV_in_homopolymer_replacement',
                                                      'INDEL_in_homopolymer',
                                                      'INDEL_in_homopolymer_pattern',
                                                      'INDEL_in_homopolymer_replacement',
                                                      'flanking_INDEL',
                                                      'flanking_INDEL_pattern',
                                                      'flanking_INDEL_replacement',
                                                      'bimera')
                
                ASVs_attributes_table_temp[['P_ij']] = as.integer(ASVs_attributes_table_temp[['P_ij']])
                ASVs_attributes_table_temp[['H_ij']] = as.integer(ASVs_attributes_table_temp[['H_ij']])
                ASVs_attributes_table_temp[['H_ijminor']] = as.integer(ASVs_attributes_table_temp[['H_ijminor']])
                
                ASVs_attributes_table_temp[['h_ij']] = as.numeric(ASVs_attributes_table_temp[['h_ij']])
                ASVs_attributes_table_temp[['h_ijminor']] = as.numeric(ASVs_attributes_table_temp[['h_ijminor']])
                
                ASVs_attributes_table_temp[['Hetob_j']] = as.numeric(ASVs_attributes_table_temp[['Hetob_j']])
                
                
                ASVs_attributes_table_temp[['nVSITES_ij']] = as.integer(ASVs_attributes_table_temp[['nVSITES_ij']])
                ASVs_attributes_table_temp[['nSNVs_ij']] = as.integer(ASVs_attributes_table_temp[['nSNVs_ij']])
                ASVs_attributes_table_temp[['nINDELs_ij']] = as.integer(ASVs_attributes_table_temp[['nINDELs_ij']])
                
                ASVs_attributes_table_temp[['dVSITES_ij']] = as.numeric(ASVs_attributes_table_temp[['dVSITES_ij']])
                ASVs_attributes_table_temp[['dSNVs_ij']] = as.numeric(ASVs_attributes_table_temp[['dSNVs_ij']])
                
                ASVs_attributes_table_temp[['SNV_in_homopolymer']] = as.logical(ASVs_attributes_table_temp[['SNV_in_homopolymer']])
                ASVs_attributes_table_temp[['INDEL_in_homopolymer']] = as.logical(ASVs_attributes_table_temp[['INDEL_in_homopolymer']])
                
                ASVs_attributes_table_temp[['flanking_INDEL']] = as.logical(as.integer(ASVs_attributes_table_temp[['flanking_INDEL']]))
                
                ASVs_attributes_table_temp[['flanking_INDEL_pattern']] = unlist(sapply(ASVs_attributes_table_temp[['flanking_INDEL_pattern']],
                                                                                       function(pattern){
                                                                                         
                                                                                         if(pattern == 'NA'|is.na(pattern)){NA}else{
                                                                                           pattern
                                                                                         }
                                                                                         
                                                                                       }))
                
                ASVs_attributes_table_temp[['INDEL_in_homopolymer_replacement']] = unlist(sapply(ASVs_attributes_table_temp[['INDEL_in_homopolymer_replacement']],
                                                                                                 function(replacement){
                                                                                                   
                                                                                                   if(replacement == 'NA'|is.na(replacement)){NA}else{
                                                                                                     replacement
                                                                                                   }
                                                                                                   
                                                                                                 }))
                
                ASVs_attributes_table_temp[['bimera']] = as.logical(ASVs_attributes_table_temp[['bimera']])
                
                ASVs_attributes_table_temp[['p_ij']] = ASVs_attributes_table_temp[['P_ij']]/nrow(gt)
                
                ASVs_attributes_table_temp[['MHap']] = colnames(gt)[mhap]
                
                
                # Identify alleles below thresholds
                
                if(grepl('SNV_in_homopolymer',mask_formula)){
                  
                  replaced_alleles = ASVs_attributes_table_temp[
                    eval(parse(text = mask_formula)),][['SNV_in_homopolymer_pattern']]
                  
                  if(length(replaced_alleles) > 0){
                    
                    mask_formula2 = str_extract(mask_formula,
                                                "ASVs_attributes_table_temp\\[\\['SNV_in_homopolymer'\\]\\] (=|!)+ (TRUE|FALSE)")
                    
                    replaced_alleles = ASVs_attributes_table_temp[
                      eval(parse(text = mask_formula2)),][['SNV_in_homopolymer_pattern']]
                    
                    
                    replacement_alleles = ASVs_attributes_table_temp[
                      eval(parse(text = mask_formula2)),][['SNV_in_homopolymer_replacement']]
                    
                    
                    for(replaced_allele in 1:length(replaced_alleles)){
                      
                      temp_replaced_allele = replaced_alleles[replaced_allele]
                      temp_replacement_allele = replacement_alleles[replaced_allele]
                      
                      asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                  !is.na(asv_table[['Amplicon']]) &
                                  asv_table[['CIGAR_masked']] == temp_replaced_allele
                                ,][['CIGAR_masked']] = temp_replacement_allele
                      
                      for(sample in 1:nrow(gt_masked)){
                        
                        replaced_pattern = str_extract(gt_masked[sample,mhap], paste0('(^|_)',temp_replaced_allele, ':'))
                        
                        if(!is.na(replaced_pattern)){
                          
                          replacement_pattern = paste0(temp_replacement_allele, ':')
                          gt_masked[sample,mhap] = gsub(replaced_pattern, replacement_pattern, gt_masked[sample,mhap])
                        }
                        
                      }
                      
                    }
                    
                    for(sample in 1:nrow(gt_masked)){
                      
                      sample_alleles = gsub(':.+','',str_split(gt_masked[sample, mhap], '_')[[1]])
                      
                      if(sum(!is.na(sample_alleles)) > 1){
                        
                        sample_allele_readdepth = gsub('.+:','',str_split(gt_masked[sample, mhap], '_')[[1]])
                        
                        sample_alleles = data.frame(sample_allele = sample_alleles,
                                                    sample_allele_readdepth = as.integer(sample_allele_readdepth))
                        
                        sample_alleles %<>% summarise(sample_allele_readdepth = sum(sample_allele_readdepth),
                                                      .by = sample_allele)
                        
                        sample_alleles %<>% arrange(desc(sample_allele_readdepth))
                        
                        gt_masked[sample, mhap] = paste(paste(sample_alleles$sample_allele, sample_alleles$sample_allele_readdepth, sep = ':'), collapse = '_')
                        
                      }
                      
                    }
                    
                  }
                  
                  
                }else if(grepl('INDEL_in_homopolymer',mask_formula)){
                  
                  replaced_alleles = ASVs_attributes_table_temp[
                    eval(parse(text = mask_formula)),][['INDEL_in_homopolymer_pattern']]
                  
                  if(length(replaced_alleles) > 0){
                    
                    mask_formula2 = str_extract(mask_formula,
                                                "ASVs_attributes_table_temp\\[\\['INDEL_in_homopolymer'\\]\\] (=|!)+ (TRUE|FALSE)")
                    
                    replaced_alleles = ASVs_attributes_table_temp[
                      eval(parse(text = mask_formula2)),][['INDEL_in_homopolymer_pattern']]
                    
                    
                    replacement_alleles = ASVs_attributes_table_temp[
                      eval(parse(text = mask_formula2)),][['INDEL_in_homopolymer_replacement']]
                    
                    
                    for(replaced_allele in 1:length(replaced_alleles)){
                      
                      temp_replaced_allele = replaced_alleles[replaced_allele]
                      temp_replacement_allele = replacement_alleles[replaced_allele]
                      
                      asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                  !is.na(asv_table[['Amplicon']]) &
                                  asv_table[['CIGAR_masked']] == temp_replaced_allele
                                ,][['CIGAR_masked']] = temp_replacement_allele
                      
                      for(sample in 1:nrow(gt_masked)){
                        
                        replaced_pattern = str_extract(gt_masked[sample,mhap], paste0('(^|_)',temp_replaced_allele, ':'))
                        
                        if(!is.na(replaced_pattern)){
                          
                          replacement_pattern = paste0(temp_replacement_allele, ':')
                          gt_masked[sample,mhap] = gsub(replaced_pattern, replacement_pattern, gt_masked[sample,mhap])
                        }
                        
                      }
                      
                    }
                    
                    for(sample in 1:nrow(gt_masked)){
                      
                      sample_alleles = gsub(':.+','',str_split(gt_masked[sample, mhap], '_')[[1]])
                      
                      if(sum(!is.na(sample_alleles)) > 1){
                        
                        sample_allele_readdepth = gsub('.+:','',str_split(gt_masked[sample, mhap], '_')[[1]])
                        
                        sample_alleles = data.frame(sample_allele = sample_alleles,
                                                    sample_allele_readdepth = as.integer(sample_allele_readdepth))
                        
                        sample_alleles %<>% summarise(sample_allele_readdepth = sum(sample_allele_readdepth),
                                                      .by = sample_allele)
                        
                        sample_alleles %<>% arrange(desc(sample_allele_readdepth))
                        
                        gt_masked[sample, mhap] = paste(paste(sample_alleles$sample_allele, sample_alleles$sample_allele_readdepth, sep = ':'), collapse = '_')
                        
                      }
                      
                    }
                    
                  }
                  
                  
                }else if(grepl('flanking_INDEL',mask_formula)){
                  
                  replaced_alleles = ASVs_attributes_table_temp[
                    eval(parse(text = mask_formula)),][['flanking_INDEL_pattern']]
                  
                  if(length(replaced_alleles) > 0){
                    
                    mask_formula2 = str_extract(mask_formula,
                                                "ASVs_attributes_table_temp\\[\\['flanking_INDEL'\\]\\] (=|!)+ (TRUE|FALSE)")
                    
                    replaced_alleles = ASVs_attributes_table_temp[
                      eval(parse(text = mask_formula2)),][['flanking_INDEL_pattern']]
                    
                    
                    replacement_alleles = ASVs_attributes_table_temp[
                      eval(parse(text = mask_formula2)),][['flanking_INDEL_replacement']]
                    
                    
                    replacement_alleles = replacement_alleles[!duplicated(replaced_alleles)]
                    replaced_alleles = replaced_alleles[!duplicated(replaced_alleles)]
                    
                    for(replaced_allele in 1:length(replaced_alleles)){
                      
                      temp_replaced_allele = unlist(strsplit(replaced_alleles[replaced_allele], '\\|\\|'))
                      temp_replacement_allele = unlist(strsplit(replacement_alleles[replaced_allele], '\\|\\|'))
                      
                      asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                  !is.na(asv_table[['Amplicon']]) &
                                  grepl(temp_replaced_allele, asv_table[['CIGAR_masked']])
                                ,][['CIGAR_masked']] = gsub(temp_replaced_allele,
                                                            temp_replacement_allele, 
                                                            asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                                                        !is.na(asv_table[['Amplicon']]) &
                                                                        grepl(temp_replaced_allele, asv_table[['CIGAR_masked']])
                                                                      ,][['CIGAR_masked']])
                      
                      for(i in 1:length(temp_replaced_allele)){
                        
                        for(sample in 1:nrow(gt_masked)){
                          
                          replaced_pattern = temp_replaced_allele[i]
                          
                          if(grepl('[ATGC]', substr(replaced_pattern, nchar(replaced_pattern),nchar(replaced_pattern)))){
                            
                            replaced_pattern = paste0(replaced_pattern, ':')
                            replacement_pattern = paste0(temp_replacement_allele[i], ':')
                            
                          }else{
                            
                            replacement_pattern = temp_replacement_allele[i]
                            
                          }
                          
                          # Mask alleles below threshold
                          gt_masked[sample,mhap] = gsub(replaced_pattern, replacement_pattern, gt_masked[sample,mhap])
                          
                        }
                        
                      }
                      
                    }
                    
                    for(sample in 1:nrow(gt_masked)){
                      
                      sample_alleles = gsub(':.+','',str_split(gt_masked[sample, mhap], '_')[[1]])
                      
                      if(sum(!is.na(sample_alleles)) > 1){
                        
                        sample_allele_readdepth = gsub('.+:','',str_split(gt_masked[sample, mhap], '_')[[1]])
                        
                        sample_alleles = data.frame(sample_allele = sample_alleles,
                                                    sample_allele_readdepth = as.integer(sample_allele_readdepth))
                        
                        sample_alleles %<>% summarise(sample_allele_readdepth = sum(sample_allele_readdepth),
                                                      .by = sample_allele)
                        
                        sample_alleles %<>% arrange(desc(sample_allele_readdepth))
                        
                        gt_masked[sample, mhap] = paste(paste(sample_alleles$sample_allele, sample_alleles$sample_allele_readdepth, sep = ':'), collapse = '_')
                        
                      }
                      
                    }
                    
                  }
                  
                  
                }else{
                  
                  removed_alleles = ASVs_attributes_table_temp[
                    eval(parse(text = mask_formula)),][['Allele']]
                  
                  if(length(removed_alleles) > 0){
                    removed_pattern =paste('_?(', paste(removed_alleles, collapse = '|'), '):\\d+_?', sep = '')
                    # Mask alleles below threshold
                    gt_masked[,mhap] = gsub(removed_pattern, '', gt_masked[,mhap])
                    
                    if(length(gt_masked[gt_masked[,mhap] == '' &
                                        !is.na(gt_masked[,mhap])
                                        , mhap]) > 0){
                      
                      gt_masked[gt_masked[,mhap] == '' &
                                  !is.na(gt_masked[,mhap])
                                , mhap] = NA
                      
                    }
                    
                    if(length(asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                        !is.na(asv_table[['Amplicon']]) &
                                        asv_table[['CIGAR_masked']] %in% removed_alleles
                                        ,][['CIGAR_masked']]) > 0){
                      
                      asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                  !is.na(asv_table[['Amplicon']]) &
                                  asv_table[['CIGAR_masked']] %in% removed_alleles
                                ,][['CIGAR_masked']] = NA
                      
                    }
                    
                    
                    
                  }
                  
                }
                
                
              }
              
            }
            
            obj@gt = gt_masked
            obj@asv_table = asv_table
            
            return(obj)
            
          }
)

# Sample, loci or alleles (ASV) performance metrics----

## sample_TotalReadDepth----

sample_ReadDepth = function(ampseq_object, stat = c('sum', 'mean', 'median', 'sd', 'iqr', 'q25', 'q75')){
  
  coverage = sapply(1:ncol(ampseq_object@gt), function(mhap){
    
    temp_mhap = strsplit(ampseq_object@gt[,mhap], '_')
    
    sapply(1:length(temp_mhap), function(sample){
      
      sum(as.numeric(gsub('[^*]+:','', temp_mhap[[sample]])))
      
    })
  })
  
  colnames(coverage) = colnames(ampseq_object@gt)
  rownames(coverage) = rownames(ampseq_object@gt)
  
  output = data.frame(Sample_id = rownames(coverage))
  
  if('sum' %in% stat){
    output[['Total_ReadDepth']] = rowSums(coverage, na.rm = T)
  }
  
  if('mean' %in% stat){
    output[['Mean_ReadDepth']] = rowMeans(coverage, na.rm = T)
  }
  
  if('sd' %in% stat){
    output[['SD_ReadDepth']] = apply(coverage, 1, function(x){sd(x, na.rm = T)})
  }
  
  if('median' %in% stat){
    output[['Median_ReadDepth']] = apply(coverage, 1, function(x){median(x, na.rm = T)})
  }
  
  if('iqr' %in% stat){
    output[['IQR_ReadDepth']] = apply(coverage, 1, function(x){IQR(x, na.rm = T)})
  }
  
  if(sum(grepl('^q\\d+$',stat)) > 0){
    
    q = gsub('q','',stat[grepl('^q\\d+$',stat)])
    
    for(q_i in q){
      
      output[[paste0('q', q_i, '_ReadDepth')]] = apply(coverage, 1, function(x){quantile(x, probs = as.numeric(paste0('.', q_i)), na.rm = T)})
    }
    
    
  }
  
  return(output)
  
}

## get_ReadDepth_coverage----

get_ReadDepth_coverage = function(ampseq_object, variable){
  coverage = sapply(1:ncol(ampseq_object@gt), function(mhap){
    
    temp_mhap = strsplit(ampseq_object@gt[,mhap], '_')
    
    sapply(1:length(temp_mhap), function(sample){
      
      sum(as.numeric(gsub('[^*]+:','', temp_mhap[[sample]])))
      
    })
  })
  
  if(sum(!is.na(coverage)) == 0){
    coverage[is.na(coverage)] = 0
  }
  
  colnames(coverage) = colnames(ampseq_object@gt)
  rownames(coverage) = rownames(ampseq_object@gt)
  
  coverage = as.data.frame(coverage)
  coverage$Sample_id = rownames(coverage)
  
  if(!is.null(variable)){
    
    coverage$var = ampseq_object@metadata[[variable]]
    
  }else{
    
    coverage$var = 'Population'
    
  }
  
  coverage %<>% filter(!is.na(var))
  
  coverage %<>% pivot_longer(cols = all_of(colnames(ampseq_object@gt)),
                             names_to = 'Markers',
                             values_to = 'Read_depth')
  
  sample_order = coverage %>%
    group_by(Sample_id) %>%
    summarise(Read_depth = sum(Read_depth, na.rm = T)) %>%
    arrange(Read_depth) %>%
    select(Sample_id) %>% unlist
  
  coverage$Markers = factor(coverage$Markers,
                            levels = ampseq_object@markers$amplicon)
  
  plot_read_depth_heatmap = coverage %>% 
    ggplot(aes(x = Markers, y = factor(Sample_id, levels = sample_order), fill = log10(Read_depth+1)))+
    geom_tile()+
    # scale_fill_gradient(low="white", high="red",
    #                     breaks = 1:ceiling(log10(max(coverage$Read_depth, na.rm = T))),
    #                     labels = 10^(1:ceiling(log10(max(coverage$Read_depth, na.rm = T)))))+
    facet_wrap(var~., scales = 'free_y', ncol = 1)+
    labs(y = 'Samples',
         fill = "Read depth")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 6),
          axis.text.y = element_text(size = 6))
  

  
  plot_read_depth_violin = coverage %>% 
    ggplot(aes(x = var, y = log10(Read_depth + 1), color = var))+
    geom_violin()+
    geom_jitter(alpha = .1, width = .2)+
    # scale_y_continuous(breaks = log10(10^(0:ceiling(log10(max(coverage$Read_depth, na.rm = T))))+1),
    #                    labels = 10^(0:ceiling(log10(max(coverage$Read_depth, na.rm = T)))))+
    labs(x = variable, y = 'Read depth')+
    theme_bw()+
    theme(legend.position = 'none')
  
  if(max(coverage$Read_depth, na.rm = T) > 0) {
    plot_read_depth_heatmap = plot_read_depth_heatmap + 
      scale_fill_gradient(low="white", high="red",
                          breaks = 1:ceiling(log10(max(coverage$Read_depth, na.rm = T))),
                          labels = 10^(1:ceiling(log10(max(coverage$Read_depth, na.rm = T)))))
    
    plot_read_depth_violin = plot_read_depth_violin + 
      scale_y_continuous(breaks = log10(10^(0:ceiling(log10(max(coverage$Read_depth, na.rm = T))))+1),
                         labels = 10^(0:ceiling(log10(max(coverage$Read_depth, na.rm = T)))))
  }
  
  
  coverage_by_sample = coverage %>% 
    summarise(Total_Read_depth = sum(Read_depth, na.rm = T), .by = c(var, Sample_id))
  
  plot_read_depth_violin_by_sample = coverage_by_sample %>% 
    ggplot(aes(x = var, y = log10(Total_Read_depth + 1), color = var))+
    geom_violin()+
    geom_jitter(alpha = .1, width = .2)+
    # scale_y_continuous(breaks = log10(10^(0:ceiling(log10(max(coverage_by_sample$Total_Read_depth, na.rm = T))))+1),
    #                    labels = 10^(0:ceiling(log10(max(coverage_by_sample$Total_Read_depth, na.rm = T)))))+
    labs(x = variable, y = 'Sample Read depth')+
    theme_bw()+
    theme(legend.position = 'none')
  
  if(max(coverage_by_sample$Total_Read_depth, na.rm = T ) > 0) {
    plot_read_depth_violin_by_sample = plot_read_depth_violin_by_sample + 
      scale_y_continuous(breaks = log10(10^(0:ceiling(log10(max(coverage_by_sample$Total_Read_depth, na.rm = T))))+1),
                         labels = 10^(0:ceiling(log10(max(coverage_by_sample$Total_Read_depth, na.rm = T)))))
  }
  
  
  return(list(plot_read_depth_heatmap = plot_read_depth_heatmap,
              plot_read_depth_violin = plot_read_depth_violin,
              plot_read_depth_violin_by_sample = plot_read_depth_violin_by_sample))
  
}


## locus_amplification_rate----

locus_amplification_rate = function(ampseq_object, threshold = .65, update_loci = TRUE, strata = NULL, based_on_strata = FALSE){
  
  # , chr_lengths = c(640851,
  #                   947102,
  #                   1067971,
  #                   1200490,
  #                   1343557,
  #                   1418242,
  #                   1445207,
  #                   1472805,
  #                   1541735,
  #                   1687656,
  #                   2038340,
  #                   2271494,
  #                   2925236,
  #                   3291936)
  
  ampseq_loci_abd_table = ampseq_object@gt
  
  loci_performance = data.frame(loci = colnames(ampseq_loci_abd_table),
                                loci_ampl_rate_Total = apply(ampseq_loci_abd_table, 2, function(x) 1 - sum(is.na(x))/length(x)))
  
  
  if(!is.null(strata)){
    
    ampseq_object@metadata[['Strata']] = ampseq_object@metadata[[strata]]
    
    ampseq_object@metadata %<>% mutate(Strata = case_when(
      !is.na(Strata) ~ Strata,
      is.na(Strata) ~ 'Missing data'
    ))
    
    for(pop in unique(ampseq_object@metadata[['Strata']])){
      
      temp_ampseq_object = filter_samples(ampseq_object,
                                          v = (ampseq_object@metadata[['Strata']] == pop)
      )
      
      if(sum((ampseq_object@metadata[['Strata']] == pop)) == 1){
        
        temp_loci_ampl_rate = as.integer(is.na(temp_ampseq_object@gt))
        
      }else if(sum((ampseq_object@metadata[['Strata']] == pop)) >1){
        
        temp_loci_ampl_rate = apply(temp_ampseq_object@gt,
                                    2,
                                    function(x){1 - sum(is.na(x))/length(x)})
        
      }
      
      
      
      loci_performance = cbind(loci_performance,
                               temp_loci_ampl_rate
      )
      
      colnames(loci_performance) = c(colnames(loci_performance)[-ncol(loci_performance)],
                                     paste0('loci_ampl_rate_', pop))
      
    }
    
    
  }
  
  loci_pass_thres = loci_performance %>% 
    pivot_longer(cols = starts_with('loci_ampl_rate'), names_to = 'Strata', values_to = 'loci_ampl_rate') %>%
    mutate(Strata = gsub('loci_ampl_rate_', '', Strata),
           Strata = factor(Strata, 
                           levels = c(sort(unique(ampseq_object@metadata[['Strata']])), 'Total'))) %>%
    summarise(pass_loci = sum(loci_ampl_rate > threshold), .by = Strata)
  
  
  all_loci_performance_plot = loci_performance %>% 
    pivot_longer(cols = starts_with('loci_ampl_rate'), names_to = 'Strata', values_to = 'loci_ampl_rate') %>%
    mutate(Strata = gsub('loci_ampl_rate_', '', Strata),
           Strata = factor(Strata, 
                           levels = c(sort(unique(ampseq_object@metadata[['Strata']])), 'Total')))%>%
    ggplot(aes(loci_ampl_rate)) +
    geom_histogram(binwidth = 0.05, alpha = .7, boundary = 1) +
    geom_text(data = loci_pass_thres, aes(x = threshold - .1, y = max(pass_loci)/2, label = paste0(pass_loci, ' loci > ', threshold)))+
    labs(x = "Amplification rate",
         y = "# of Loci")+
    facet_grid(Strata~.)+
    geom_vline(xintercept = threshold) +
    theme_bw()
  
  # amplification_rate_per_locus = ggplot()+
  #   geom_segment(data = data.frame(x = 0, y = 1:14, xend =chr_lengths, yend = 1:14),
  #                aes(x=x, y=y, xend=xend, yend=yend), alpha = .5)+
  #   geom_point(data = data.frame(Position = ampseq_object@markers$pos,
  #                                Amplification_rate = 
  #                                  loci_performance$loci_ampl_rate,
  #                                Chr = as.character(gsub("(Pf3D7_|_v3|PvP01_|_v1)", "", ampseq_object@markers$chromosome))),
  #              aes(x = Position, y = Chr, color = Amplification_rate),
  #              pch = "|", size = 5)+
  #   theme_bw()+
  #   labs(title = "Amplification rate by locus",
  #        x = "Chromosome position",
  #        y = "Chromosome",
  #        color = NULL)+
  #   theme(legend.position = c(.9,.4))+
  #   scale_y_continuous(breaks = 1:14)+
  #   scale_color_continuous(type = "viridis")
  
  if(update_loci){
    
    if(!is.null(strata) & based_on_strata){
      
      discarded_loci = loci_performance[apply(loci_performance[, grepl('loci_ampl_rate', colnames(loci_performance))], 
                                              1,
                                              function(x){
                                                sum(x >= threshold) == 0
                                              }),][["loci"]]
      keeped_loci = loci_performance[apply(loci_performance[, grepl('loci_ampl_rate', colnames(loci_performance))], 
                                           1,
                                           function(x){
                                             sum(x >= threshold) > 0
                                           }),][["loci"]]
      
      ampseq_loci_abd_table_discarded_loci =
        ampseq_loci_abd_table[, colnames(ampseq_loci_abd_table) %in% discarded_loci]
      
      ampseq_loci_abd_table = 
        ampseq_loci_abd_table[, colnames(ampseq_loci_abd_table) %in% keeped_loci]
      
      markers = ampseq_object@markers
      
      discarded_markers = markers[markers[['amplicon']] %in% discarded_loci,]
      markers = markers[markers[['amplicon']] %in% keeped_loci,]
      
      markers[["distance"]] = Inf
      
      for(chromosome in levels(as.factor(markers[["chromosome"]]))){
        for(amplicon in 1:(nrow(markers[markers[["chromosome"]] == chromosome,])-1)){
          markers[
            markers[["chromosome"]] == chromosome,
          ][amplicon, "distance"] = 
            markers[markers[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - 
            markers[markers[["chromosome"]] == chromosome,][amplicon, "pos"]
        }
      }
      
      loci_performance_complete = loci_performance
      loci_performance = loci_performance[keeped_loci,]
      
      ampseq_object@gt = ampseq_loci_abd_table
      ampseq_object@markers = markers
      ampseq_object@loci_performance = loci_performance
      ampseq_object@discarded_loci = list(gt = ampseq_loci_abd_table_discarded_loci,
                                          loci_performance = loci_performance_complete,
                                          markers = discarded_markers)
      ampseq_object@plots[["all_loci_amplification_rate"]] = all_loci_performance_plot
      #ampseq_object@plots[["amplification_rate_per_locus"]] = amplification_rate_per_locus
      return(ampseq_object)
      
    }else if((is.null(strata) | !based_on_strata)){
      
      discarded_loci = loci_performance[loci_performance[["loci_ampl_rate_Total"]] <= threshold,][["loci"]]
      keeped_loci = loci_performance[loci_performance[["loci_ampl_rate_Total"]] > threshold,][["loci"]]
      
      ampseq_loci_abd_table_discarded_loci =
        ampseq_loci_abd_table[, colnames(ampseq_loci_abd_table) %in% discarded_loci]
      ampseq_loci_abd_table = 
        ampseq_loci_abd_table[, colnames(ampseq_loci_abd_table) %in% keeped_loci]
      
      markers = ampseq_object@markers
      
      discarded_markers = markers[markers[['amplicon']] %in% discarded_loci,]
      markers = markers[markers[['amplicon']] %in% keeped_loci,]
      
      markers[["distance"]] = Inf
      
      for(chromosome in levels(as.factor(markers[["chromosome"]]))){
        for(amplicon in 1:(nrow(markers[markers[["chromosome"]] == chromosome,])-1)){
          markers[
            markers[["chromosome"]] == chromosome,
          ][amplicon, "distance"] = 
            markers[markers[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - 
            markers[markers[["chromosome"]] == chromosome,][amplicon, "pos"]
        }
      }
      
      loci_performance_complete = loci_performance
      loci_performance = loci_performance[keeped_loci,]
      
      ampseq_object@gt = ampseq_loci_abd_table
      ampseq_object@markers = markers
      ampseq_object@loci_performance = loci_performance
      ampseq_object@discarded_loci = list(gt = ampseq_loci_abd_table_discarded_loci,
                                          loci_performance = loci_performance_complete,
                                          markers = discarded_markers)
      ampseq_object@plots[["all_loci_amplification_rate"]] = all_loci_performance_plot
      #ampseq_object@plots[["amplification_rate_per_locus"]] = amplification_rate_per_locus
      return(ampseq_object)
      
    }
    
    
  }else{
    
    return(all_loci_performance_plot)
    
  }
  
  
  
}

## sample_amplification_rate----

sample_amplification_rate = function(ampseq_object, threshold = .8, update_samples = TRUE, strata = NULL){
  
  metadata = ampseq_object@metadata
  ampseq_loci_abd_table = ampseq_object@gt
  loci_performance = ampseq_object@loci_performance
  
  metadata[["sample_ampl_rate"]] = apply(ampseq_loci_abd_table, 1, function(x) 1 - sum(is.na(x))/length(x))
  
  if(!is.null(strata)){
    
    
    sample_pass_thres = metadata %>% 
      summarise(pass_samples = sum(sample_ampl_rate >= threshold), .by = all_of(strata))
    
    all_samples_performance_plot = metadata %>% ggplot(aes(sample_ampl_rate))+
      geom_histogram(binwidth = 0.05, alpha = .7, boundary = 1)+
      geom_text(data = sample_pass_thres, aes(x = threshold - .1, 
                                              y = max(pass_samples)/3, 
                                              label = paste0(pass_samples, ' samples > ', threshold))) +
      labs(title = paste(sum(metadata[["sample_ampl_rate"]] > threshold), "Samples have an amplification rate >",threshold),
           x = "Amplification rate",
           y = "# of samples")+
      geom_vline(xintercept = threshold) +
      facet_grid(get(strata)~.)+
      theme_bw()
    
  }else{
    
    all_samples_performance_plot = metadata %>% ggplot(aes(sample_ampl_rate))+
      geom_histogram(binwidth = 0.05, alpha = .7, boundary = 1)+
      labs(title = paste(sum(metadata[["sample_ampl_rate"]] > threshold), "Samples have an amplification rate >",threshold),
           x = "Amplification rate",
           y = "# of samples")+
      geom_vline(xintercept = threshold) +
      theme_bw()
    
  }
  
  if(update_samples){
    
    ampseq_loci_abd_table_discarded_samples = ampseq_loci_abd_table[!(rownames(ampseq_loci_abd_table) %in% metadata[metadata[["sample_ampl_rate"]] > threshold ,][["Sample_id"]]),]
    ampseq_loci_abd_table = ampseq_loci_abd_table[rownames(ampseq_loci_abd_table) %in% metadata[metadata[["sample_ampl_rate"]] > threshold ,][["Sample_id"]],]
    
    metadata_complete = metadata
    metadata = metadata[metadata[["sample_ampl_rate"]] > threshold ,]
    
    
    # loci_performance[["loci_ampl_rate2"]] = apply(ampseq_loci_abd_table, 2, function(x) 1- sum(is.na(x))/length(x))
    # 
    # 
    # loci_performance_plot = loci_performance %>% ggplot(aes(loci_ampl_rate2))+
    #   geom_histogram(binwidth = 0.05, alpha = .7, boundary = 1)+
    #   labs(title = "Loci amplification rate",
    #        x = "Amplification rate",
    #        y = "# of Loci")+
    #   geom_vline(xintercept = threshold) +
    #   theme_bw()
    
    
    ampseq_object@gt = ampseq_loci_abd_table
    ampseq_object@metadata = metadata
    # ampseq_object@loci_performance = loci_performance
    
    ampseq_object@discarded_samples = list(gt = ampseq_loci_abd_table_discarded_samples,
                                           metadata = metadata_complete)
    
    # ampseq_object@plots[["loci_amplification_rate"]] = loci_performance_plot
    ampseq_object@plots[["samples_amplification_rate"]] = all_samples_performance_plot
    
    return(ampseq_object)
    
  }else{
    
    return(all_samples_performance_plot)
    
  }
  
}

## Get ASV's attributes----
setGeneric("get_ASVs_attributes", function(obj = NULL, ref_fasta = NULL, homopolymer_length = 5
) standardGeneric("get_ASVs_attributes"))

setMethod("get_ASVs_attributes", signature(obj = "ampseq"),
          
          function(obj = NULL, ref_fasta = NULL, homopolymer_length = 5){
            
            gt = obj@gt
            
            mhaps = obj@markers
            
            asv_table = obj@asv_table
            
            mhaps$homopolymer_regions = NA
            
            ref_sequences = readDNAStringSet(ref_fasta)
            
            homopolymer_pattern = '(A{length,}|T{length,}|G{length,}|C{length,})'
            
            homopolymer_pattern = gsub('length', homopolymer_length, homopolymer_pattern)
            
            for(mhap in mhaps$amplicon){
              
              homopolymers = unlist(str_extract_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern))
              homopolymers_location = str_locate_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern)
              
              mhaps[mhaps$amplicon == mhap, ][['homopolymer_regions']] = paste(paste(homopolymers,
                                                                                     paste(homopolymers_location[[1]][,'start'],
                                                                                           homopolymers_location[[1]][,'end'], sep = '-'), sep = ':'),
                                                                               collapse = ',')
              
            }
            
            mhaps %<>% mutate(homopolymer_regions = case_when(
              homopolymer_regions == '' ~ NA,
              homopolymer_regions != '' ~ homopolymer_regions
            ))
            
            alt = sapply(colnames(gt), function(mhap){
              alt = unique(unlist(strsplit(gsub(':\\d+', '',gt[,mhap]), '_')))
              
              alt = paste(alt[!is.na(alt) & alt != '.'], collapse = ',')
            })
            
            
            gt = gsub(':\\d+', '',gt)
            
            # Heterozygous positions
            HetPos = matrix(grepl('_', gt), ncol = ncol(gt), nrow = nrow(gt))
            
            ASVs_attributes_table = NULL
            
            for(mhap in 1:ncol(gt)){
              temp_gts = gt[,mhap] # genotypes observed in that site
              alleles = strsplit(alt[mhap], ',')[[1]] # alternative alleles observed in that site
              
              if(length(alleles) != 0){
                
                # Vector of presence or absence of each alternative allele
                h_ij = t(sapply(alleles,
                                function(allele){
                                  P_ij = grepl(allele, temp_gts)
                                  
                                  H_ijminor = grepl(paste0('_',allele), temp_gts)
                                  
                                  # Samples where alternative alleles are present and the site is heterozygous
                                  H_ij = (P_ij == 1 & HetPos[,mhap] == 1)    
                                  
                                  H_ijminor = (H_ijminor == 1 & HetPos[,mhap] == 1)
                                  
                                  Hetob_j = sum(HetPos[,mhap])/nrow(HetPos)
                                  
                                  # Number and Density of variant sites in the allele i of the mhap j
                                  nVSITES_ij = length(str_extract_all(allele, '\\d+')[[1]])
                                  dVSITES_ij = nVSITES_ij/mhaps[mhap,][['length']]
                                  
                                  VSITES_ij = unlist(strsplit(allele, '\\d+'))[-1]
                                  
                                  # Number and Density of SNPs and INDELs in the allele i of the mhap j
                                  nSNVs_ij = sum(nchar(VSITES_ij) == 1)
                                  nINDELs_ij = sum(nchar(VSITES_ij) != 1)
                                  
                                  dSNVs_ij = nSNVs_ij/mhaps[mhap,][['length']]
                                  dINDELs_ij = nINDELs_ij/mhaps[mhap,][['length']]
                                  
                                  # SNV or INDELs in homopolymer regions
                                  
                                  positions = as.integer(unlist(strsplit(allele, '(\\.|[ATCG]+|(I|D)=[ATCG]+)')))
                                  polymorphisms = unlist(strsplit(gsub('^\\d+', '', allele), '\\d+'))
                                  type_polymorphism = sapply(polymorphisms, function(polymorphism){
                                    if(nchar(polymorphism) == 1){
                                      'SNV'}else if(nchar(polymorphism) > 1){
                                        'INDEL'
                                      }
                                  })
                                  
                                  # SNV_in_homopolymer
                                  
                                  if('SNV' %in% type_polymorphism & !is.na(mhaps[mhap,][['homopolymer_regions']])){
                                    
                                    snv_positions = positions[type_polymorphism == 'SNV']
                                    snv_polymorphisms = polymorphisms[type_polymorphism == 'SNV']
                                    
                                    start = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '\\d+-'))))
                                    end = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '-\\d+'))))
                                    
                                    masked_region = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos], end[pos], 1)
                                    })))
                                    
                                    SNV_in_homopolymer = sum(snv_positions %in% masked_region) > 0
                                    
                                  }else{
                                    SNV_in_homopolymer = FALSE
                                  }
                                  
                                  
                                  # SNV_in_homopolymer pattern and replacement
                                  
                                  if(SNV_in_homopolymer == 1){
                                    
                                    SNV_in_homopolymer_pattern = allele
                                    
                                    snv_positions_remined = snv_positions[!(snv_positions %in% masked_region)]
                                    snv_polymorphisms_remined = snv_polymorphisms[!(snv_positions %in% masked_region)]
                                    
                                    indel_postions_to_add = positions[type_polymorphism == 'INDEL']
                                    indel_polymorphisms_to_add = polymorphisms[type_polymorphism == 'INDEL']
                                    
                                    SNV_in_homopolymer_replacement = data.frame(positions = c(snv_positions_remined,
                                                                                              indel_postions_to_add),
                                                                                polymorphisms = c(snv_polymorphisms_remined,
                                                                                                  indel_polymorphisms_to_add)
                                    )
                                    
                                    SNV_in_homopolymer_replacement = SNV_in_homopolymer_replacement[order(SNV_in_homopolymer_replacement$positions),]
                                    
                                    
                                    SNV_in_homopolymer_replacement = gsub(' ',
                                                                          '',
                                                                          paste(apply(SNV_in_homopolymer_replacement, 
                                                                                      1, 
                                                                                      function(pos){paste(pos, collapse = '')}), 
                                                                                collapse = ''))
                                    
                                    if(SNV_in_homopolymer_replacement == ''){
                                      SNV_in_homopolymer_replacement = '.'
                                    }
                                    
                                  }else{
                                    
                                    SNV_in_homopolymer_pattern = NA
                                    SNV_in_homopolymer_replacement = NA
                                  }
                                  
                                  
                                  # INDEL_in_homopolymer
                                  
                                  if('INDEL' %in% type_polymorphism & !is.na(mhaps[mhap,][['homopolymer_regions']])){
                                    
                                    indel_positions = positions[type_polymorphism == 'INDEL']
                                    indel_polymorphisms = polymorphisms[type_polymorphism == 'INDEL']
                                    
                                    start = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '\\d+-'))))
                                    end = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '-\\d+'))))
                                    
                                    masked_region = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos], end[pos], 1)
                                    })))
                                    
                                    INDEL_in_homopolymer = sum(indel_positions %in% masked_region) > 0
                                    
                                  }else{
                                    INDEL_in_homopolymer = FALSE
                                  }
                                  
                                  
                                  # INDEL_in_homopolymer pattern and replacement
                                  
                                  if(INDEL_in_homopolymer == 1){
                                    
                                    INDEL_in_homopolymer_pattern = allele
                                    
                                    indel_positions_remined = indel_positions[!(indel_positions %in% masked_region)]
                                    indel_polymorphisms_remined = indel_polymorphisms[!(indel_positions %in% masked_region)]
                                    
                                    snv_postions_to_add = positions[type_polymorphism == 'SNV']
                                    snv_polymorphisms_to_add = polymorphisms[type_polymorphism == 'SNV']
                                    
                                    INDEL_in_homopolymer_replacement = data.frame(positions = c(indel_positions_remined,
                                                                                                snv_postions_to_add),
                                                                                  polymorphisms = c(indel_polymorphisms_remined,
                                                                                                    snv_polymorphisms_to_add)
                                    )
                                    
                                    INDEL_in_homopolymer_replacement = INDEL_in_homopolymer_replacement[order(INDEL_in_homopolymer_replacement$positions),]
                                    
                                    
                                    INDEL_in_homopolymer_replacement = gsub(' ',
                                                                            '',
                                                                            paste(apply(INDEL_in_homopolymer_replacement, 
                                                                                        1, 
                                                                                        function(pos){paste(pos, collapse = '')}), 
                                                                                  collapse = ''))
                                    
                                    if(INDEL_in_homopolymer_replacement == ''){
                                      INDEL_in_homopolymer_replacement = '.'
                                    }
                                    
                                  }else{
                                    
                                    INDEL_in_homopolymer_pattern = NA
                                    INDEL_in_homopolymer_replacement = NA
                                  }
                                  
                                  
                                  # INDELs in flanking regions
                                  
                                  flanking_INDEL = as.integer(grepl(paste0(paste('(^1', mhaps[mhap,][['length']] + 1, sep = '|'),')(I|D)(=|\\.)[ATGC]+') ,allele))
                                  
                                  if(flanking_INDEL == 1){
                                    flanking_INDEL_pattern = str_extract_all(allele, paste0(paste('(^1', paste0('[ATGC]?', mhaps[mhap,][['length']] + 1), sep = '|'),')(I|D)(=|\\.)[ATGC]+\\d?'))[[1]]
                                    
                                    flanking_INDEL_replacement = gsub(paste0(paste('(^1', mhaps[mhap,][['length']] + 1, sep = '|'),')(I|D)(=|\\.)[ATGC]+'), '', flanking_INDEL_pattern)
                                    
                                    flanking_INDEL_replacement = ifelse(flanking_INDEL_replacement == '', '.', flanking_INDEL_replacement)
                                  }else{
                                    
                                    flanking_INDEL_pattern = NA
                                    
                                    flanking_INDEL_replacement = NA
                                  }
                                  
                                  # Bimeras
                                  
                                  bimera = unique(asv_table[asv_table[['Amplicon']] == colnames(gt)[mhap] &
                                                              !is.na(asv_table[['Amplicon']]) &
                                                              asv_table[['CIGAR_masked']] == allele
                                                            ,][['bimera']])
                                  
                                  if(length(bimera) >= 2){
                                    bimera = FALSE
                                  }else if(length(bimera) == 0){
                                    bimera = FALSE
                                  }
                                  
                                  
                                  c(sum(P_ij, na.rm = T), 
                                    sum(H_ij, na.rm = T), 
                                    sum(H_ijminor, na.rm = T),
                                    sum(H_ij, na.rm = T)/sum(P_ij, na.rm = T),
                                    ifelse(is.na(sum(H_ijminor, na.rm = T)/sum(H_ij, na.rm = T)), 
                                           0,
                                           sum(H_ijminor, na.rm = T)/sum(H_ij, na.rm = T)),
                                    
                                    Hetob_j,
                                    
                                    nVSITES_ij,
                                    dVSITES_ij,
                                    
                                    nSNVs_ij,
                                    dSNVs_ij,
                                    
                                    nINDELs_ij,
                                    dINDELs_ij,
                                    
                                    mhaps[mhap,][['homopolymer_regions']],
                                    
                                    SNV_in_homopolymer,
                                    SNV_in_homopolymer_pattern,
                                    SNV_in_homopolymer_replacement,
                                    
                                    INDEL_in_homopolymer,
                                    INDEL_in_homopolymer_pattern,
                                    INDEL_in_homopolymer_replacement,
                                    
                                    flanking_INDEL,
                                    paste(flanking_INDEL_pattern, collapse = '||'),
                                    paste(flanking_INDEL_replacement, collapse = '||'),
                                    
                                    bimera
                                    
                                  )
                                  
                                }, simplify = T))
                
                
                ASVs_attributes_table_temp = as.data.frame(cbind(alleles, h_ij))
                
                names(ASVs_attributes_table_temp) = c('Allele',
                                                      'P_ij',
                                                      'H_ij',
                                                      'H_ijminor',
                                                      'h_ij',
                                                      'h_ijminor',
                                                      'Hetob_j',
                                                      'nVSITES_ij',
                                                      'dVSITES_ij',
                                                      'nSNVs_ij',
                                                      'dSNVs_ij',
                                                      'nINDELs_ij',
                                                      'dINDELs_ij',
                                                      'homopolymer_regions',
                                                      'SNV_in_homopolymer',
                                                      'SNV_in_homopolymer_pattern',
                                                      'SNV_in_homopolymer_replacement',
                                                      'INDEL_in_homopolymer',
                                                      'INDEL_in_homopolymer_pattern',
                                                      'INDEL_in_homopolymer_replacement',
                                                      'flanking_INDEL',
                                                      'flanking_INDEL_pattern',
                                                      'flanking_INDEL_replacement',
                                                      'bimera')
                
                ASVs_attributes_table_temp[['P_ij']] = as.integer(ASVs_attributes_table_temp[['P_ij']])
                ASVs_attributes_table_temp[['H_ij']] = as.integer(ASVs_attributes_table_temp[['H_ij']])
                ASVs_attributes_table_temp[['H_ijminor']] = as.integer(ASVs_attributes_table_temp[['H_ijminor']])
                
                ASVs_attributes_table_temp[['h_ij']] = as.numeric(ASVs_attributes_table_temp[['h_ij']])
                ASVs_attributes_table_temp[['h_ijminor']] = as.numeric(ASVs_attributes_table_temp[['h_ijminor']])
                
                ASVs_attributes_table_temp[['Hetob_j']] = as.numeric(ASVs_attributes_table_temp[['Hetob_j']])
                
                ASVs_attributes_table_temp[['nVSITES_ij']] = as.integer(ASVs_attributes_table_temp[['nVSITES_ij']])
                ASVs_attributes_table_temp[['nSNVs_ij']] = as.integer(ASVs_attributes_table_temp[['nSNVs_ij']])
                ASVs_attributes_table_temp[['nINDELs_ij']] = as.integer(ASVs_attributes_table_temp[['nINDELs_ij']])
                
                ASVs_attributes_table_temp[['dVSITES_ij']] = as.numeric(ASVs_attributes_table_temp[['dVSITES_ij']])
                ASVs_attributes_table_temp[['dSNVs_ij']] = as.numeric(ASVs_attributes_table_temp[['dSNVs_ij']])
                
                ASVs_attributes_table_temp[['SNV_in_homopolymer']] = as.logical(ASVs_attributes_table_temp[['SNV_in_homopolymer']])
                ASVs_attributes_table_temp[['INDEL_in_homopolymer']] = as.logical(ASVs_attributes_table_temp[['INDEL_in_homopolymer']])
                
                ASVs_attributes_table_temp[['flanking_INDEL']] = as.logical(as.integer(ASVs_attributes_table_temp[['flanking_INDEL']]))
                
                ASVs_attributes_table_temp[['bimera']] = as.logical(ASVs_attributes_table_temp[['bimera']])
                
                ASVs_attributes_table_temp[['p_ij']] = ASVs_attributes_table_temp[['P_ij']]/nrow(gt)
                
                ASVs_attributes_table_temp[['MHap']] = colnames(gt)[mhap]
                
                
                ASVs_attributes_table = rbind(ASVs_attributes_table, ASVs_attributes_table_temp)
                
              }
              
              
            }
            
            ASVs_attributes_table %<>% mutate(SNV_in_homopolymer_pattern = case_when(
              SNV_in_homopolymer_pattern == 'NA' ~ NA,
              SNV_in_homopolymer_pattern != 'NA' ~ SNV_in_homopolymer_pattern),
              
              SNV_in_homopolymer_replacement = case_when(
                SNV_in_homopolymer_replacement == 'NA' ~ NA,
                SNV_in_homopolymer_replacement != 'NA' ~ SNV_in_homopolymer_replacement),
              
              INDEL_in_homopolymer_pattern = case_when(
                INDEL_in_homopolymer_pattern == 'NA' ~ NA,
                INDEL_in_homopolymer_pattern != 'NA' ~ INDEL_in_homopolymer_pattern),
              
              INDEL_in_homopolymer_replacement = case_when(
                INDEL_in_homopolymer_replacement == 'NA' ~ NA,
                INDEL_in_homopolymer_replacement != 'NA' ~ INDEL_in_homopolymer_replacement),
              
              flanking_INDEL_pattern = case_when(
                flanking_INDEL_pattern == 'NA' ~ NA,
                flanking_INDEL_pattern != 'NA' ~ flanking_INDEL_pattern),
              
              flanking_INDEL_replacement = case_when(
                flanking_INDEL_replacement == 'NA' ~ NA,
                flanking_INDEL_replacement != 'NA' ~ flanking_INDEL_replacement)
              
            )
            
            
            ASVs_attributes_table %<>% dplyr::select(MHap,
                                                     Allele,
                                                     P_ij,
                                                     H_ij,
                                                     H_ijminor,
                                                     p_ij,
                                                     h_ij,
                                                     h_ijminor,
                                                     Hetob_j,
                                                     nVSITES_ij,
                                                     dVSITES_ij,
                                                     nSNVs_ij,
                                                     dSNVs_ij,
                                                     nINDELs_ij,
                                                     dINDELs_ij,
                                                     homopolymer_regions,
                                                     SNV_in_homopolymer,
                                                     SNV_in_homopolymer_pattern,
                                                     SNV_in_homopolymer_replacement,
                                                     INDEL_in_homopolymer,
                                                     INDEL_in_homopolymer_pattern,
                                                     INDEL_in_homopolymer_replacement,
                                                     flanking_INDEL,
                                                     flanking_INDEL_pattern,
                                                     flanking_INDEL_replacement,
                                                     bimera)
            
            return(ASVs_attributes_table)
            
          }
)





# Functions for tertiary analysis----

## Functions for DRS and other variants of interest----
### haplotypes_respect_to_reference----

haplotypes_respect_to_reference = function(ampseq_object,
                                           gene_names = c('PfDHFR',
                                                          'PfMDR1',
                                                          'PfDHPS',
                                                          'PfKelch13',
                                                          'PF3D7_1447900'),
                                           gene_ids = c('PF3D7_0417200',
                                                        'PF3D7_0523000',
                                                        'PF3D7_0810800',
                                                        'PF3D7_1343700',
                                                        'PF3D7_1447900'),
                                           gff_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7.gff",
                                           fasta_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7_Genome.fasta",
                                           plot_haplo_freq = FALSE,
                                           variables = c('Sample_id', 'Population'),
                                           filters = NULL,
                                           na.var.rm = FALSE){
  library(ape)
  library(Biostrings)
  
  # Call 3D7 genome and gff---
  # Call 3D7 reference genome and its corresponding annotation in the gff file
  reference_gff = read.gff(gff_file)
  reference_genome = readDNAStringSet(fasta_file)
  
  # Call set of MHAPs---
  # Set of amplicon markers, This table contains the list mhap, their chromosome location in the 3D7 strain
  markers = ampseq_object@markers
  
  ## Filter Drug resistant markers---
  
  markers_of_interest = markers[grep(paste(gene_names, collapse = "|"),markers$amplicon),]
  
  ### Rename gene common name to gene_ID---
  
  markers_of_interest[['gene_ids']] = NA
  
  for(gene in 1:length(gene_names)){
    markers_of_interest[grepl(gene_names[gene], markers_of_interest[['amplicon']]),][['gene_ids']] = gene_ids[gene]
  }
  
  # Calculates the start and end position of each drugR marker on the CDS of each gene---
  
  ## Start and end position in 3D7 CDSs---
  markers_of_interest[['start_cds']] = NA
  markers_of_interest[['end_cds']] = NA
  markers_of_interest[['strand']] = NA
  markers_of_interest[['ref_length']] = NA
  
  for(amplicon in markers_of_interest$amplicon){ # for each drugR marker
    
    # Gene where the drugR marker is located
    gene = markers_of_interest[markers_of_interest$amplicon == amplicon,][['gene_ids']]
    
    # Filter the genes attributes of the CDS
    temp_gff = reference_gff[grepl(gene, reference_gff$attributes, ignore.case = TRUE) &
                               reference_gff$type == 'CDS',]
    
    # start counting from 0 taking into account the CDSs (exons) where the drugR marker is located
    start_cds = 0
    
    # check how many CDSs does the gene have
    if(nrow(temp_gff) > 1){ # if the gene has 2 or more CDSs then
      
      # Identify on which CDS the drugR marker is located
      cds = which(temp_gff$start < markers_of_interest[markers_of_interest$amplicon== amplicon,]$start &
                    temp_gff$end > markers_of_interest[markers_of_interest$amplicon== amplicon,]$start)
      
      if(cds > 1){ # if the drugR marker is located in CDS 2 or above then
        
        for(pre_cds in 1:(cds - 1)){ # sum the length of all previous CDSs
          start_cds = start_cds + temp_gff[pre_cds,][['end']] - temp_gff[pre_cds,][['start']] + 1
        }
        
        # sum the position on the located CDS
        start_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[cds,][['start']] + 1
        end_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$length - 1
        
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[cds,][['strand']])
        
        
      }else{ # if the drugR marker is located in CDS 1
        
        # Take the position on the located CDS
        
        start_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[cds,][['start']] + 1
        end_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$length - 1
        
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[cds,][['strand']])
        
      }
    }else{# if the gene has 1 CDS then
      
      # Take the position on the CDS 1
      start_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[['start']] + 1
      end_cds = start_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$length - 1
      
      markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds
      markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
      markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[['strand']])
      
    }
  }
  
  
  # Generate a reference sequence of each gene---
  # Genes located in the negative strand will not be transformed to their reverse complement yet
  
  # Reference sequences from 3D7---
  ref_seqs = NULL
  
  for(gene in unique(markers_of_interest$gene_ids)){
    
    temp_gff = reference_gff[grepl(gene, reference_gff$attributes, ignore.case = T)&
                               reference_gff$type == 'CDS',]
    
    temp_refseq = NULL
    
    if(nrow(temp_gff)>1){
      for(cds in 1:nrow(temp_gff)){
        temp_refseq = paste0(temp_refseq,
                             subseq(reference_genome[grep(temp_gff[cds,][['seqid']], names(reference_genome))],
                                    start = temp_gff[cds,][['start']],
                                    end = temp_gff[cds,][['end']]))
      }
      
      ref_seqs = c(ref_seqs, as.character(temp_refseq))
      
    }else{
      temp_refseq = subseq(reference_genome[grep(temp_gff[['seqid']], names(reference_genome))],
                           start = temp_gff[['start']],
                           end = temp_gff[['end']])
      
      ref_seqs = c(ref_seqs, as.character(temp_refseq))
    }
    
    markers_of_interest[markers_of_interest$gene_ids == gene,'ref_length'] = nchar(temp_refseq)
    
  }
  
  names(ref_seqs) = unique(markers_of_interest$gene_ids)
  ref_seqs = DNAStringSet(ref_seqs)
  
  
  # Convert cigar format to standard mutation nomenclature PMC1867422---
  
  ## Filter drugR markers---
  
  if(length(gene_ids) > 1){
    moi_loci_abd_table = ampseq_object@gt[,grep(paste(gene_names, collapse = "|"), colnames(ampseq_object@gt))]  
  }else{
    moi_loci_abd_table = matrix(ampseq_object@gt[,grep(paste(gene_names, collapse = "|"), colnames(ampseq_object@gt))],
                                ncol = 1,
                                dimnames = list(rownames(ampseq_object@gt),
                                                gene_names))
  }
  
  
  ## Remove read abundace---
  moi_loci_abd_table = gsub(":[0-9]+", "", moi_loci_abd_table)
  
  moi_loci_dna_table = moi_loci_abd_table
  moi_loci_aa_table = moi_loci_abd_table
  
  for(amplicon in colnames(moi_loci_abd_table)){ # For each amplicon in columns
    for(sample in 1:nrow(moi_loci_abd_table)){ # For each sample in rows
      
      locus = moi_loci_abd_table[sample, amplicon] # Get the genotype in the locus
      
      if(is.na(locus)){ # if the locus is NULL complete the cell with NA
        
        moi_loci_dna_table[sample, amplicon] = NA
        moi_loci_aa_table[sample, amplicon] = NA
        
      }else{
        
        clones = unlist(strsplit(locus, '_')) # get all different clones at that locus
        
        dna_clones = NULL
        aa_clones = NULL
        
        for(clone in clones){ # for each clone get alleles of all different SNPs
          
          alleles = paste0(unlist(strsplit(gsub('[ATCGDI//.]+$', '', clone), '[ATCGDI//.]+')),
                           unlist(strsplit(gsub('^[0-9]+', '', clone), '[0-9]+')))
          
          dna_alleles = NULL
          aa_alleles = NULL
          
          
          if(clone == "."){ # if allele equals to wild type
            
            codons = "."
            
          }else{
            
            codons = data.frame(alleles = gsub('[0-9]','',alleles), cds_position = sapply(alleles, function(allele){
              
              # position in the mhap
              mhap_position = as.integer(gsub('[ATCGDI//.]+', '', allele))
              
              # variant or nucleotide found
              mhap_variant = gsub('[0-9]+', '', allele)
              
              # calculate position in the CDS
              
              cds_position = markers_of_interest[markers_of_interest$amplicon == amplicon,'start_cds'] + mhap_position - 1
              
            }))
            
            
            # Identify nucleotide in the reference strain
            codons$ref_variant = sapply(1:nrow(codons), function(x){
              as.character(subseq(ref_seqs[which(names(ref_seqs) ==
                                                   markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
                                  start = codons[x, 'cds_position'],
                                  end = codons[x,'cds_position']))})
            
            # calculate the aminoacid position
            codons$aa_position = ceiling(codons$cds_position/3)
            codons$first_nucleotide = 3*codons$aa_position - 2
            codons$last_nucleotide = 3*codons$aa_position
            
          }
          
          if(length(codons) > 1){
            
            list_of_codons = unique(codons$aa_position)
            
          }else{
            list_of_codons = codons
            
          }
          
          
          for(codon in list_of_codons){ # for each allele
            
            if(codon == "."){ # if allele equals to wild type
              
              dna_alleles = c(dna_alleles, 'c.(=)')
              aa_alleles = c(aa_alleles, 'p.(=)')
              
            }else{ # if allele different than the wild type
              
              # get the reference amino acid variant
              ref_aa_variant = ifelse(markers_of_interest[markers_of_interest$amplicon == amplicon,'strand'] == "+",
                                      as.character(Biostrings::translate(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
                                                                                start = unique(codons[codons$aa_position == codon,][['first_nucleotide']]),
                                                                                end = unique(codons[codons$aa_position == codon,][['last_nucleotide']])))),
                                      as.character(Biostrings::translate(reverseComplement(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
                                                                                                  start = unique(codons[codons$aa_position == codon,][['first_nucleotide']]),
                                                                                                  end = unique(codons[codons$aa_position == codon,][['last_nucleotide']]))))))
              
              
              codon_positons = unique(codons[codons$aa_position == codon,][['first_nucleotide']]):unique(codons[codons$aa_position == codon,][['last_nucleotide']])
              
              positions_from_ref = codon_positons[!(unique(codons[codons$aa_position == codon,][['first_nucleotide']]):unique(codons[codons$aa_position == codon,][['last_nucleotide']]) %in% 
                                                      codons[codons$aa_position == codon,][['cds_position']])]
              
              positions_from_ref_df = NULL
              
              for(positions_from_ref1 in  positions_from_ref){
                positions_from_ref_df = rbind(positions_from_ref_df,
                                              data.frame(alleles = as.character(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_ids'])],
                                                                                       start = positions_from_ref1,
                                                                                       end = positions_from_ref1)),
                                                         cds_position = positions_from_ref1))
              }
              
              
              sample_codon = rbind(codons[codons$aa_position == codon,c('alleles', 'cds_position')],
                                   positions_from_ref_df
              )
              
              sample_codon = sample_codon[order(sample_codon$cds_position),]
              
              mhap_aa_variant = ifelse(markers_of_interest[markers_of_interest$amplicon == amplicon,'strand'] == "+",
                                       as.character(Biostrings::translate(DNAString(paste0(sample_codon$alleles, collapse = "")))),
                                       as.character(Biostrings::translate(reverseComplement(DNAString(paste0(sample_codon$alleles, collapse = ""))))))
              
              
              if(markers_of_interest[markers_of_interest$amplicon == amplicon,'strand'] == "+"){
                dna_alleles = c(dna_alleles, paste(paste0('c.',codons[codons$aa_position == codon,][['cds_position']],
                                                          codons[codons$aa_position == codon,][['ref_variant']],
                                                          '>',
                                                          codons[codons$aa_position == codon,][['alleles']]), collapse = " "))
                aa_alleles = c(aa_alleles, paste0(ref_aa_variant, codon, mhap_aa_variant))
              }else{
                dna_alleles = c(dna_alleles, paste(paste0('c.',
                                                          markers_of_interest[markers_of_interest$amplicon == amplicon,'ref_length'] - 
                                                            codons[codons$aa_position == codon,][['cds_position']] + 1,
                                                          reverseComplement(DNAStringSet(codons[codons$aa_position == codon,][['ref_variant']])), '>',
                                                          reverseComplement(DNAStringSet(codons[codons$aa_position == codon,][['alleles']]))), collapse = ' '))
                aa_alleles = c(aa_alleles, paste0(ref_aa_variant,
                                                  ceiling((markers_of_interest[markers_of_interest$amplicon == amplicon,'ref_length'] - 3*codon + 1)/3),
                                                  mhap_aa_variant))
                
              }
              
              
            }
            
            
          }
          
          dna_clones[[clone]] = paste(dna_alleles, collapse = " ")
          aa_clones[[clone]] = paste(aa_alleles, collapse = " ")
          
        }
        
        moi_loci_dna_table[sample, amplicon] = paste(dna_clones, collapse = " / ")
        moi_loci_aa_table[sample, amplicon] = paste(aa_clones, collapse = " / ")
        
        
      }
    }
  }
  
  
  if(plot_haplo_freq){
    
    
    # Empty table to fill cigar outputs
    aacigar_table = matrix(NA,
                           nrow = nrow(moi_loci_aa_table),
                           ncol = length(gene_names),
                           dimnames = list(rownames(moi_loci_aa_table),
                                           gene_names))
    
    
    for(gene in 1:length(gene_names)){
      
      if(length(gene_names) > 1){
        
        #gene_aa = moi_loci_aa_table[,grepl(gene_names[gene], colnames(moi_loci_aa_table))]
        
        gene_aa = matrix(moi_loci_aa_table[,grepl(gene_names[gene], colnames(moi_loci_aa_table))],
                         ncol = sum(grepl(gene_names[gene], colnames(moi_loci_aa_table))),
                         dimnames = list(rownames(moi_loci_aa_table),
                                         colnames(moi_loci_aa_table)[grepl(gene_names[gene], colnames(moi_loci_aa_table))]))
        
      }else{
        gene_aa = matrix(moi_loci_aa_table[,grepl(gene_names[gene], colnames(moi_loci_aa_table))],
                         ncol = 1,
                         dimnames = list(rownames(moi_loci_aa_table),
                                         gene_names))
      }
      
      # filter amplicons for the gene of interest
      gene_of_interest_info = markers_of_interest[markers_of_interest[['gene_ids']] == gene_ids[gene],]
      strand = unique(gene_of_interest_info[['strand']])
      
      # Pick the list of amplicons
      if(strand == '+'){
        
        # if the gene is located in the positive strand
        # sort the amplicons in ascending order
        gene_of_interest_info %<>% arrange(pos)
        amplicons = gene_of_interest_info$amplicon
        
      }else if(strand == '-'){
        
        # if the gene is located in the negative strand
        # sort the amplicons in descending order
        gene_of_interest_info %<>% arrange(desc(pos))
        amplicons = gene_of_interest_info$amplicon
        
      }
      
      # for each amplicon in the gene
      for(amplicon in amplicons){
        
        # get all observed polymorphic positions in the population sorted
        
        #positions = stringr::str_extract(unique(gene_aa[,amplicon]), '\\d+')
        positions = unlist(stringr::str_extract(gene_aa[,amplicon], '\\d+'))
        positions = positions[!is.na(positions)]
        if(sum(positions != "p.(=)") > 0){
          positions = unique(positions)
          positions = sort(positions)
        }else{
          positions = NULL
        }
        
        # for each sample
        for(sample in 1:nrow(gene_aa)){
          
          # get the full genotype (combination of haplotypes for each observed clone)
          sample_clones = gene_aa[sample, amplicon]
          
          
          if(!is.na(sample_clones)){# if the genotype is not NA
            
            # Split the observed haplotypes of each clone
            clones = unlist(strsplit(sample_clones, ' / ')) 
            
            # empty vector to write the final genotype
            clone_alleles = NULL
            
            # for each observed haplotype
            for(clone in clones){
              
              # if the haplotype is not equal to the reference
              if(clone != 'p.(=)'){
                
                # get the splited polymorphism of the clone in the sample
                sample_obs_alleles = unlist(strsplit(clone, ' '))
                
                # get the sorted polymorphism of the clone in the sample
                sample_alleles = NULL
                
                # for each polymorphic position in the population
                for(position in positions){
                  
                  # get the allele of the sample at the polymorphic position
                  temp_sample_allele = sample_obs_alleles[grepl(position, sample_obs_alleles)]
                  
                  # if the position was polymorphic
                  if(length(temp_sample_allele) > 0){
                    sample_alleles = c(sample_alleles, temp_sample_allele) 
                  }else{
                    
                    # if the position wasn't polymorphic
                    sample_allele = stringr::str_extract(unique(gene_aa[,amplicon]), '.\\d+.') #get all polymorphic sites
                    sample_allele = sample_allele[!is.na(sample_allele)]# remove empty data and haplotypes equals to the reference (p.(=))
                    sample_allele = unique(sample_allele) # remove duplicates
                    sample_allele = sample_allele[grepl(position, sample_allele)] # filter the allele at the desire position
                    sample_allele = gsub('\\d+.','',sample_allele) # get the reference allele at the polymorphic position
                    
                    sample_alleles = c(sample_alleles, paste0(sample_allele, position, sample_allele)) # define the sample allele equals to the reference
                    
                  }
                }
                
              }else{
                
                # if the full haplotype was equals to the reference
                sample_positions = NULL
                sample_alleles = NULL
                
                for(position in positions){
                  
                  sample_allele = stringr::str_extract(unique(gene_aa[,amplicon]), '.\\d+.')#get all polymorphic sites
                  sample_allele = sample_allele[!is.na(sample_allele)]# remove empty data and haplotypes equals to the reference (p.(=))
                  sample_allele = unique(sample_allele) # remove duplicates
                  sample_allele = sample_allele[grepl(position, sample_allele)] # filter the allele at the desire position
                  sample_allele = gsub('\\d+.','',sample_allele)# get the reference allele at the polymorphic position
                  
                  sample_alleles = c(sample_alleles, paste0(sample_allele, position, sample_allele)) # define the sample allele equals to the reference
                  
                }
              }
              
              clone_alleles = rbind(clone_alleles, sample_alleles) # bind haplotypes of each clone 
              
            }
            
            if(length(clones) > 1){ #if there are multiple clones
              
              # paste each clone allele with |
              sample_alleles = apply(clone_alleles, 2, function(x) ifelse(x[1] == x[2], x[1], paste(x[1], x[2], sep = '|')))
              
            }else{ # if there is only one clone
              
              sample_alleles = clone_alleles
              
            }
            
          }else{ # if the full haplotye was missing
            
            sample_positions = NULL
            sample_alleles = NULL
            
            for(position in positions){
              
              sample_allele = stringr::str_extract(unique(gene_aa[,amplicon]), '.\\d+.')# get all polymorphic sites
              sample_allele = sample_allele[!is.na(sample_allele)]# remove empty data and haplotypes equals to the reference (p.(=))
              sample_allele = unique(sample_allele) # remove duplicates
              sample_allele = sample_allele[grepl(position, sample_allele)]# filter the allele at the desire position
              sample_allele = gsub('\\d+.','',sample_allele)# get the reference allele at the polymorphic position
              
              sample_alleles = c(sample_alleles, paste0(sample_allele, position, "?"))# define the sample allele equals to "?"
              
            }
            
          }
          
          # if no data has been previously added to the cell in the aacigar_table
          if(is.na(aacigar_table[sample, gene])){
            
            aacigar_table[sample, gene] = paste(sample_alleles, collapse = ' ')
            
          }else{
            
            aacigar_table[sample, gene] = paste(aacigar_table[sample, gene], paste(sample_alleles, collapse = ' '), sep = ' ')
            
          } 
        }
      }
    }
    
    if(length(gene_names) > 1){
      mon_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) == 0),]
      poly_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) != 0),]  
    }else{
      mon_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) == 0),]
      mon_aacigar_table = matrix(mon_aacigar_table, ncol = 1,
                                 dimnames = list(
                                   names(mon_aacigar_table),
                                   gene_names
                                 ))
      poly_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) != 0),]
      
      if(!isEmpty(poly_aacigar_table)){
        poly_aacigar_table = matrix(poly_aacigar_table, ncol = 1,
                                    dimnames = list(
                                      names(mon_aacigar_table),
                                      gene_names
                                    ))
      }
    }
    
    
    if(isEmpty(poly_aacigar_table)){
      poly1 = gsub('\\|([a-z])', '', poly_aacigar_table, ignore.case = T)
      poly2 = gsub('([a-z]\\|)', '', poly_aacigar_table, ignore.case = T)
      
      extended_aacigar_table = rbind(mon_aacigar_table,
                                     poly1,
                                     poly2)
    }else{
      extended_aacigar_table = mon_aacigar_table
    }
    
    
    extended_aacigar_table = data.frame(samples = rownames(extended_aacigar_table), extended_aacigar_table)
    
    metadata = ampseq_object@metadata[,variables]
    
    
    if(length(variables)==2){
      names(metadata) = c('samples', 'var1')
    }else if(length(variables)==3){
      names(metadata) = c('samples', 'var1', 'var2')
    }else if(length(variables) > 3){print('This function allows up to 2 variables plus the sample_id variable')
    }else if(length(variables) <= 1){print('You have to privide at least one and no more than 2 variables plus you variable that identify the samples ids')}
    
    
    extended_aacigar_table = merge(extended_aacigar_table,
                                   metadata,
                                   by = 'samples',
                                   all.x = TRUE)
    
    extended_aacigar_table %<>% pivot_longer(cols = all_of(gene_names),
                                             names_to = 'gene_names',
                                             values_to = 'haplotype')
    
    extended_aacigar_table$gene_ids = NA
    
    for(gene in 1:length(gene_names)){
      extended_aacigar_table[extended_aacigar_table$gene_names == gene_names[gene],][['gene_ids']] = gene_ids[gene]
      
    }
    
    
    if(length(variables)==2){
      if(na.var.rm){
        extended_aacigar_table = extended_aacigar_table[!is.na(extended_aacigar_table$var1),]
      }else{
        extended_aacigar_table %<>% mutate(var1 = case_when(
          is.na(var1) ~ paste(variables[2], 'missing'),
          !is.na(var1) ~ var1))
      }
      
      haplotype_counts = extended_aacigar_table %>% group_by(gene_names, var1, haplotype)%>%
        summarise(count = n())
      
    }else if(length(variables)==3){
      if(na.var.rm){
        extended_aacigar_table = extended_aacigar_table[!is.na(extended_aacigar_table$var1)|
                                                          (!is.na(extended_aacigar_table$var2)|
                                                             grepl('NA',extended_aacigar_table$var2)),]
      }else{
        extended_aacigar_table %<>% mutate(var1 = case_when(
          is.na(var1) ~ paste(variables[2], 'missing'),
          !is.na(var1) ~ var1),
          var2 = case_when(
            is.na(var2)|grepl('NA',var2) ~ paste(variables[3], 'missing'),
            (!is.na(var2))&(!grepl('NA',var2)) ~ var2))
      }
      
      haplotype_counts = extended_aacigar_table %>% group_by(gene_names, var1, var2, haplotype)%>%
        summarise(count = n())
      
    }else if(length(variables) > 3){print('This function allows up to 2 variables plus the sample_id variable')
    }else if(length(variables) <= 1){print('You have to privide at least one and no more than 2 variables plus you variable that identify the samples ids')}
    
    
    if(!is.null(filters)){
      filters = strsplit(filters,';')
      for(temp_filter in 1:length(filters)){
        if(which(variables == filters[[temp_filter]][1]) == 2){
          haplotype_counts %<>% filter(var1 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
        }else if(which(variables == filters[[temp_filter]][1]) == 3){
          haplotype_counts %<>% filter(var2 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
        }
      }
    }
    
    
    haplotype_counts$freq = NA
    
    for(gene in levels(as.factor(haplotype_counts$gene_names))){
      for(Pop in levels(as.factor(haplotype_counts[haplotype_counts$gene_names == gene,][['var1']]))){
        
        if(length(variables) == 3){
          for(date in levels(as.factor(haplotype_counts[haplotype_counts$gene_names == gene&
                                                        haplotype_counts$var1 == Pop,][['var2']]))){
            
            haplotype_counts[haplotype_counts$gene_names == gene&
                               haplotype_counts$var1 == Pop&
                               haplotype_counts$var2 == date,][['freq']] = haplotype_counts[haplotype_counts$gene_names == gene&
                                                                                              haplotype_counts$var1 == Pop&
                                                                                              haplotype_counts$var2 == date,][['count']]/
              sum(haplotype_counts[haplotype_counts$gene_names == gene&
                                     haplotype_counts$var1 == Pop&
                                     haplotype_counts$var2 == date,][['count']])
            
          }
          
        }else{
          
          haplotype_counts[haplotype_counts$gene_names == gene&
                             haplotype_counts$var1 == Pop,][['freq']] = haplotype_counts[haplotype_counts$gene_names == gene&
                                                                                           haplotype_counts$var1 == Pop,][['count']]/
            sum(haplotype_counts[haplotype_counts$gene_names == gene&
                                   haplotype_counts$var1 == Pop,][['count']])
          
        }
        
        
      }
      
    }
    
    
    haplotype_counts %<>% mutate(gene_haplo = paste(gene_names, haplotype, sep = ": "))
    
    nhaplo = haplotype_counts %>% group_by(gene_names)%>% summarise(nhaplo = nlevels(as.factor(haplotype))) 
    
    if(max(nhaplo$nhaplo) <= 12){
      
      max_ncolors = brewer.pal(max(nhaplo$nhaplo), 'Set3')
      
    }else{
      
      max_ncolors = c(brewer.pal(12, 'Set3'), brewer.pal(max(nhaplo$nhaplo) - 12, 'Accent'))
      
    }
    
    colors = NULL
    
    for(n in 1:nrow(nhaplo)){
      colors = c(colors, max_ncolors[1:nhaplo[n,][['nhaplo']]])
    }
    
    
    if(length(variables) == 2){
      
      haplo_freq_plot = haplotype_counts%>%
        ggplot(aes(y = freq, x = var1, fill  = gene_haplo)) +
        geom_col()+
        facet_grid(.~ gene_names)+
        scale_fill_manual(values = colors)+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
        labs(y = 'Frequency in population',
             x = 'Population',
             fill = 'Gene: Haplotype')
      
      names(haplotype_counts) = c(names(haplotype_counts)[1], variables[2], names(haplotype_counts)[-1:-2])
      
    }else if(length(variables) == 3){
      
      haplo_freq_plot = haplotype_counts%>%
        ggplot(aes(y = freq, x = var2, fill  = gene_haplo)) +
        geom_col()+
        facet_grid(var1 ~ gene_names)+
        scale_fill_manual(values = colors)+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
        labs(y = 'Frequency in population',
             x = 'Date of Collection',
             fill = 'Gene: Haplotype')
      
      names(haplotype_counts) = c(names(haplotype_counts)[1], variables[2:3], names(haplotype_counts)[-1:-3])
      
    }
    
    
    haplotypes_respect_to_reference = list(markers_of_interest = markers_of_interest,
                                           loci_dna_table = moi_loci_dna_table,
                                           loci_aa_table = moi_loci_aa_table,
                                           ref_seqs = ref_seqs,
                                           aacigar_table = aacigar_table,
                                           haplotype_counts = haplotype_counts,
                                           haplo_freq_plot = haplo_freq_plot
    )
    
  }else{
    
    haplotypes_respect_to_reference = list(markers_of_interest = markers_of_interest,
                                           loci_dna_table = moi_loci_dna_table,
                                           loci_aa_table = moi_loci_aa_table,
                                           ref_seqs = ref_seqs)
    
    
  }
  
  
  return(haplotypes_respect_to_reference)
  
}

### drug_resistant_haplotypes----

drug_resistant_haplotypes = function(ampseq_object,
                                     reference_alleles = 'drugR_alleles.csv',
                                     gene_names = c('PfDHFR',
                                                    'PfMDR1',
                                                    'PfDHPS',
                                                    'PfKelch13',
                                                    'PF3D7_1447900'),
                                     gene_ids = c('PF3D7_0417200',
                                                  'PF3D7_0523000',
                                                  'PF3D7_0810800',
                                                  'PF3D7_1343700',
                                                  'PF3D7_1447900'),
                                     drugs = c('Artemisinine', 'Chloroquine', 'Pyrimethamine', 'Sulfadoxine', 'Lumefantrine', 'Mefloquine'),
                                     gff_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7.gff",
                                     fasta_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7_Genome.fasta",
                                     variables = c('samples', 'Population', 'quarter_of_collection'),
                                     Longitude = 'Longitude',
                                     Latitude = 'Latitude',
                                     filters = c('Population;Buenaventura,Quibdo,Guapi',
                                                 'quarter_of_collection;2020-Q4,2021-Q1,2021-Q2,2021-Q3,2021-Q4'),
                                     na.var.rm = FALSE,
                                     na.hap.rm = TRUE,
                                     hap_color_palette = 'random'){
  
  
  # Call reference alleles
  print("Uploading Resistant and Sensitive Alleles")
  drugR_reference_alleles = read.csv(reference_alleles)
  
  drugR_reference_alleles['reference'] = gsub('[0-9]+[A-Z]$', '', drugR_reference_alleles$Mutation)
  drugR_reference_alleles['mutant'] = gsub('^[A-Z][0-9]+', '', drugR_reference_alleles$Mutation)
  drugR_reference_alleles['position'] = gsub('^[A-Z]', '', gsub('[A-Z]$', '', drugR_reference_alleles$Mutation))
  
  print("Removing undesired categories based on var_filter")
  if(!is.null(filters)){
    filters = strsplit(filters,';')
    for(temp_filter in 1:length(filters)){
      
      ampseq_object = filter_samples(ampseq_object,
                                     ampseq_object@metadata[[filters[[temp_filter]][1]]] %in% strsplit(filters[[temp_filter]][2],',')[[1]])
    }
  }
  
  # Define haplotypes respect to a reference genome
  print("Defining haplotypes respect to a reference genome")
  haps_respect_to_ref = haplotypes_respect_to_reference(ampseq_object,
                                                        gene_names,
                                                        gene_ids,
                                                        gff_file,
                                                        fasta_file)
  
  loci_aa_table = haps_respect_to_ref[['loci_aa_table']]
  
  drug_markers = haps_respect_to_ref$markers_of_interest
  
  ref_seqs = haps_respect_to_ref$ref_seqs
  
  # Empty table to fill cigar outputs
  aacigar_table = matrix(NA,
                         nrow = nrow(loci_aa_table),
                         ncol = length(gene_ids),
                         dimnames = list(rownames(loci_aa_table),
                                         gene_ids))
  
  phenotype_table = matrix(NA,
                           nrow = nrow(loci_aa_table),
                           ncol = length(gene_ids),
                           dimnames = list(rownames(loci_aa_table),
                                           gene_ids))
  
  print("Defining aa haplotypes and phenotype respect to the presence of resistant alleles")
  for(gene in gene_ids){ # For each gene
    
    for(amplicon in drug_markers[drug_markers$gene_ids == gene, 'amplicon']){ # for each amplicon in the gene
      
      if(amplicon %in% colnames(loci_aa_table)){ # Apply only for amplicons present in loci table (some amplicons might be removed in previous steps)
        
        # Get all polymorphic positions in the population
        found_positions = loci_aa_table[,amplicon]
        
        found_positions = as.character(unlist(strsplit(as.character(unlist(strsplit(found_positions, " / "))), ' ')))
        
        
        found_positions = unique(as.character(gsub('^[A-Z]',
                                                   '',
                                                   gsub('[A-Z]$',
                                                        '',
                                                        found_positions))))
        
        found_positions = found_positions[!grepl('p',found_positions)]
        found_positions = found_positions[!is.na(found_positions)]
        
        # logic vector of amplified amplicon
        #tested_amplicon = drug_markers$amplicon %in% amplicon
        
        # For each sample compare their amino acid genotype respect to the table of reference
        for(sample in rownames(loci_aa_table)){ # For each sample
          # Identify all polymorphic positions in reference table
          
          positions = unique(c(found_positions, drugR_reference_alleles[drugR_reference_alleles$Gene_Id == gene, 'position']))
          
          # Identify tested positions
          
          if(drug_markers[drug_markers$amplicon == amplicon,'strand'][1] == '+'){
            
            tested_positions = seq(
              floor(
                drug_markers[drug_markers$amplicon == amplicon,
                             'start_cds']/3),
              ceiling(drug_markers[drug_markers$amplicon == amplicon,
                                   'end_cds']/3), 1)
            
          }else{
            
            tested_positions = seq(
              floor(
                (drug_markers[drug_markers$amplicon == amplicon,
                              'ref_length'] - 
                   drug_markers[drug_markers$amplicon == amplicon,
                                'end_cds'] + 1)/3),
              ceiling((drug_markers[drug_markers$amplicon == amplicon,
                                    'ref_length'] - 
                         drug_markers[drug_markers$amplicon == amplicon,
                                      'start_cds'] + 1)/3), 1)
            
          }
          
          
          # Keep polymorphic positions in the reference table that have been tested
          positions = positions[as.integer(positions) %in% tested_positions]
          
          # Select amplicons for defined sample and gene
          sample_clones = loci_aa_table[sample,
                                        amplicon]
          
          # if Genotype is not missing, Get aminoacid allele for each clone and each tested position in the sample 
          if(!is.na(sample_clones)){ 
            
            clones = unlist(strsplit(sample_clones, ' / '))
            
            
            clone_alleles = NULL
            
            for(clone in clones){
              
              if(clone != 'p.(=)'){
                
                sample_positions = gsub('^[A-Z]',
                                        '',
                                        gsub('[A-Z]$',
                                             '',
                                             unlist(strsplit(clone, ' '))))
                
                sample_alleles = gsub('^[A-Z][0-9]+',
                                      '',
                                      unlist(strsplit(clone, ' ')))
                
                remaining_positions = positions[!(positions %in% sample_positions)]
                
                for(position in remaining_positions){
                  
                  sample_positions = c(sample_positions, position)
                  
                  if(drug_markers[drug_markers$amplicon == amplicon,'strand'][1] == '+'){
                    
                    sample_alleles = c(sample_alleles,
                                       as.character(
                                         Biostrings::translate(
                                           subseq(
                                             ref_seqs[grep(gene,names(ref_seqs))],
                                             start = as.numeric(position)*3 - 2,
                                             end = as.numeric(position)*3))))
                    
                  }else{
                    
                    sample_alleles = c(sample_alleles,
                                       as.character(
                                         Biostrings::translate(
                                           reverseComplement(
                                             subseq(
                                               ref_seqs[grep(gene,names(ref_seqs))],
                                               start = drug_markers[drug_markers$amplicon == amplicon,
                                                                    'ref_length'] - as.numeric(position)*3 + 1,
                                               end = drug_markers[drug_markers$amplicon == amplicon,
                                                                  'ref_length'] - as.numeric(position)*3 + 3)))))
                    
                  }
                }
                
              }else{
                
                sample_positions = NULL
                sample_alleles = NULL
                
                for(position in positions){
                  
                  sample_positions = c(sample_positions, position)
                  
                  if(drug_markers[drug_markers$amplicon == amplicon,'strand'][1] == '+'){
                    
                    sample_alleles = c(sample_alleles,
                                       as.character(
                                         Biostrings::translate(
                                           subseq(
                                             ref_seqs[grep(gene,names(ref_seqs))],
                                             start = as.numeric(position)*3 - 2,
                                             end = as.numeric(position)*3))))
                    
                  }else{
                    
                    sample_alleles = c(sample_alleles,
                                       as.character(
                                         Biostrings::translate(
                                           reverseComplement(
                                             subseq(
                                               ref_seqs[grep(gene,names(ref_seqs))],
                                               start = drug_markers[drug_markers$amplicon == amplicon,
                                                                    'ref_length'] - as.numeric(position)*3 + 1,
                                               end = drug_markers[drug_markers$amplicon == amplicon,
                                                                  'ref_length'] - as.numeric(position)*3 + 3)))))
                    
                  }
                }
                
              }
              
              clone_alleles = rbind(clone_alleles, sample_alleles)
              
            }
            
            if(length(clones) > 1){
              sample_alleles = apply(clone_alleles, 2, function(x) ifelse(x[1] == x[2], x[1], paste(x[1], x[2], sep = '|')))
              
            }else{
              
              sample_alleles = clone_alleles
              
            }
            
            
          }else{ # else amino acid alleles are ?
            
            sample_positions = NULL
            sample_alleles = NULL
            
            for(position in positions){
              
              sample_positions = c(sample_positions, position)
              
              sample_alleles = c(sample_alleles, '?')
              
            }
            
          }
          
          # Get amino acid allele for each tested position in the reference table
          reference_alleles = NULL
          
          for(position in sample_positions){
            
            # If position is in reference table
            if(nrow(drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                            drugR_reference_alleles[['position']] == position,]) != 0){
              
              reference_alleles = c(reference_alleles, unique(drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                                                        drugR_reference_alleles[['position']] == position,'reference']))
            }else{
              
              if(drug_markers[drug_markers$amplicon == amplicon,'strand'][1] == '+'){
                
                reference_alleles = c(reference_alleles,
                                      as.character(
                                        Biostrings::translate(
                                          subseq(
                                            ref_seqs[grep(gene,names(ref_seqs))],
                                            start = as.numeric(position)*3 - 2,
                                            end = as.numeric(position)*3))))
                
              }else{
                
                reference_alleles = c(reference_alleles,
                                      as.character(
                                        Biostrings::translate(
                                          reverseComplement(
                                            subseq(
                                              ref_seqs[grep(gene,names(ref_seqs))],
                                              start = drug_markers[drug_markers$amplicon == amplicon,
                                                                   'ref_length'] - as.numeric(position)*3 + 1,
                                              end = drug_markers[drug_markers$amplicon == amplicon,
                                                                 'ref_length'] - as.numeric(position)*3 + 3)))))
                
              }
              
              
            }
            
          }
          
          
          # Write the haplotype for the amplicon (with one or more polymorphic sites)
          aacigar_haplotype = data.frame(position = as.integer(sample_positions),
                                         aacigar_haplotype = paste(reference_alleles, sample_positions, sample_alleles, sep = ""))
          
          aacigar_haplotype[['aacigar_haplotype']] = ifelse(aacigar_haplotype[['position']] %in% drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene,'position'],
                                                            aacigar_haplotype[['aacigar_haplotype']],
                                                            tolower(aacigar_haplotype[['aacigar_haplotype']]))
          
          # Get the phenotype for each position in the amplicon
          sample_phenotype = NULL
          
          for(position in sort(aacigar_haplotype[['position']])){
            
            sample_clones = gsub('[A-z][0-9]+', '', aacigar_haplotype[aacigar_haplotype[['position']] == position, 'aacigar_haplotype'])
            
            sample_clones = strsplit(sample_clones, '\\|')[[1]]
            
            if(length(sample_clones) > 1){# If sample is heterozygous for the position
              
              clone_phenotype = NULL
              
              for(sample_allele in sample_clones){ # For each allele in each clone in the position
                
                
                # if position in reference table
                if(nrow(drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                drugR_reference_alleles[['position']] == position,]) != 0){
                  
                  
                  # Get reference sensitive alleles
                  reference_sensitive_alleles = drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                                          drugR_reference_alleles[['position']] == position,'reference']
                  
                  # Get reference resistant alleles
                  reference_resistant_alleles = drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                                          drugR_reference_alleles[['position']] == position,'mutant']
                  
                  # If observed allele in reference sensitive alleles, the phenotype is sensitive
                  if(sample_allele %in% reference_sensitive_alleles){
                    
                    temp_clone_phenotype = 'Sensitive phenotype' 
                    
                  }else if(sample_allele %in% reference_resistant_alleles){
                    # If observed allele in reference resistant alleles, the phenotype is resistant
                    
                    phenotype = drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                          drugR_reference_alleles[['position']] == position  &
                                                          drugR_reference_alleles[['mutant']] == sample_allele,'Annotation']
                    
                    temp_clone_phenotype = phenotype
                    
                  }else if(sample_allele != '?'){
                    
                    temp_clone_phenotype = paste0(sample_allele, ' variant unreported for position ', position, ' in gene', gene)
                    
                  }
                  
                }else{# if position is not in reference table
                  
                  reference_allele = gsub('[0-9]+([A-z]|\\?|[A-z]\\|[A-z])', '', aacigar_haplotype[aacigar_haplotype[['position']] == position,'aacigar_haplotype'], ignore.case = T)
                  
                  if(reference_allele == sample_allele){
                    
                    temp_clone_phenotype = 'Sensitive phenotype'
                    
                  }else if(reference_allele != sample_allele){
                    
                    temp_clone_phenotype = paste0(reference_allele, position, sample_allele,
                                                  ' polymorphism in gene ',
                                                  gene,
                                                  ' respect to Reference Strain')
                    
                  }
                  
                }
                
                if(is.null(clone_phenotype)){
                  
                  clone_phenotype = temp_clone_phenotype
                  
                }else{
                  clone_phenotype = paste(clone_phenotype, temp_clone_phenotype, sep = "|")
                }
                
              }
              
              sample_phenotype = c(sample_phenotype, clone_phenotype)
              
            }else{# If sample is homozygous for the position
              
              sample_allele = sample_clones
              
              
              # if position in reference table
              if(nrow(drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                              drugR_reference_alleles[['position']] == position,]) != 0){
                
                
                reference_sensitive_alleles = drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                                        drugR_reference_alleles[['position']] == position,'reference']
                
                if(sample_allele %in% reference_sensitive_alleles){
                  
                  sample_phenotype = c(sample_phenotype, 'Sensitive phenotype')  
                  
                }else{
                  
                  reference_resistant_alleles = drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                                          drugR_reference_alleles[['position']] == position,'mutant']
                  
                  if(sum(strsplit(sample_allele, "\\|")[[1]] %in% reference_resistant_alleles) >= 1){
                    
                    if(grepl('|', sample_allele)){
                      
                      sample_allele = strsplit(sample_allele, "\\|")[[1]]
                      
                      sample_allele = sample_allele[sample_allele %in% reference_resistant_alleles]
                      
                    }
                    
                    phenotype = drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                                          drugR_reference_alleles[['position']] == position  &
                                                          drugR_reference_alleles[['mutant']] == sample_allele,'Annotation']
                    
                    sample_phenotype = c(sample_phenotype, phenotype)
                    
                  }else{
                    
                    if(sample_allele != '?'){
                      
                      sample_phenotype = c(sample_phenotype, paste0(sample_allele, ' variant unreported for position ', position, ' in gene', gene))
                      
                    }else{
                      
                      sample_phenotype = c(sample_phenotype, paste0(amplicon, ' amplicon did not amplify'))
                      
                    }
                    
                  }
                  
                }
                
              }else{# if position is not in reference table
                
                reference_allele = gsub('[0-9]+([A-z]|\\?)', '', aacigar_haplotype[aacigar_haplotype[['position']] == position,'aacigar_haplotype'], ignore.case = T)
                
                if(reference_allele == sample_allele){
                  
                  sample_phenotype = c(sample_phenotype, 'Sensitive phenotype')  
                  
                }else if(reference_allele != sample_allele & !grepl('\\?',sample_allele)){
                  
                  sample_phenotype = c(sample_phenotype, paste0(aacigar_haplotype[aacigar_haplotype[['position']] == position,'aacigar_haplotype'], ' polymorphism in gene ', gene,' respect to Reference Strain'))  
                  
                }else if(grepl('\\?',sample_allele)){
                  
                  sample_phenotype = c(sample_phenotype, paste0(amplicon, ' amplicon did not amplify'))
                  
                }
                
              }
            }
            
            
          }
          
          sample_phenotype = paste(sample_phenotype, collapse = '; ')
          
          aacigar_haplotype = paste(aacigar_haplotype[order(aacigar_haplotype$position),'aacigar_haplotype'], collapse = " ")
          
          
          # Update for each amplicon in the corresponding gene
          if(is.na(aacigar_table[sample, gene])){
            
            aacigar_table[sample, gene] = aacigar_haplotype
            phenotype_table[sample, gene] = sample_phenotype
            
          }else{
            
            aacigar_table[sample, gene] = paste(aacigar_table[sample, gene], aacigar_haplotype, sep = ' ')
            phenotype_table[sample, gene] = paste(phenotype_table[sample, gene], sample_phenotype, sep = '; ')
            
          }
          
        }
        
      }
      
    }
  }
  
  
  # Remove duplicates from phenotype_table
  print("Summarizing phenotype_table")
  for(sample in rownames(phenotype_table)){
    
    for(gene in colnames(phenotype_table)){
      
      phenotype = strsplit(phenotype_table[sample, gene], '; ')[[1]]
      
      if(sum(grepl('\\|', phenotype)) > 0){# if sample is heterozygous for the gene
        
        phenotype_clone1 = gsub('\\|.+','',phenotype)
        phenotype_clone2 = gsub('.+\\|','',phenotype)
        
        phenotypes = rbind(phenotype_clone1, phenotype_clone2)
        
        unique_phenotypes = NULL
        
        for(phenotype in 1:nrow(phenotypes)){
          
          phenotype = phenotypes[phenotype,]
          temp_unique_phenotypes = unique(phenotype)
          
          
          if(sum(grepl('resistance', temp_unique_phenotypes))>0){
            
            resistance_phenotypes = temp_unique_phenotypes[grepl('resistance', temp_unique_phenotypes)]
            for(temp_resistance_phenotype in resistance_phenotypes){
              number_of_DR_mutations = sum(grepl(temp_resistance_phenotype, phenotype))
              
              if(number_of_DR_mutations == 1){
                resistance_phenotype = paste0(number_of_DR_mutations,
                                              ' ',
                                              temp_unique_phenotypes[grepl(temp_resistance_phenotype, temp_unique_phenotypes)],
                                              ' mutation')
                temp_unique_phenotypes[grepl(temp_resistance_phenotype, temp_unique_phenotypes)] = resistance_phenotype
                
              }else if(number_of_DR_mutations > 1){
                resistance_phenotype = paste0(number_of_DR_mutations,
                                              ' ',
                                              temp_unique_phenotypes[grepl(temp_resistance_phenotype, temp_unique_phenotypes)],
                                              ' mutations')
                temp_unique_phenotypes[grepl(temp_resistance_phenotype, temp_unique_phenotypes)] = resistance_phenotype
              }
              
            }
          }
          
          
          if(sum(grepl('did not amplify',temp_unique_phenotypes)) >= 1 &
             sum(grepl('did not amplify',temp_unique_phenotypes)) < length(temp_unique_phenotypes)){
            
            not_amplified_amplicons = sort(gsub(' amplicon did not amplify', '', temp_unique_phenotypes[grepl('did not amplify',temp_unique_phenotypes)]))
            
            temp_unique_phenotypes[grepl('did not amplify',temp_unique_phenotypes)] = paste0('Amplicon(s) ', paste(not_amplified_amplicons, collapse = ','), ' for gene ', gene, " did not amplify" )
            
            temp_unique_phenotypes = unique(temp_unique_phenotypes)
            
          }else if(sum(grepl('did not amplify',temp_unique_phenotypes)) == length(temp_unique_phenotypes)){
            temp_unique_phenotypes = paste0('Gene ', gene, ' did not amplify')
          }
          
          if(length(temp_unique_phenotypes) > 1 & sum(grepl('Sensitive', temp_unique_phenotypes)) > 0){
            temp_unique_phenotypes = temp_unique_phenotypes[!grepl('Sensitive', temp_unique_phenotypes)]
          }
          
          
          temp_unique_phenotypes = paste(temp_unique_phenotypes, collapse = '; ')
          
          if(is.null(unique_phenotypes)){
            unique_phenotypes = temp_unique_phenotypes
          }else{
            unique_phenotypes = paste(unique_phenotypes, temp_unique_phenotypes, sep = "|")
          }
          
          
        }
        
        phenotype_table[sample, gene] = unique_phenotypes
        
      }else{# if sample is homozygous for the gene
        
        
        unique_phenotypes = unique(phenotype)
        
        
        if(sum(grepl('resistance', unique_phenotypes))>0){
          
          resistance_phenotypes = unique_phenotypes[grepl('resistance', unique_phenotypes)]
          for(temp_resistance_phenotype in resistance_phenotypes){
            number_of_DR_mutations = sum(grepl(temp_resistance_phenotype, phenotype))
            
            if(number_of_DR_mutations == 1){
              resistance_phenotype = paste0(number_of_DR_mutations,
                                            ' ',
                                            unique_phenotypes[grepl(temp_resistance_phenotype, unique_phenotypes)],
                                            ' mutation')
              unique_phenotypes[grepl(temp_resistance_phenotype, unique_phenotypes)] = resistance_phenotype
              
            }else if(number_of_DR_mutations > 1){
              resistance_phenotype = paste0(number_of_DR_mutations,
                                            ' ',
                                            unique_phenotypes[grepl(temp_resistance_phenotype, unique_phenotypes)],
                                            ' mutations')
              unique_phenotypes[grepl(temp_resistance_phenotype, unique_phenotypes)] = resistance_phenotype
            }
            
          }
        }
        
        
        if(sum(grepl('did not amplify',unique_phenotypes)) >= 1 &
           sum(grepl('did not amplify',unique_phenotypes)) < length(unique_phenotypes)){
          
          not_amplified_amplicons = sort(gsub(' amplicon did not amplify', '', unique_phenotypes[grepl('did not amplify',unique_phenotypes)]))
          
          unique_phenotypes[grepl('did not amplify',unique_phenotypes)] = paste0('Amplicon(s) ', paste(not_amplified_amplicons, collapse = ','), ' for gene ', gene, " did not amplify" )
          
          unique_phenotypes = unique(unique_phenotypes)
          
        }else if(sum(grepl('did not amplify',unique_phenotypes)) == length(unique_phenotypes)){
          unique_phenotypes = paste0('Gene ', gene, ' did not amplify')
        }
        
        if(length(unique_phenotypes) > 1 & sum(grepl('Sensitive', unique_phenotypes)) > 0){
          unique_phenotypes = unique_phenotypes[!grepl('Sensitive', unique_phenotypes)]
        }
        
        
        unique_phenotypes = paste(unique_phenotypes, collapse = '; ')
        
        phenotype_table[sample, gene] = unique_phenotypes
        
        
        
      }
      
      
      
    }
    
    
  }
  
  # Match Genotypes with Phenotypes
  
  print("Match Genotypes with Phenotypes")
  genotype_phenotype_table = data.frame(Sample = rownames(aacigar_table),
                                        Gene = rep(colnames(aacigar_table), each = nrow(aacigar_table)),
                                        Genotype = c(aacigar_table),
                                        Phenotype = c(phenotype_table))
  
  
  genotype_phenotype_match = genotype_phenotype_table[!duplicated(genotype_phenotype_table$Genotype),-1]
  
  
  genotype_phenotype_match_mono = genotype_phenotype_match[!grepl('\\|',genotype_phenotype_match$Genotype),]
  
  genotype_phenotype_match_poly1 = genotype_phenotype_match[grepl('\\|',genotype_phenotype_match$Genotype),]
  genotype_phenotype_match_poly1$Genotype = gsub('\\|[A-z]', '', genotype_phenotype_match_poly1$Genotype)
  genotype_phenotype_match_poly1$Phenotype = gsub('\\|.+', '', genotype_phenotype_match_poly1$Phenotype)
  
  
  genotype_phenotype_match_poly2 = genotype_phenotype_match[grepl('\\|',genotype_phenotype_match$Genotype),]
  genotype_phenotype_match_poly2$Genotype = gsub('[A-z]\\|', '', genotype_phenotype_match_poly2$Genotype)
  genotype_phenotype_match_poly2$Phenotype = gsub('.+\\|', '', genotype_phenotype_match_poly2$Phenotype)
  
  
  genotype_phenotype_match = rbind(genotype_phenotype_match_mono,
                                   genotype_phenotype_match_poly1,
                                   genotype_phenotype_match_poly2)
  
  
  genotype_phenotype_match = genotype_phenotype_match[!duplicated(genotype_phenotype_match$Genotype),]
  
  print("Handeling polyclonal samples for haplotype count")
  
  mon_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) == 0),]
  
  poly_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) != 0),]
  
  poly1 = gsub('\\|([a-z])', '', poly_aacigar_table, ignore.case = T)
  
  poly2 = gsub('([a-z]\\|)', '', poly_aacigar_table, ignore.case = T)
  
  extended_aacigar_table = rbind(mon_aacigar_table,
                                 poly1,
                                 poly2)
  
  extended_aacigar_table = data.frame(Sample_id = rownames(extended_aacigar_table), extended_aacigar_table)
  
  print("Adding metadata to haplotype counts")
  
  if(!is.null(Longitude) & !is.null(Latitude)){
    metadata = ampseq_object@metadata[,c(variables, Longitude, Latitude)]
    names(metadata) = c('Sample_id', 'var1', 'var2', 'Longitude', 'Latitude')
    
  }else{
    metadata = ampseq_object@metadata[,c(variables)]
    names(metadata) = c('Sample_id', 'var1', 'var2')
  }
  
  metadata[['var1']] = as.factor(metadata[['var1']])
  metadata[['var2']] = as.factor(metadata[['var2']])
  
  extended_aacigar_table = merge(extended_aacigar_table,
                                 metadata,
                                 by = 'Sample_id',
                                 all.x = TRUE)
  
  extended_aacigar_table %<>% pivot_longer(cols = all_of(gene_ids),
                                           names_to = 'gene_ids',
                                           values_to = 'haplotype')
  
  extended_aacigar_table$gene_names = NA
  
  for(gene in 1:length(gene_ids)){
    extended_aacigar_table[extended_aacigar_table$gene_ids == gene_ids[gene],][['gene_names']] = gene_names[gene]
    
  }
  
  
  if(na.var.rm){
    extended_aacigar_table = extended_aacigar_table[!is.na(extended_aacigar_table$var1) &
                                                      (!is.na(extended_aacigar_table$var2) &
                                                         !grepl('NA',extended_aacigar_table$var2)),]
  }else{
    extended_aacigar_table %<>% mutate(var1 = case_when(
      is.na(var1) ~ paste(variables[2], 'missing'),
      !is.na(var1) ~ var1),
      var2 = case_when(
        is.na(var2)|grepl('NA',var2) ~ paste(variables[3], 'missing'),
        (!is.na(var2))&(!grepl('NA',var2)) ~ var2))
  }
  
  
  # samples_pop_quarter = extended_aacigar_table %>%
  #   summarise(count = nlevels(as.factor(Sample_id)), .by = c(var1, var2)) 
  
  print("Summarizing haplotype counts")
  
  haplotype_counts = extended_aacigar_table %>%
    summarise(count = n(), .by = c(gene_names, var1, var2, haplotype))
  
  # 
  # print("Removing undesired categories based on var_filter")
  # if(!is.null(filters)){
  #   filters = strsplit(filters,';')
  #   for(temp_filter in 1:length(filters)){
  #     if(which(variables == filters[[temp_filter]][1]) == 2){
  #       haplotype_counts %<>% filter(var1 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
  #       samples_pop_quarter %<>% filter(var1 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
  #     }else if(which(variables == filters[[temp_filter]][1]) == 3){
  #       haplotype_counts %<>% filter(var2 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
  #       samples_pop_quarter %<>% filter(var2 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
  #     }
  #   }
  # }
  
  print('Calculating haplotype frequencies')
  haplotype_counts$freq = NA
  haplotype_counts$freq_lower = NA
  haplotype_counts$freq_upper = NA
  
  
  for(gene in levels(as.factor(haplotype_counts$gene_names))){
    for(Pop in levels(as.factor(haplotype_counts[haplotype_counts$gene_names == gene,][['var1']]))){
      for(date in levels(as.factor(haplotype_counts[haplotype_counts$gene_names == gene&
                                                    haplotype_counts$var1 == Pop,][['var2']]))){
        
        temp_freq = binconf(haplotype_counts[haplotype_counts$gene_names == gene&
                                               haplotype_counts$var1 == Pop&
                                               haplotype_counts$var2 == date,][['count']],
                            sum(haplotype_counts[haplotype_counts$gene_names == gene&
                                                   haplotype_counts$var1 == Pop&
                                                   haplotype_counts$var2 == date,][['count']]),
                            method = 'exact'
        )
        
        haplotype_counts[haplotype_counts$gene_names == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['freq']] = temp_freq[,1]
        
        haplotype_counts[haplotype_counts$gene_names == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['freq_lower']] = temp_freq[,2]
        
        haplotype_counts[haplotype_counts$gene_names == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['freq_upper']] = temp_freq[,3]
        
      }
    }
  }
  
  haplotype_counts %<>% mutate(gene_haplo = paste(gene_names, haplotype, sep = ": "))
  
  haplotype_counts$phenotype = NA
  
  for(haplotype in haplotype_counts$haplotype){
    haplotype_counts[haplotype_counts$haplotype == haplotype, ][['phenotype']] =
      genotype_phenotype_match[genotype_phenotype_match$Genotype == haplotype, ][['Phenotype']]
  }
  
  if(hap_color_palette == 'auto'){
    
    print('Assigning colors to haplotypes based on their phenotype')
    blue = brewer.pal(9, 'Blues')[7]
    reds = brewer.pal(9, 'Reds')[3:8]
    orange = brewer.pal(9, 'YlOrRd')[4]
    
    genotype_phenotype_match$color_pal = sapply(genotype_phenotype_match$Phenotype, function(Phenotype){
      
      if(grepl('Delayed clearance', Phenotype)){
        reds[1]
      }else if(grepl('Linked', Phenotype)){
        reds[1]
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| )+ resistance')), '\\d+')))) == 1){
        reds[2]
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| )+ resistance')), '\\d+')))) == 2){
        reds[3]
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| )+ resistance')), '\\d+')))) == 3){
        reds[5]
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| )+ resistance')), '\\d+')))) >= 4){
        reds[6]
      }else if(grepl('Sensitive', Phenotype)){
        blue
      }else if(grepl('polymorphism', Phenotype) &
               !grepl('resistance', Phenotype)){
        'gold3'
      }else if(grepl('variant unreported', Phenotype) & !grepl('resistance', Phenotype)){
        orange
      }else if(grepl("^Amplicon.+amplify; $", Phenotype)){
        'gray70'
      }else if(grepl('Gene .+ did not amplify', Phenotype)){
        'gray30'
      }
      
    }, simplify = T)
    
    
    genotype_phenotype_match$transparency = sapply(genotype_phenotype_match$Phenotype, function(Phenotype){
      transparency = 0.5 + 0.1*length(unlist(str_extract_all(Phenotype, '(variant unreported|polymorphism in gene)')))
      
      if(transparency > 1){
        
        transparency = 1
        
      }
      
      return(transparency)
      
    }, simplify = T)
    
    # Sort haplotypes based on gene and phenotype
    
    genotype_phenotype_match_sorted = NULL
    
    for(gene in unique(genotype_phenotype_match$Gene)){
      
      for(color_pal in c(reds[6:1], orange, 'gold3', blue, 'gray70','gray30')){
        
        genotype_phenotype_match_sorted = rbind(genotype_phenotype_match_sorted,
                                                genotype_phenotype_match[genotype_phenotype_match$Gene == gene &
                                                                           genotype_phenotype_match$color_pal == color_pal,])
        
      }
      
    }
    
  }else if(hap_color_palette == 'random'){
    
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    #pie(rep(1,n), col=sample(col_vector, n))
    
    
    genotype_phenotype_match$haplo_order = unlist(sapply(genotype_phenotype_match$Phenotype, function(Phenotype){
      
      if(grepl('Delayed clearance', Phenotype)){
        5
      }else if(grepl('Linked', Phenotype)){
        5
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| )+ resistance')), '\\d+')))) == 1){
        4
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| )+ resistance')), '\\d+')))) == 2){
        3
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| )+ resistance')), '\\d+')))) == 3){
        2
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| )+ resistance')), '\\d+')))) >= 4){
        1
      }else if(grepl('Sensitive', Phenotype)){
        8
      }else if(grepl('polymorphism', Phenotype) &
               !grepl('resistance', Phenotype)){
        7
      }else if(grepl('variant unreported', Phenotype) & !grepl('resistance', Phenotype)){
        6
      }else if(grepl("^Amplicon.+amplify(; )?$", Phenotype)){
        9
      }else if(grepl('Gene .+ did not amplify', Phenotype)){
        10
      }
      
    }, simplify = T))
    
    genotype_phenotype_match$color_pal = sample(col_vector, nrow(genotype_phenotype_match))
    
    genotype_phenotype_match$transparency = 1
    
    genotype_phenotype_match_sorted = NULL
    
    for(gene in unique(genotype_phenotype_match$Gene)){
      
      for(position in sort(unlist(unique(genotype_phenotype_match$haplo_order)))){
        
        genotype_phenotype_match_sorted = rbind(genotype_phenotype_match_sorted,
                                                genotype_phenotype_match[genotype_phenotype_match$Gene == gene &
                                                                           genotype_phenotype_match$haplo_order == position,])
        
      }
      
    }
    
    genotype_phenotype_match_sorted$haplo_order = NULL
    
  }
  
  
  genotype_phenotype_match_sorted$Gene_name = NA
  
  for(gene in gene_ids){
    
    genotype_phenotype_match_sorted[genotype_phenotype_match_sorted$Gene == gene, ][['Gene_name']] = gene_names[gene_ids == gene]
    
  }
  
  
  genotype_phenotype_match_sorted %<>% mutate(
    gene_haplo = paste(Gene_name, Genotype, sep = ': ')
  )
  
  
  genotype_phenotype_match_sorted = 
    genotype_phenotype_match_sorted[genotype_phenotype_match_sorted$gene_haplo %in%
                                      unique(haplotype_counts$gene_haplo ),]
  
  
  haplotype_counts = left_join(haplotype_counts,
                               genotype_phenotype_match_sorted[,c("gene_haplo",
                                                                  "transparency")],
                               by = 'gene_haplo')
  
  
  haplotype_counts$gene_haplo = factor(haplotype_counts$gene_haplo,
                                       levels = genotype_phenotype_match_sorted$gene_haplo)
  
  print('haplotype_freq_barplot')
  haplotype_freq_barplot = haplotype_counts %>%
    ggplot(aes(y = freq, x = var2, fill  = gene_haplo, alpha = gene_haplo)) +
    geom_col(color = 'gray75')+
    facet_grid(var1 ~ gene_names)+
    scale_fill_manual(values = genotype_phenotype_match_sorted$color_pal)+
    scale_alpha_manual(breaks = genotype_phenotype_match_sorted$gene_haplo, values = genotype_phenotype_match_sorted$transparency)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = 'bottom') +
    labs(y = 'Frequency in population',
         x = 'Date of Collection',
         fill = 'Gene: Haplotype')+
    guides(alpha = "none",
           fill=guide_legend(ncol=3))
  
  
  print('haplotypes_freq_lineplot')
  haplotypes_freq_lineplot = haplotype_counts %>%
    ggplot(aes(y = freq, x = var2, group  = gene_haplo, color = gene_haplo)) +
    geom_point()+
    geom_errorbar(aes(ymin = freq_lower, ymax = freq_upper), alpha = .5, width = .2)+
    geom_line()+
    facet_grid(var1 ~ gene_names)+
    scale_color_manual(values = genotype_phenotype_match_sorted$color_pal)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = 'bottom') +
    labs(y = 'Frequency in population',
         x = 'Date of Collection',
         color = 'Gene: Haplotype')+
    guides(alpha = "none",
           color=guide_legend(ncol=3))
  
  #names(haplotype_counts) = c(names(haplotype_counts)[1], variables[2:3], names(haplotype_counts)[-1:-3])
  
  print('Defining phenotypes based to drug of interest')
  
  drug_phenotype_table = data.frame(Sample_id = rownames(aacigar_table))
  
  for(drug in drugs){
    
    genes = unique(drugR_reference_alleles[grepl(drug, drugR_reference_alleles$Annotation),][['Gene_Id']])
    
    genes = genes[genes %in% gene_ids]
    
    temp_drug_phenotype_table = as.data.frame(drug_phenotype_table[,'Sample_id'])
    
    for(gene in genes){
      
      temp_drug_phenotype_table = cbind(temp_drug_phenotype_table,
                                        genotype_phenotype_table[genotype_phenotype_table$Gene == gene, ][['Phenotype']])
      
    }
    
    colnames(temp_drug_phenotype_table) = c('Sample_id', genes)
    
    
    rownames(temp_drug_phenotype_table) = drug_phenotype_table$Sample
    
    temp_drug_phenotype_table[[drug]] = NA
    
    for(sample in rownames(temp_drug_phenotype_table)){
      
      phenotype = temp_drug_phenotype_table[sample,genes]
      
      # Check for resistance associated mutants
      
      if(sum(grepl(paste0('\\d+ ([A-z]+ and )?', drug, '( and [A-z]+)? resistance mutations?'), phenotype)) > 0){
        
        resistance_phenotype_genes = which(grepl(paste0('\\d+ ([A-z]+ and )?', drug, '( and [A-z]+)? resistance mutations?'), phenotype))
        
        for(resistante_gene in resistance_phenotype_genes){
          
          gene_phenotype = str_extract(temp_drug_phenotype_table[sample,][[genes[resistante_gene]]], paste0('\\d+ ([A-z]+ and )?', drug, '( and [A-z]+)? resistance mutations?'))
          
          gene_phenotype = paste0(gene_phenotype, ' in gene ', genes[resistante_gene])
          
          if(is.na(temp_drug_phenotype_table[sample,][[drug]])){
            
            temp_drug_phenotype_table[sample,][[drug]] = gene_phenotype
            
          }else{
            
            temp_drug_phenotype_table[sample,][[drug]] = paste(temp_drug_phenotype_table[sample,][[drug]],
                                                               gene_phenotype, sep = '; ')
            
          }
          
        }
        
      }
      
      # Check for not reported mutants
      
      if(sum(grepl('polymorphism', phenotype)) > 0){
        
        unreported_mutations = which(grepl('polymorphism', phenotype))
        
        for(unreported_mutation in unreported_mutations){
          
          gene_phenotype = str_extract(temp_drug_phenotype_table[sample,][[genes[unreported_mutation]]], paste0('[a-z]\\d+[a-z] polymorphism in gene ', genes[unreported_mutation]))
          
          
          if(is.na(temp_drug_phenotype_table[sample,][[drug]])){
            
            temp_drug_phenotype_table[sample,][[drug]] = gene_phenotype
            
          }else{
            
            temp_drug_phenotype_table[sample,][[drug]] = paste(temp_drug_phenotype_table[sample,][[drug]],
                                                               gene_phenotype, sep = '; ')
            
          }
          
        }
        
      }
      
      # Check for sensitive phenotypes
      
      if(sum(grepl('Sensitive', phenotype)) == length(phenotype)){
        
        temp_drug_phenotype_table[sample,][[drug]] = paste0('Sensitive phenotype based on gene(s) ',
                                                            ifelse(length(genes) > 1,
                                                                   paste(paste(genes[-length(genes)], collapse = ', '), genes[length(genes)], sep = ' and '),
                                                                   genes)
        )
        
      }
      
      # Check for complete missing data
      
      if(sum(grepl('Gene .+ did not amplified', phenotype)) == length(phenotype)){
        
        temp_drug_phenotype_table[sample,][[drug]] = paste0('Gene(s) ',
                                                            ifelse(length(genes) > 1,
                                                                   paste(paste(genes[-length(genes)], collapse = ', '), genes[length(genes)], sep = ' and '),
                                                                   genes),
                                                            ' did not amplified'
        )
        
      }
      
      
      # Check for partial haplotypes
      
      if(is.na(temp_drug_phenotype_table[sample,][[drug]])){
        
        
        partial_haplotypes = phenotype[which(!grepl('Sensitive', phenotype))]
        
        
        temp_drug_phenotype_table[sample,][[drug]] = paste0('Partial haplotype: ',
                                                            ifelse(length(partial_haplotypes) > 1,
                                                                   paste(paste(partial_haplotypes[-length(partial_haplotypes)], collapse = ', '), partial_haplotypes[length(partial_haplotypes)], sep = ' and '),
                                                                   partial_haplotypes))
        
        
      }
      
      
      
    }
    
    drug_phenotype_table[[drug]] = temp_drug_phenotype_table[[drug]]
    
  }
  
  
  drug_phenotype_table %<>% pivot_longer(cols = all_of(drugs), names_to = 'Drug', values_to = 'Phenotype')
  
  
  drug_phenotype_summary = drug_phenotype_table %>% mutate(
    Phenotype = case_when(
      grepl(paste0('\\d+ ([A-z]+ and )?(', paste(drugs, collapse = '|'), ')( and [A-z]+)? resistance mutations?'), Phenotype) ~ "Mutation(s) associated with a resistant phenotype",
      grepl(paste0('^[a-z]\\d+[a-z] polymorphism in gene (', paste(gene_ids, collapse = '|'), ')'), Phenotype) ~ "Polymorphism(s) respect to reference strain",
      grepl('Sensitive phenotype', Phenotype) ~ "Sensitive Phenotype",
      grepl('Partial haplotype:', Phenotype) ~ "Partial Haplotype",
      grepl('Gene\\(s\\).+did not amplified', Phenotype) ~ "Missing data"
    ))
  
  
  drug_phenotype_summary = merge(drug_phenotype_summary,
                                 metadata,
                                 by = 'Sample_id',
                                 all.x = TRUE)
  
  
  if(na.var.rm){
    drug_phenotype_summary = drug_phenotype_summary[!is.na(drug_phenotype_summary$var1) &
                                                      (!is.na(drug_phenotype_summary$var2) &
                                                         !grepl('NA',drug_phenotype_summary$var2)),]
  }else{
    drug_phenotype_summary %<>% mutate(var1 = case_when(
      is.na(var1) ~ paste(variables[2], 'missing'),
      !is.na(var1) ~ var1),
      var2 = case_when(
        is.na(var2)|grepl('NA',var2) ~ paste(variables[3], 'missing'),
        (!is.na(var2))&(!grepl('NA',var2)) ~ var2))
  }
  
  
  if(na.hap.rm){
    drug_phenotype_summary = drug_phenotype_summary[!grepl('(Partial|Missing)',drug_phenotype_summary$Phenotype),]
  }
  
  
  # samples_pop_quarter = drug_phenotype_summary %>%
  #   summarise(count = nlevels(as.factor(Sample_id)), .by = c(var1, var2)) 
  
  if(!is.null(Longitude) & !is.null(Latitude)){
    
    drug_phenotype_summary = drug_phenotype_summary %>%
      summarise(count = n(),
                Longitude = mean(Longitude),
                Latitude = mean(Latitude),
                .by = c(Drug, var1, var2, Phenotype))
    
  }else{
    
    drug_phenotype_summary = drug_phenotype_summary %>%
      summarise(count = n(),
                .by = c(Drug, var1, var2, Phenotype))
    
  }
  
  
  
  # if(!is.null(filters)){
  #   #filters = strsplit(filters,';')
  #   for(temp_filter in 1:length(filters)){
  #     if(which(variables == filters[[temp_filter]][1]) == 2){
  #       drug_phenotype_summary %<>% filter(var1 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
  #       samples_pop_quarter %<>% filter(var1 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
  #     }else if(which(variables == filters[[temp_filter]][1]) == 3){
  #       drug_phenotype_summary %<>% filter(var2 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
  #       samples_pop_quarter %<>% filter(var2 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
  #     }
  #   }
  # }
  
  drug_phenotype_summary$freq = NA
  drug_phenotype_summary$freq_lower = NA
  drug_phenotype_summary$freq_upper = NA
  
  
  for(drug in levels(as.factor(drug_phenotype_summary$Drug))){
    for(Pop in levels(as.factor(drug_phenotype_summary[drug_phenotype_summary$Drug == drug,][['var1']]))){
      for(date in levels(as.factor(drug_phenotype_summary[drug_phenotype_summary$Drug == drug&
                                                          drug_phenotype_summary$var1 == Pop,][['var2']]))){
        
        temp_freq = binconf(drug_phenotype_summary[drug_phenotype_summary$Drug == drug&
                                                     drug_phenotype_summary$var1 == Pop&
                                                     drug_phenotype_summary$var2 == date,][['count']],
                            sum(drug_phenotype_summary[drug_phenotype_summary$Drug == drug&
                                                         drug_phenotype_summary$var1 == Pop&
                                                         drug_phenotype_summary$var2 == date,][['count']]),
                            method = 'exact'
        )
        
        drug_phenotype_summary[drug_phenotype_summary$Drug == drug&
                                 drug_phenotype_summary$var1 == Pop&
                                 drug_phenotype_summary$var2 == date,][['freq']] = temp_freq[,1]
        
        drug_phenotype_summary[drug_phenotype_summary$Drug == drug&
                                 drug_phenotype_summary$var1 == Pop&
                                 drug_phenotype_summary$var2 == date,][['freq_lower']] = temp_freq[,2]
        
        drug_phenotype_summary[drug_phenotype_summary$Drug == drug&
                                 drug_phenotype_summary$var1 == Pop&
                                 drug_phenotype_summary$var2 == date,][['freq_upper']] = temp_freq[,3]
        
        
      }
    }
  }
  
  print('drug_phenotype_barplot')
  drug_phenotype_barplot = drug_phenotype_summary %>%
    ggplot(aes(y = freq, x = var2, fill  = factor(Phenotype,
                                                  levels = c("Mutation(s) associated with a resistant phenotype",
                                                             "Polymorphism(s) respect to reference strain",
                                                             "Sensitive Phenotype",
                                                             "Partial Haplotype",
                                                             "Missing data")))) +
    geom_col()+
    facet_grid(var1 ~ Drug)+
    scale_fill_manual(values = c('firebrick3', 'gold3', 'dodgerblue3', 'gray70', 'gray30'))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = 'bottom') +
    labs(y = 'Pop. Frequency',
         x = 'Date of Collection',
         fill = 'Phenotype')
  
  
  print('drug_phenotyope_lineplot')
  drug_phenotyope_lineplot = drug_phenotype_summary %>%
    ggplot(aes(y = freq, x = var2, group  = factor(Phenotype,
                                                   levels = c("Mutation(s) associated with a resistant phenotype",
                                                              "Polymorphism(s) respect to reference strain",
                                                              "Sensitive Phenotype",
                                                              "Partial Haplotype",
                                                              "Missing data")),
               color = factor(Phenotype,
                              levels = c("Mutation(s) associated with a resistant phenotype",
                                         "Polymorphism(s) respect to reference strain",
                                         "Sensitive Phenotype",
                                         "Partial Haplotype",
                                         "Missing data")))) +
    geom_point()+
    geom_errorbar(aes(ymin = freq_lower, ymax = freq_upper), alpha = .5, width = .2)+
    geom_line()+
    facet_grid(var1 ~ Drug)+
    scale_color_manual(values = c('firebrick3', 'gold3', 'dodgerblue3', 'gray70', 'gray30'))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = 'bottom') +
    labs(y = 'Pop. Frequency',
         x = 'Date of Collection',
         color = 'Phenotype')
  
  
  
  print("Estimating frequency for drug resistant phenotypes")
  
  if(!is.null(Longitude) & !is.null(Latitude)){
    
    drug_phenotype_summary_sdf = drug_phenotype_summary %>% 
      summarise(Longitude = mean(Longitude),
                Latitude = mean(Latitude),
                count = sum(count), .by = c(Drug, Phenotype, var1))
    
    
    drug_phenotype_summary_sdf$ssize = NA
    drug_phenotype_summary_sdf$freq = NA
    drug_phenotype_summary_sdf$freq_lower = NA
    drug_phenotype_summary_sdf$freq_upper = NA
    
    
    for(drug in levels(as.factor(drug_phenotype_summary_sdf$Drug))){
      for(Pop in levels(as.factor(drug_phenotype_summary_sdf[drug_phenotype_summary_sdf$Drug == drug,][['var1']]))){
        
        ssize = sum(drug_phenotype_summary_sdf[drug_phenotype_summary_sdf$Drug == drug &
                                                 drug_phenotype_summary_sdf$var1 == Pop,][['count']])
        
        temp_freq = binconf(drug_phenotype_summary_sdf[drug_phenotype_summary_sdf$Drug == drug &
                                                         drug_phenotype_summary_sdf$var1 == Pop,][['count']],
                            ssize,
                            method = 'exact'
        )
        
        drug_phenotype_summary_sdf[drug_phenotype_summary_sdf$Drug == drug &
                                     drug_phenotype_summary_sdf$var1 == Pop,][['ssize']] = ssize
        
        drug_phenotype_summary_sdf[drug_phenotype_summary_sdf$Drug == drug &
                                     drug_phenotype_summary_sdf$var1 == Pop,][['freq']] = round(temp_freq[,1], 2)
        
        drug_phenotype_summary_sdf[drug_phenotype_summary_sdf$Drug == drug &
                                     drug_phenotype_summary_sdf$var1 == Pop,][['freq_lower']] = round(temp_freq[,2], 2)
        
        drug_phenotype_summary_sdf[drug_phenotype_summary_sdf$Drug == drug &
                                     drug_phenotype_summary_sdf$var1 == Pop,][['freq_upper']] = round(temp_freq[,3], 2)
        
        
      }
    }
    
    
    drug_phenotype_summary_sdf %<>% filter(!is.na(Longitude), !is.na(Latitude))
    
    drug_phenotype_summary_sdf %<>% mutate(logssize = log(ssize, 1.3))
    
    drug_phenotype_summary_sdf %<>% filter(Phenotype == "Mutation(s) associated with a resistant phenotype")
    
    
    print("Transforming data to spatial points")
    drug_phenotype_summary_sdf = SpatialPointsDataFrame(coords = drug_phenotype_summary_sdf[,c("Longitude", "Latitude")],
                                                        data = drug_phenotype_summary_sdf,
                                                        proj4string = CRS("+init=epsg:4326"))
    
    
    tmap_mode('view')
    print('i_drug_map')
    i_drug_map = tm_shape(drug_phenotype_summary_sdf)+
      tm_dots(size = "logssize", style="pretty", col = "freq")+
      tm_text("freq", size=1)+
      tm_facets(by = "Drug")+
      tm_scale_bar()
    
    
    drug_resistant_hap_list = list(aa_mutations = haps_respect_to_ref$loci_aa_table,
                                   dna_mutations = haps_respect_to_ref$loci_dna_table,
                                   genotype_phenotype_table = genotype_phenotype_table,
                                   drug_phenotype_table = drug_phenotype_table,
                                   drug_phenotyope_lineplot = drug_phenotyope_lineplot,
                                   drug_phenotype_barplot = drug_phenotype_barplot,
                                   drug_phenotype_summary_sdf = drug_phenotype_summary_sdf,
                                   i_drug_map = i_drug_map,
                                   haplotypes_freq_lineplot = haplotypes_freq_lineplot,
                                   haplotype_freq_barplot = haplotype_freq_barplot)
    
  }else{
    
    drug_resistant_hap_list = list(aa_mutations = haps_respect_to_ref$loci_aa_table,
                                   dna_mutations = haps_respect_to_ref$loci_dna_table,
                                   genotype_phenotype_table = genotype_phenotype_table,
                                   drug_phenotype_table = drug_phenotype_table,
                                   drug_phenotyope_lineplot = drug_phenotyope_lineplot,
                                   drug_phenotype_barplot = drug_phenotype_barplot,
                                   haplotypes_freq_lineplot = haplotypes_freq_lineplot,
                                   haplotype_freq_barplot = haplotype_freq_barplot)
  }
  
  
  
  return(drug_resistant_hap_list)
  
}



get_Fws = function(ampseq_object = NULL){
  
  gt = ampseq_object@gt
  
  ExpHet = get_loci_diversity(loci_abd_table = ampseq_object@gt, variance = F)
  
  ExpHet = ExpHet[['Hexp']]
  
  Hw = sapply(1:nrow(gt), function(sample){
    
    samp_alleles= gsub(':\\d+', '', gt[sample,])
    samp_allcounts = gsub('([A-Z]|\\d|\\.)+:', '', gt[sample,])
    
    samp_alleles1 = gsub('_([A-Z]|\\d|\\.)+$', '', samp_alleles)
    samp_alleles2 = gsub('^([A-Z]|\\d|\\.)+_', '', samp_alleles)
    
    samp_check = samp_alleles1 != samp_alleles2
    
    samp_alleles2[!samp_check] = NA
    
    samp_allcounts1 = gsub('_\\d+$', '', samp_allcounts)
    samp_allcounts2 = gsub('^\\d+_', '', samp_allcounts)
    
    samp_allcounts2[!samp_check] = NA
    
    samp_allcountsT = rowSums(cbind(as.integer(samp_allcounts1), as.integer(samp_allcounts2)), na.rm = T)
    
    samp_allfreq = cbind(as.integer(samp_allcounts1), as.integer(samp_allcounts2))/samp_allcountsT
    
    Hw = 1 - rowSums(samp_allfreq^2, na.rm = T)
    Hw[Hw==1] = NA
    Hw
  })
  
  Fws = 1 - (Hw/ExpHet)
  
  Fws = colMeans(Fws, na.rm = T)
  
  Fws[Fws < 0] = 0 
  
  return(Fws)
  
}



## Functions for Complexity of infection----

### get_polygenomic----

get_polygenomic = function(ampseq_object, strata = NULL, update_popsummary = T, na.rm = FALSE, filters = NULL, poly_quantile = .75, poly_formula = "NHetLoci >= 1 & Fws < 1"){
  
  
  if(grepl("(NHetLoci|Frac_HetLoci|max_nAlleles|Fws)(<|>|!|=)+", poly_formula)){
    stop("All mathematical and logical operators must be separated by blank spaces in mask_formula")
  }
  
  # modify poly_formula
  
  
  if(grepl("NHetLoci ", poly_formula)){
    
    mask_filter = str_extract(poly_formula, "NHetLoci (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|FracHetLoci_quantile|Fws_quantile)")
    
    if(!is.na(mask_filter)){
      print(paste0('Filter ', str_extract(poly_formula, "NHetLoci (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|FracHetLoci_quantile|Fws_quantile)"), ' will be applied'))
      poly_formula = gsub("NHetLoci ", "polygenomic[['NHetLoci']] ", poly_formula)
    }else{
      stop("Filter NHetLoci is been called but there are spelling issues in this part of the formula")
    }
  }
  
  
  if(grepl("Frac_HetLoci ", poly_formula)){
    
    mask_filter = str_extract(poly_formula, "Frac_HetLoci (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|FracHetLoci_quantile|Fws_quantile)")
    
    if(!is.na(mask_filter)){
      print(paste0('Filter ', str_extract(poly_formula, "Frac_HetLoci (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|FracHetLoci_quantile|Fws_quantile)"), ' will be applied'))
      poly_formula = gsub("Frac_HetLoci ", "polygenomic[['Frac_HetLoci']] ", poly_formula)
    }else{
      stop("Filter Frac_HetLoci is been called but there are spelling issues in this part of the formula")
    }
  }
  
  
  if(grepl("max_nAlleles ", poly_formula)){
    
    mask_filter = str_extract(poly_formula, "max_nAlleles (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|FracHetLoci_quantile|Fws_quantile)")
    
    if(!is.na(mask_filter)){
      print(paste0('Filter ', str_extract(poly_formula, "max_nAlleles (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|FracHetLoci_quantile|Fws_quantile)"), ' will be applied'))
      poly_formula = gsub("max_nAlleles ", "polygenomic[['max_nAlleles']] ", poly_formula)
    }else{
      stop("Filter max_nAlleles is been called but there are spelling issues in this part of the formula")
    }
  }
  
  
  if(grepl("Fws ", poly_formula)){
    
    mask_filter = str_extract(poly_formula, "Fws (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|FracHetLoci_quantile|Fws_quantile)")
    
    if(!is.na(mask_filter)){
      print(paste0('Filter ', str_extract(poly_formula, "Fws (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|FracHetLoci_quantile|Fws_quantile)"), ' will be applied'))
      poly_formula = gsub("Fws ", "polygenomic[['Fws']] ", poly_formula)
    }else{
      stop("Filter Fws is been called but there are spelling issues in this part of the formula")
    }
  }
  
  
  mask_formula_check = str_split(poly_formula, "&|\\|")[[1]]
  mask_formula_check  = mask_formula_check[!grepl("polygenomic", mask_formula_check)]
  
  
  if(length(mask_formula_check) > 0){
    for(wrong_filter in mask_formula_check){
      print(paste0("Spelling error with filter ", wrong_filter))
    }
    stop("Execution halted, revise mask_filter argument.\nPossible filters are:\nNHetLoci, Frac_HetLoci, max_nAlleles, Fws")
  }
  
  library(Hmisc)

  
  if(!is.null(strata)){
    if(na.rm){
      ampseq_object = filter_samples(ampseq_object, v = !(is.na(ampseq_object@metadata[[strata]]) | grepl('NA',ampseq_object@metadata[[strata]])))
    }else if(length(metadata[is.na(metadata[[strata]]) | grepl('NA',metadata[[strata]]),][[strata]])>0){
      ampseq_object@metadata[is.na(ampseq_object@metadata[[strata]]) | grepl('NA',ampseq_object@metadata[[strata]]),][[strata]] = 'missing data'
    }
    
    if(!is.null(filters)){
      ampseq_object = filter_samples(ampseq_object, v = grepl(filters,ampseq_object@metadata[[strata]]))
    }
    
  }
  
  gt = ampseq_object@gt
  metadata = ampseq_object@metadata
  loci_performance = ampseq_object@loci_performance
  
  polygenomic = NULL
  
  for(sample in rownames(gt)){
    polygenomic = rbind(polygenomic, data.frame(NHetLoci = sum(grepl("_",gt[sample, ])),
                                                Frac_HetLoci = sum(grepl("_",gt[sample, ]))/sum(!is.na(gt[sample, ])),
                                                Polyloci = paste(names(gt[sample, ])[which(grepl("_",gt[sample, ]))], collapse = "/"),
                                                alleles_at_loci = paste(gt[sample, ][which(grepl("_",gt[sample, ]))], collapse = "/"),
                                                nalleles_per_loci = paste(
                                                  sapply(
                                                    sapply(
                                                      gt[sample, ][which(grepl("_", gt[sample, ]))], function(x){strsplit(x, "_")}),
                                                    function(x) length(x)),
                                                  collapse = "/"),
                                                max_nAlleles = ifelse(sum(grepl("_",gt[sample, ])) == 0,
                                                             1,
                                                             max(unlist(sapply(
                                                               sapply(
                                                                 gt[sample, ][which(grepl("_",gt[sample, ]))], function(x){strsplit(x, "_")}),
                                                               function(x) length(x)))))))
  }
  
  
  
  
  polygenomic[['Fws']] = get_Fws(ampseq_object)
  
  
  FracHetLoci_quantile = quantile(polygenomic$Frac_HetLoci, poly_quantile, na.rm = T)
  Fws_quantile = quantile(polygenomic$Fws, probs = 1 - poly_quantile, na.rm = T)

  plot_fracHet_vs_Fws = ggdraw()+
    draw_plot(polygenomic %>%
                ggplot(aes(x = Frac_HetLoci, y = Fws,
                           color = as.character(max_nAlleles)
                ))+
                geom_point(alpha = .5)+
                geom_vline(xintercept = FracHetLoci_quantile, linetype = 2)+
                geom_hline(yintercept = Fws_quantile, linetype = 2)+
                scale_color_viridis_d()+
                #scale_x_continuous(limits = c(-0.05,1.05))+
                #scale_y_continuous(limits = c(-0.05,1.05))+
                labs(color = 'max_nAlleles')+
                theme_bw()+
                theme(axis.title = element_blank(),
                      axis.text = element_blank(),
                      legend.position = c(.8,.8)),
              x = .3,
              width = .7,
              y = .3,
              height = .7
              )+
    draw_plot(polygenomic %>%
                ggplot(aes(x = Frac_HetLoci))+
                geom_histogram(binwidth = .01)+
                geom_vline(xintercept = c(FracHetLoci_quantile), linetype = 2)+
                #scale_x_continuous(limits = c(-0.05,1.05))+
                theme_bw(),
              x = .225,
              width = .775,
              y = 0,
              height = .3)+
    draw_plot(polygenomic %>%
                ggplot(aes(x = Fws))+
                geom_histogram(binwidth = .01)+
                geom_vline(xintercept = c(Fws_quantile), linetype = 2)+
                #scale_x_continuous(limits = c(-0.05,1.05))+
                theme_bw()+
                coord_flip(),
              x = 0,
              width = .3,
              y = .25,
              height = .75)
  
  
  
  
  polyclonals = which(eval(parse(text = poly_formula)))
  
  polygenomic$Clonality = NA
  polygenomic[polyclonals,][['Clonality']] = 'Polyclonal'
  polygenomic[-1*polyclonals,][['Clonality']] = 'Monoclonal'
  
  pop_summary = data.frame(
    pop = "Total",
    n = nrow(polygenomic),
    #mean_coi = mean(polygenomic[["coi"]]),
    n_poly = nrow(polygenomic[polygenomic[["Clonality"]] == 'Polyclonal',]),
    prop_poly = binconf(x = nrow(polygenomic[polygenomic[["Clonality"]] == 'Polyclonal',]),
                        n = nrow(polygenomic),
                        method = "exact")[1],
    prop_poly_lower = binconf(x = nrow(polygenomic[polygenomic[["Clonality"]] == 'Polyclonal',]),
                              n = nrow(polygenomic),
                              method = "exact")[2],
    prop_poly_upper = binconf(x = nrow(polygenomic[polygenomic[["Clonality"]] == 'Polyclonal',]),
                              n = nrow(polygenomic),
                              method = "exact")[3])
  
  
  if(!is.null(strata)){
    
    for(pop in levels(as.factor(metadata[[strata]]))){
      pop_summary = rbind(pop_summary,
                          data.frame(
                            pop = pop,
                            n = nrow(polygenomic[metadata[[strata]] == pop,]),
                            #mean_coi = mean(polygenomic[metadata[[strata]] == pop,][["coi"]]),
                            n_poly = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["Clonality"]] == 'Polyclonal',]),
                            prop_poly = binconf(x = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["Clonality"]] == 'Polyclonal',]),
                                                n = nrow(polygenomic[metadata[[strata]] == pop,]),
                                                method = "exact")[1],
                            prop_poly_lower = binconf(x = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["Clonality"]] == 'Polyclonal',]),
                                                      n = nrow(polygenomic[metadata[[strata]] == pop,]),
                                                      method = "exact")[2],
                            prop_poly_upper = binconf(x = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["Clonality"]] == 'Polyclonal',]),
                                                      n = nrow(polygenomic[metadata[[strata]] == pop,]),
                                                      method = "exact")[3]))
    }
    
  }
  
  
  
  # Defining the proportion of polygenomic infections detected by loci --
  
  prop_poly_detected = apply(gt, 2, function(x) sum(grepl("_",x))/length(x))
  
  loci_performance[['prop_poly_detected']] = prop_poly_detected
  
  if(update_popsummary){
    ampseq_object@metadata = cbind(ampseq_object@metadata, polygenomic)
    ampseq_object@loci_performance = loci_performance
    ampseq_object@pop_summary = if(is.null(ampseq_object@pop_summary)){pop_summary}else{cbind(ampseq_object@pop_summary, pop_summary)}
    ampseq_object@plots[['plot_fracHet_vs_Fws']] = plot_fracHet_vs_Fws
    return(ampseq_object)
  }else{
    return(list(coi_bySample = data.frame(Sample_id = rownames(gt),polygenomic),
                pop_summary = pop_summary,
                coi_byLoci = data.frame(locus = names(prop_poly_detected), prop_poly_detected = prop_poly_detected),
                plot_fracHet_vs_Fws = plot_fracHet_vs_Fws))
  }
  
}


### draw_haplotypes----

draw_haplotypes = function(ampseq_object = NULL, comparison_list = NULL){
  
  coded_haplotypes = NULL
  
  for(comparison in 1:length(comparison_list)){
    
    temp_ampseq = filter_samples(ampseq_object, comparison_list[[comparison]])
    
    gt = temp_ampseq@gt
    gt = gsub(':\\d+', '', gt)
    
    for(mhap in 1:ncol(gt)){
      
      alleles = gt[,mhap]
      alleles = alleles[!is.na(alleles)]
      alleles = unlist(strsplit(alleles, '_')) 
      alleles = unique(alleles)
      
      if(length(alleles) > 0){
        alleles_labels = 1:length(alleles)
        
        for(allele in alleles_labels){
          
          if(alleles[allele] != '.'){
            
            for(sample in 1:nrow(gt)){
              
              pattern_to_replace = paste0('(^|_)',alleles[allele], '(_|$)')
              
              if(grepl(pattern_to_replace, gt[sample, mhap])){
                
                pattern_to_replace = str_extract(gt[sample, mhap], pattern_to_replace)
                replacement_pattern = gsub(alleles[allele], allele, pattern_to_replace)
                
                gt[sample, mhap] = gsub(pattern_to_replace, replacement_pattern, gt[sample, mhap])
                
              }
              
            }
            
          }else{
            gt[, mhap] = gsub('\\.', allele, gt[, mhap])
          }
          
        }
      }
    }
    
    clonality = sapply(1:nrow(gt),function(sample){
      
      mhaps = strsplit(gt[sample,], '_')
      max(sapply(mhaps, length))
      
    })
    
    monoclonals = clonality == 1
    polyclonals = clonality > 1
    
    gt_final = NULL
    
    for(sample in comparison_list[[comparison]]){
      
      if(monoclonals[which(comparison_list[[comparison]] == sample)]){
        temp_gt_polyclonal = matrix(gt[sample,],
                                    ncol = ncol(gt),
                                    nrow = 1,
                                    dimnames = list(
                                      sample,
                                      colnames(gt)
                                    ))
        
      }else if(polyclonals[which(comparison_list[[comparison]] == sample)]){
        poly_sample = which(comparison_list[[comparison]] == sample)
        
        temp_gt_polyclonal = matrix(NA, 
                                    nrow = clonality[poly_sample],
                                    ncol = ncol(gt),
                                    dimnames = list(
                                      paste0(comparison_list[[comparison]][poly_sample], '_C',1:clonality[poly_sample]),
                                      colnames(gt)
                                    ))
        
        alleles = strsplit(gt[poly_sample,], '_')
        
        for(clone in 1:clonality[poly_sample]){
          
          temp_gt_polyclonal[clone,] = unlist(sapply(alleles, function(mhap){
            
            if(length(mhap) < clone){
              mhap[1]
            }else{
              mhap[clone]
            }
            
          }))
          
        }
        
      }
      gt_final = rbind(gt_final, temp_gt_polyclonal)
      
    }
    
    coded_haplotypes_temp = data.frame(Marker = colnames(gt_final),
                                       t(gt_final))
    
    colnames(coded_haplotypes_temp) = c('Marker', rownames(gt_final))
    
    coded_haplotypes_temp = left_join(coded_haplotypes_temp, 
                                      ampseq_object@markers[,c('amplicon', 'chromosome', 'pos')],
                                      by = join_by('Marker' == 'amplicon'))
    
    coded_haplotypes_temp %<>% pivot_longer(cols = all_of(rownames(gt_final)), names_to = 'Haplotype', values_to = 'Allele')
    
    
    
    coded_haplotypes_temp$Haplotype = 
      factor(coded_haplotypes_temp$Haplotype,
             levels = rownames(gt_final))
    
    coded_haplotypes_temp$Marker = 
      factor(coded_haplotypes_temp$Marker,
             levels = unique(coded_haplotypes_temp$Marker))
    
    coded_haplotypes_temp$chromosome = gsub('(^[A-z]+\\d+_|_v\\d+$)', '', coded_haplotypes_temp$chromosome)
    
    
    coded_haplotypes_temp$Allele = as.character(coded_haplotypes_temp$Allele)
    
    coded_haplotypes_temp$Comparison = paste0('Comparison ', comparison)
    
    coded_haplotypes = rbind(coded_haplotypes, coded_haplotypes_temp)
    
  }
  
  haplotypes_plot = coded_haplotypes %>%
    ggplot(aes(x = Marker, y = Haplotype, fill = Allele)) +
    geom_tile()+
    scale_fill_brewer(palette = 'Set1')+
    facet_grid(Comparison~chromosome, scales = 'free', space = 'free')+
    theme(axis.text.x = element_blank())
  
  return(haplotypes_plot)
  
}



## Functions for IBD and other population structure metrics----
### fs_checks----

# Function to check if frequencies are ordered, sum to one, in [0,1], and informative.

# In hmmloglikelihood.cpp any frequency greater than 1e-20 is considered non-zero
fs_checks <- function(fs, warn,
                      non_zero_fs_lb = 0, #1e-20, # Magic number
                      max_dev_lim = 0, #1e-5, # Magic number
                      do_return = FALSE){
  
  if (!warn) {
    # turn off new warn option while returning old
    op <- options(warn = -1)
    on.exit(options(op))
  }
  
  # Convert fs into numeric logic for zero/non-zero
  fs01 <- 1 * (fs > non_zero_fs_lb)
  
  # Check frequencies are in [0,1]
  # Those that exceed one do not break the hmmloglikelihood.cpp code so are liable to go undetected
  # Those that are negative break the hmmloglikelihood.cpp code
  # Those that are negative also trigger error below due to cummin so this check precedes one below
  if(any(fs < 0) | any(fs > 1)) {
    stop("Some frequencies are not in [0,1].", call. = FALSE)
  }
  
  # Check for disordered frequencies.
  # Specifically, check fs ordered s.t. non-zero precede zero for all markers.
  # Disordered frequencies break the hmmloglikelihood.cpp code.
  if(!all(apply(fs01, 1, function(x) all(x == cummin(x))))) {
    stop ("Disordered fs. Per row, all non-zero frequencies should precede all zero frequencies.",
          call. = FALSE)
  }
  
  # Check per-marker frequencies sum to one
  # Those that do not do not break the hmmloglikelihood.cpp code so are liable to go undetected
  fs_sum <- rowSums(fs)
  if(any(fs_sum != 1)) {
    max_dev <- max(abs(fs_sum - 1))
    max_dev
    if(max_dev > max_dev_lim) {
      stop(sprintf("Some markers have frequencies whose sum deviates from one by up to %s.", max_dev),
           call. = FALSE)
    } else {
      warning(sprintf("Some markers have frequencies whose sum deviates from one by up to %s.", max_dev),
              call. = FALSE)
    }
  }
  
  
  
  # Check for uninformative markers
  # These markers do not break the hmmloglikelihood.cpp code but could be omitted.
  if (any(fs == 1)) {
    warning("Some markers are uninformative (have allele frequencies equal to one).",
            call. = FALSE)
  }
  
  if(do_return) { # For use in Ys check in estimate_r_and_k()
    return(list(fs01 = fs01,
                non_zero_fs_lb = non_zero_fs_lb))
  }
}



### estimate_r_and_k----

estimate_r_and_k <- function(fs, ds, Ys, epsilon = 0.001, rho = 7.4 * 10 ^ (-7),
                             kinit = 50, rinit = 0.5, warn_fs = TRUE) {
  
  # Convert to a into matrix if not already (loglikelihood_cpp expects a matrix)
  # and perform some checks
  if (!is.matrix(fs)) fs <- as.matrix(fs)
  # If no errors, returns numeric logic for fs zero/non-zero
  # fs_checks_return <- fs_checks(fs, warn = warn_fs, do_return = TRUE)
  # 
  # # Check Ys for NAs
  # if(any(is.na(Ys))) stop("Missing values detected in Ys.\n  Please remove and recompute ds accordingly.")
  # 
  # # Check for alleles that are permissible only if epsilon exceeds zero.
  # # These alleles will break the hmmloglikelihood.cpp code if epsilon is zero.
  # problem <- "Some per-marker allele counts exceed per-marker non-zero allele frequencies."
  # fs_checks_return$fs01 <- 1 * (fs > fs_checks_return$non_zero_fs_lb)
  # Kts <- rowSums(fs_checks_return$fs01)
  # if (any(sapply(1:nrow(Ys), function(i) Ys[i] > (Kts[i]-1)))) {
  #   if (epsilon > 0) warning (paste0(problem, " Data are permissible due to non-zero epsilon."))
  #   else if (epsilon == 0)  stop(paste0(problem, " Data are incompatible with zero-valued epsilon."))
  # }
  
  # Define the function to pass to optim()
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, fs, ds, epsilon, rho)
  
  # Optimise the negative log likelihood
  optimization <- optim(par = c(kinit, rinit), fn = function(x) - ll(x[1], x[2]))
  
  # Extract and name estimates
  rkhats <- c("khat" = optimization$par[1], "rhat" = optimization$par[2])
  
  if (all(rkhats == c(kinit, rinit))) {
    warning("optimization has returned initial parameter values. Data are possibly uniformative.")
  }
  
  # End of function
  return(rkhats)
}


### pairwise_hmmIBD----

pairwise_hmmIBD = function(ampseq_object, parallel = TRUE, w = 1, n = 1){
  library(parallel)
  library(doMC)
  library(svMisc)
  
  loci_object = ampseq2loci(ampseq_object)
  
  loci_table = loci_object[["loci_table"]]
  freq_table = loci_object[["freq_table"]]
  markers = loci_object[["markers"]]
  
  pairs = as.data.frame(t(combn(rownames(loci_table), 2)))
  
  s = round(seq(1,nrow(pairs)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  pairs = pairs[low:high,]
  
  fx_get_relatedness = function(Yi, Yj, freq_table, markers){
    
    Ys = cbind(Yi, Yj)
    
    rownames(Ys) = colnames(loci_table)
    
    Ys =Ys[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    
    fs = freq_table[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    
    no_na_markers = markers[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    
    no_na_markers[["distance"]] = Inf
    
    for(chromosome in levels(as.factor(no_na_markers[["chromosome"]]))){
      for(amplicon in 1:(nrow(no_na_markers[no_na_markers[["chromosome"]] == chromosome,])-1)){
        no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon, "distance"] = 
          no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon, "pos"]
      }
    }
    
    ds = no_na_markers$distance
    
    estimate = estimate_r_and_k(fs = fs,
                                ds = ds,
                                Ys = Ys,
                                warn_fs = F)
    
    return(estimate)
  }
  
  if(parallel){
    registerDoMC(detectCores())
    pairwise_df = foreach(pair = 1:nrow(pairs), .combine = 'rbind') %dopar% {
      
      # If both samples are polyclonal
      if(sum(grepl("_",loci_table[pairs[pair, 1],])) != 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) != 0){
        Yi = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 1],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 1],])))
        Yj = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 2],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 2],])))
        
        estimate = NULL
        for(haplotype_i in 1:ncol(Yi)){
          for(haplotype_j in 1:ncol(Yj)){
            estimate = rbind(estimate, fx_get_relatedness(Yi[,haplotype_i], Yj[,haplotype_j], freq_table, markers))
          }
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If only the second sample is polyclonal 
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) == 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) != 0){
        Yi = as.numeric(loci_table[pairs[pair, 1],])
        Yj = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 2],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 2],])))
        
        estimate = NULL
        for(haplotype_j in 1:ncol(Yj)){
          estimate = rbind(estimate, fx_get_relatedness(Yi, Yj[,haplotype_j], freq_table, markers))
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If only the first sample is polyclonal
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) != 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) == 0){
        Yi = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 1],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 1],])))
        Yj = as.numeric(loci_table[pairs[pair, 2],])
        
        estimate = NULL
        for(haplotype_i in 1:ncol(Yi)){
          estimate = rbind(estimate, fx_get_relatedness(Yi[,haplotype_i], Yj, freq_table, markers))
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If both samples are monoclonal
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) == 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) == 0){
        
        Yi = as.numeric(loci_table[pairs[pair, 1],])
        Yj = as.numeric(loci_table[pairs[pair, 2],])
        
        estimate = fx_get_relatedness(Yi, Yj, freq_table, markers)
        
      }
      
      data.frame(Yi = pairs[pair, 1], Yj = pairs[pair, 2], t(estimate))
      
    }
    
  }else{
    
    pairwise_df = NULL
    
    for(pair in 1:nrow(pairs)){
      
      # If both samples are polyclonal
      if(sum(grepl("_",loci_table[pairs[pair, 1],])) != 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) != 0){
        Yi = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 1],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 1],])))
        Yj = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 2],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 2],])))
        
        estimate = NULL
        for(haplotype_i in 1:ncol(Yi)){
          for(haplotype_j in 1:ncol(Yj)){
            estimate = rbind(estimate, fx_get_relatedness(Yi[,haplotype_i], Yj[,haplotype_j], freq_table, markers))
          }
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If only the second sample is polyclonal 
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) == 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) != 0){
        Yi = as.numeric(loci_table[pairs[pair, 1],])
        Yj = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 2],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 2],])))
        
        estimate = NULL
        for(haplotype_j in 1:ncol(Yj)){
          estimate = rbind(estimate, fx_get_relatedness(Yi, Yj[,haplotype_j], freq_table, markers))
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If only the first sample is polyclonal
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) != 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) == 0){
        Yi = cbind(as.numeric(gsub("_.+","",loci_table[pairs[pair, 1],])), as.numeric(gsub(".+_","",loci_table[pairs[pair, 1],])))
        Yj = as.numeric(loci_table[pairs[pair, 2],])
        
        estimate = NULL
        for(haplotype_i in 1:ncol(Yi)){
          estimate = rbind(estimate, fx_get_relatedness(Yi[,haplotype_i], Yj, freq_table, markers))
        }
        
        estimate = estimate[which.max(estimate[,'rhat']),]
        
        # If both samples are monoclonal
      }else if(sum(grepl("_",loci_table[pairs[pair, 1],])) == 0 & sum(grepl("_",loci_table[pairs[pair, 2],])) == 0){
        
        Yi = as.numeric(loci_table[pairs[pair, 1],])
        Yj = as.numeric(loci_table[pairs[pair, 2],])
        
        estimate = fx_get_relatedness(Yi, Yj, freq_table, markers)
        
      }
      
      progress(round(100*pair/nrow(pairs)))
      
      pairwise_df = rbind(pairwise_df, data.frame(Yi = pairs[pair, 1], Yj = pairs[pair, 2], t(estimate)))
      
    }
    
  }
  
  return(pairwise_df)
  
}

### plot_relatedness_distribution----

plot_relatedness_distribution = function(pairwise_relatedness = pairwise_relatedness,
                                            metadata = ampseq@metadata,
                                            Population = 'Subnational_level2',
                                            fill_color = c("firebrick3", "firebrick1", "dodgerblue1", "dodgerblue3", "gold1", "gold3"),
                                            type_pop_comparison = 'within', # c('within', 'between', 'both')
                                            ncol = 5,
                                            pop_levels = NULL
){
  
  pairwise_relatedness_l = pairwise_relatedness
  
  pairwise_relatedness_l %<>% filter(Yi %in% metadata$Sample_id,
                                     Yj %in% metadata$Sample_id)
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population)], by.x = 'Yi', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-5], 'Yi_Population')
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population)], by.x = 'Yj', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-6], 'Yj_Population')
  
  pairwise_relatedness_l %<>% filter(!is.na(Yi_Population), !is.na(Yj_Population))
  
  pairwise_relatedness_l %<>% mutate(Pop_comparison = case_when(
    Yi_Population == Yj_Population ~ Yi_Population,
    Yi_Population < Yj_Population ~ paste(Yi_Population, Yj_Population, sep = '_vs_'),
    Yi_Population > Yj_Population ~ paste(Yj_Population, Yi_Population, sep = '_vs_')
  ))
  
  # if(sum(is.na(pairwise_relatedness_l$Pop_comparison)) + sum(grepl('NA', pairwise_relatedness_l$Pop_comparison)) >= 1){
  #   pairwise_relatedness_l %<>% mutate(Pop_comparison = case_when(
  #     is.na(Pop_comparison) | grepl('NA', Pop_comparison) ~ "missing data",
  #     !(is.na(Pop_comparison) | grepl('NA', Pop_comparison)) ~ Pop_comparison
  #   ))
  #   
  #   pairwise_relatedness_l %<>% mutate(Type_of_comparison = case_when(
  #     grepl("-",Pop_comparison) ~ "Between",
  #     !grepl("-",Pop_comparison) ~ "Within",
  #     Pop_comparison == 'missing data' ~ 'missing data'
  #   ))
  #   
  #   pairwise_relatedness_l %<>% filter(Pop_comparison != 'missing data', Type_of_comparison != 'missing data')
  # }else{
  #   
  
  pairwise_relatedness_l %<>% mutate(Type_of_comparison = case_when(
    grepl("_vs_",Pop_comparison) ~ "Between",
    !grepl("_vs_",Pop_comparison) ~ "Within"
  ))
  
  # }
  
  
  if(!is.null(pop_levels)){
    pairwise_relatedness_l$Pop_comparison = factor(pairwise_relatedness_l$Pop_comparison,
                                                   levels = pop_levels)
  }else{
    pairwise_relatedness_l$Pop_comparison = factor(pairwise_relatedness_l$Pop_comparison,
                                                   levels = c(sort(unique(metadata[!is.na(metadata[[Population]]),][[Population]])),
                                                              apply(combn(unique(metadata[!is.na(metadata[[Population]]),][[Population]]),2), 2, function(x){paste(sort(x), collapse = '_vs_')})))
  }
  
  
  pairwise_relatedness_l$Type_of_comparison = factor(pairwise_relatedness_l$Type_of_comparison, levels =
                                                       c('Within', 'Between'))
  
  
  
  if(type_pop_comparison == 'within'){
    plot_pairwise_relatedness_distribution = pairwise_relatedness_l %>%
      filter(Type_of_comparison == 'Within')%>%
      ggplot(aes(x = rhat, fill = Pop_comparison)) +
      geom_histogram(position = "stack", alpha = .7)+
      geom_vline(xintercept = mean(pairwise_relatedness_l$rhat), linetype = 2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      facet_wrap(~Pop_comparison, ncol = ncol,
                 scales = "free_y")+
      labs(y = "Count",
           x = "Relatedness",
           fill = "Populations")+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
    pairwise_relatedness_l %<>% filter(Type_of_comparison == 'Within')
    
  }else if(type_pop_comparison == 'between'){
    plot_pairwise_relatedness_distribution = pairwise_relatedness_l %>%
      filter(Type_of_comparison == 'Between')%>%
      ggplot(aes(x = rhat, fill = Pop_comparison)) +
      geom_histogram(position = "stack", alpha = .7)+
      geom_vline(xintercept = mean(pairwise_relatedness_l$rhat), linetype = 2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      facet_wrap(~Pop_comparison, ncol = ncol,
                 scales = "free_y")+
      labs(y = "Count",
           x = "Relatedness",
           fill = "Populations")+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
    pairwise_relatedness_l %<>% filter(Type_of_comparison == 'Between')
    
  }else if(type_pop_comparison == 'both'){
    plot_pairwise_relatedness_distribution = pairwise_relatedness_l %>%
      ggplot(aes(x = rhat, fill = Pop_comparison)) +
      geom_histogram(position = "stack", alpha = .7)+
      geom_vline(xintercept = mean(pairwise_relatedness_l$rhat), linetype = 2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      facet_wrap(Type_of_comparison~
                   Pop_comparison,
                 scales = "free_y")+
      labs(y = "Count",
           x = "Relatedness",
           fill = "Populations")+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.position = "none")
  }
  
  
  pairwise_relatedness_distribution = list(relatedness = pairwise_relatedness_l,
                                           plot = plot_pairwise_relatedness_distribution)
  
  return(pairwise_relatedness_distribution)
  
}

### plot_frac_highly_related----

plot_frac_highly_related = function(pairwise_relatedness = pairwise_relatedness,
                                       metadata = ampseq@metadata,
                                       Population = 'Population',
                                       fill_color = c("dodgerblue3",  "firebrick3", "gold3", "gray50", "gray50", "gray50"),
                                       threshold = 0.99,
                                       type_pop_comparison = 'between',
                                       pop_levels = NULL){
  
  pairwise_relatedness_l = pairwise_relatedness
  
  pairwise_relatedness_l %<>% filter(Yi %in% metadata$Sample_id,
                                     Yj %in% metadata$Sample_id)
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population)], by.x = 'Yi', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-5], 'Yi_Population')
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population)], by.x = 'Yj', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-6], 'Yj_Population')
  

  pairwise_relatedness_l %<>% filter(!is.na(Yi_Population), !is.na(Yj_Population))
  
  pairwise_relatedness_l %<>% mutate(Pop_comparison = case_when(
    Yi_Population == Yj_Population ~ Yi_Population,
    Yi_Population < Yj_Population ~ paste(Yi_Population, Yj_Population, sep = '_vs_'),
    Yi_Population > Yj_Population ~ paste(Yj_Population, Yi_Population, sep = '_vs_')
  ))
  
  # if(sum(is.na(pairwise_relatedness_l$Pop_comparison)) + sum(grepl('NA', pairwise_relatedness_l$Pop_comparison)) >= 1){
  #   pairwise_relatedness_l %<>% mutate(Pop_comparison = case_when(
  #     is.na(Pop_comparison) | grepl('NA', Pop_comparison) ~ "missing data",
  #     !(is.na(Pop_comparison) | grepl('NA', Pop_comparison)) ~ Pop_comparison
  #   ))
  #   
  #   pairwise_relatedness_l %<>% mutate(Type_of_comparison = case_when(
  #     grepl("-",Pop_comparison) ~ "Between",
  #     !grepl("-",Pop_comparison) ~ "Within",
  #     Pop_comparison == 'missing data' ~ 'missing data'
  #   ))
  #   
  #   pairwise_relatedness_l %<>% filter(Pop_comparison != 'missing data', Type_of_comparison != 'missing data')
  # }else{
  #   
  pairwise_relatedness_l %<>% mutate(Type_of_comparison = case_when(
    grepl("_vs_",Pop_comparison) ~ "Between",
    !grepl("_vs_",Pop_comparison) ~ "Within"
  ))
  
  # }
  
  if(!is.null(pop_levels)){
    pairwise_relatedness_l$Pop_comparison = factor(pairwise_relatedness_l$Pop_comparison,
                                                   levels = pop_levels)
    
  }else{
    pairwise_relatedness_l$Pop_comparison = factor(pairwise_relatedness_l$Pop_comparison,
                                                   levels = c(sort(unique(metadata[!is.na(metadata[[Population]]),][[Population]])),
                                                              apply(combn(unique(metadata[!is.na(metadata[[Population]]),][[Population]]),2), 2, function(x){paste(sort(x), collapse = '_vs_')})))
    
  }
  
   
  
  pairwise_relatedness_l$Type_of_comparison = factor(pairwise_relatedness_l$Type_of_comparison, levels =
                                                       c('Within', 'Between'))
  
  
  highly_related_table = pairwise_relatedness_l %>%
    group_by(Type_of_comparison, Pop_comparison) %>% 
    dplyr::summarise(freq = sum(rhat >= threshold),
                     n = n()) %>% group_by(Type_of_comparison, Pop_comparison)%>%
    mutate(prop = binconf(freq,
                          n,
                          alpha = 0.05,
                          method = "exact")[1],
           lower = binconf(freq,
                           n,
                           alpha = 0.05,
                           method = "exact")[2],
           upper = binconf(freq,
                           n,
                           alpha = 0.05,
                           method = "exact")[3]) %>%
    mutate(lower = case_when(
      prop == 0 ~ 0,
      prop != 0 ~ lower),
      upper = case_when(
        prop == 0 ~ 0,
        prop != 0 ~ upper))
  
  if(type_pop_comparison == 'within'){
    
    plot_frac_highly_related = highly_related_table %>%
      filter(Type_of_comparison == 'Within') %>%
      ggplot(aes(x = Pop_comparison, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      labs(y = paste0('Frac. of samples w. IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 90, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
    highly_related_table %<>%
      filter(Type_of_comparison == 'Within')
    
  }else if(type_pop_comparison == 'between'){
    
    plot_frac_highly_related = highly_related_table %>%
      filter(Type_of_comparison == 'Between') %>%
      ggplot(aes(x = Pop_comparison, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      labs(y = paste0('Frac. of samples w. IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 90, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
    highly_related_table %<>%
      filter(Type_of_comparison == 'Between')
    
  }else if(type_pop_comparison == 'both'){
    
    plot_frac_highly_related = highly_related_table %>%
      ggplot(aes(x = Pop_comparison, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      labs(y = paste0('Frac. of samples w. IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 90, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }
  
  frac_highly_related = list(highly_related_table = highly_related_table,
                             plot = plot_frac_highly_related)
  
  return(frac_highly_related)
  
}

### plot_frac_highly_related_over_time----

plot_frac_highly_related_over_time = function(pairwise_relatedness = pairwise_relatedness,
                                                 metadata = ampseq@metadata,
                                                 Population = c('Population', 'Quarter_of_Collection'),
                                                 fill_color = c("dodgerblue3",  "firebrick3", "gold3", "gray50", "gray50", "gray50"),
                                                 threshold = 0.99,
                                                 type_pop_comparison = 'within',## c('within', 'between', 'both')
                                                 ncol = 4,
                                                 pop_levels = NULL){
  
  pairwise_relatedness_l = pairwise_relatedness
  
  pairwise_relatedness_l %<>% filter(Yi %in% metadata$Sample_id,
                                     Yj %in% metadata$Sample_id)
  
  # Giving wraning message, otherwise it is faster
  # pairwise_relatedness_l %<>% mutate(
  #   Population_Yi = metadata[grepl(Yi, metadata[['samples']]),][[Population[1]]],
  #   Population_Yj = metadata[grepl(Yj, metadata[['samples']]),][[Population[1]]],
  #   Date_Yi = metadata[grepl(Yi, metadata[['samples']]),][[Population[2]]],
  #   Date_Yj = metadata[grepl(Yj, metadata[['samples']]),][[Population[2]]])
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population[1])], by.x = 'Yi', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-5], 'Population_Yi')
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population[1])], by.x = 'Yj', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-6], 'Population_Yj')
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population[2])], by.x = 'Yi', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-7], 'Date_Yi')
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population[2])], by.x = 'Yj', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-8], 'Date_Yj')
  
  pairwise_relatedness_l %<>% filter(!is.na(Population_Yi),
                                     !is.na(Population_Yj),
                                     !is.na(Date_Yi),
                                     !grepl('NA',Date_Yi),
                                     !is.na(Date_Yj),
                                     !grepl('NA',Date_Yj))
  
  pairwise_relatedness_l %<>% mutate(Pop_Date_Yi = paste(Population_Yi, Date_Yi, sep = "_"),
                                     Pop_Date_Yj = paste(Population_Yj, Date_Yj, sep = "_"))
  
  
  pairwise_relatedness_l_sorted = NULL
  
  comparisons = cbind(rbind(
    sort(unique(c(unique(pairwise_relatedness_l$Pop_Date_Yi), unique(pairwise_relatedness_l$Pop_Date_Yj)))),
    sort(unique(c(unique(pairwise_relatedness_l$Pop_Date_Yi), unique(pairwise_relatedness_l$Pop_Date_Yj))))),
    combn(sort(unique(c(unique(pairwise_relatedness_l$Pop_Date_Yi), unique(pairwise_relatedness_l$Pop_Date_Yj)))), 2))

  
  
  for(comparison in 1:ncol(comparisons)){
    Pop_Date_Yi = comparisons[,comparison][1]
    Pop_Date_Yj = comparisons[,comparison][2]
    if(Pop_Date_Yi == Pop_Date_Yj){
      
      pairwise_relatedness_l_sorted = rbind(
        pairwise_relatedness_l_sorted,
        pairwise_relatedness_l[pairwise_relatedness_l[['Pop_Date_Yi']] == Pop_Date_Yi &
                               pairwise_relatedness_l[['Pop_Date_Yj']] == Pop_Date_Yj, ])
      
    }else{
      
      temp1 = pairwise_relatedness_l[pairwise_relatedness_l[['Pop_Date_Yi']] == Pop_Date_Yi &
                                       pairwise_relatedness_l[['Pop_Date_Yj']] == Pop_Date_Yj, ]
      
      temp2 = pairwise_relatedness_l[pairwise_relatedness_l[['Pop_Date_Yi']] == Pop_Date_Yj &
                                       pairwise_relatedness_l[['Pop_Date_Yj']] == Pop_Date_Yi, ]
      
      Yi = temp2[['Yj']]
      Population_Yi = temp2[['Population_Yj']]
      Date_Yi = temp2[['Date_Yj']]
      Pop_Date_Yi = temp2[['Pop_Date_Yj']]
      
      Yj = temp2[['Yi']]
      Population_Yj = temp2[['Population_Yi']]
      Date_Yj = temp2[['Date_Yi']]
      Pop_Date_Yj = temp2[['Pop_Date_Yi']]
      
      temp2[['Yi']] = Yi
      temp2[['Population_Yi']] = Population_Yi
      temp2[['Date_Yi']] = Date_Yi
      temp2[['Pop_Date_Yi']] = Pop_Date_Yi
      
      temp2[['Yj']] = Yj
      temp2[['Population_Yj']] = Population_Yj
      temp2[['Date_Yj']] = Date_Yj
      temp2[['Pop_Date_Yj']] = Pop_Date_Yj
      
      pairwise_relatedness_l_sorted = rbind(pairwise_relatedness_l_sorted,
                                            temp1,
                                            temp2)
      
    }
    
  }
    
  plot_IBD_correlation_matrix = pairwise_relatedness_l_sorted %>%
    group_by(Pop_Date_Yi, Pop_Date_Yj) %>% 
    dplyr::summarise(freq = sum(rhat >= threshold),
                     n = n()) %>%
    group_by(Pop_Date_Yi, Pop_Date_Yj) %>% 
    mutate(prop = binconf(freq,
                          n,
                          alpha = 0.05,
                          method = "exact")[1],
           lower = binconf(freq,
                           n,
                           alpha = 0.05,
                           method = "exact")[2],
           upper = binconf(freq,
                           n,
                           alpha = 0.05,
                           method = "exact")[3])%>%
    ggplot(aes(x = Pop_Date_Yi, y = Pop_Date_Yj, fill = prop)) + 
    geom_tile(alpha = .85)+
    scale_fill_gradient2(low = "white", high = "red", space = "Lab",
                         name=paste0('IBD >= ', threshold))+
    theme_bw()+
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, vjust = .5),
          axis.title.y = element_text(size = 12),
          strip.text = element_text(size = 12))
  
  
  pairwise_relatedness_l_sorted %<>% mutate(Pop_comparison = case_when(
    Population_Yi < Population_Yj ~ paste(Population_Yi, Population_Yj, sep = '_vs_'),
    Population_Yi > Population_Yj ~ paste(Population_Yj, Population_Yi, sep = '_vs_'),
    Population_Yi == Population_Yj ~ Population_Yi),
    Type_Pop_comparison = case_when(
      Population_Yi != Population_Yj ~ "Between",
      Population_Yi == Population_Yj ~ "Within"),
    Date_comparison = case_when(
      Date_Yi != Date_Yj ~ "Between",
      Date_Yi == Date_Yj ~ "Within")
  )
  
  if(!is.null(pop_levels)){
    pairwise_relatedness_l_sorted$Pop_comparison = factor(pairwise_relatedness_l_sorted$Pop_comparison,
                                                   levels = pop_levels)
    
  }else{
    pairwise_relatedness_l_sorted$Pop_comparison = factor(pairwise_relatedness_l_sorted$Pop_comparison,
                                                   levels = c(sort(unique(pairwise_relatedness_l_sorted$Pop_comparison)[!grepl("_vs_", unique(pairwise_relatedness_l_sorted$Pop_comparison))]),
                                                              sort(unique(pairwise_relatedness_l_sorted$Pop_comparison)[grepl("_vs_", unique(pairwise_relatedness_l_sorted$Pop_comparison))])))
  }
  
    
  pairwise_relatedness_l_sorted$Type_Pop_comparison = factor(pairwise_relatedness_l_sorted$Type_Pop_comparison, levels =
                                                        c('Within', 'Between'))
  
  pairwise_relatedness_l_sorted %<>% filter(Date_comparison == 'Within')
  
  frac_highly_related = pairwise_relatedness_l_sorted %>%
    group_by(Pop_comparison, Date_Yi, Type_Pop_comparison) %>% 
    dplyr::summarise(freq = sum(rhat >= threshold),
                     n = n()) %>% group_by(Pop_comparison, Date_Yi, Type_Pop_comparison)%>%
    mutate(prop = binconf(freq,
                          n,
                          alpha = 0.05,
                          method = "exact")[1],
           lower = binconf(freq,
                           n,
                           alpha = 0.05,
                           method = "exact")[2],
           upper = binconf(freq,
                           n,
                           alpha = 0.05,
                           method = "exact")[3])%>%
    mutate(lower = case_when(
      prop == 0 ~ 0,
      prop != 0 ~ lower),
      upper = case_when(
        prop == 0 ~ 0,
        prop != 0 ~ upper))
  
  if(type_pop_comparison == 'both'){
    
    plot_frac_highly_related = frac_highly_related %>%
      ggplot(aes(x = Date_Yi, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      facet_wrap(Type_Pop_comparison~Pop_comparison)+
      theme_bw()+
      labs(y = paste0('Frac. of samples w. IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 90, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }else if(type_pop_comparison == 'within'){
    
    plot_frac_highly_related = frac_highly_related %>%
      filter(Type_Pop_comparison == 'Within')%>%
      ggplot(aes(x = Date_Yi, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      facet_wrap(~Pop_comparison, ncol = ncol)+
      theme_bw()+
      labs(y = paste0('Frac. of samples w. IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 90, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
    frac_highly_related %<>%
      filter(Type_Pop_comparison == 'Within')
    
  }else if(type_pop_comparison == 'between'){
    
    plot_frac_highly_related = frac_highly_related %>%
      filter(Type_Pop_comparison == 'Between')%>%
      ggplot(aes(x = Date_Yi, y = prop, fill = Pop_comparison)) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      facet_wrap(~Pop_comparison, ncol = ncol)+
      theme_bw()+
      labs(y = paste0('Frac. of samples w. IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 90, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
    frac_highly_related %<>%
      filter(Type_Pop_comparison == 'Between')
    
  }
  
  frac_highly_related = list(frac_highly_related = frac_highly_related,
                             plot_frac_highly_related = plot_frac_highly_related,
                             plot_IBD_correlation_matrix = plot_IBD_correlation_matrix)
  
  return(frac_highly_related)
  
}

### plot_network----

plot_network = function(pairwise_relatedness,  
                           threshold,
                           metadata,
                           sample_id,
                           group_by,
                           levels,
                           colors,
                        vertex.size = 4){
  
  if(sum(is.na(metadata[[group_by]])) > 0){
    metadata[is.na(metadata[[group_by]]),][[group_by]] = 'missing data'
  }
  
  pairwise_relatedness_l = pairwise_relatedness
  
  pairwise_relatedness_l %<>% filter(Yi %in% metadata$Sample_id,
                                Yj %in% metadata$Sample_id)
  
  variable = 'rhat'
  cols = c("Yi", "Yj")
  
  library(igraph)
  
  edges = NULL
  
  selected_edges <- pairwise_relatedness_l[pairwise_relatedness_l[[variable]] >= threshold,][,cols]
  for(edge in 1:nrow(selected_edges)){
    edges = c(edges,unlist(selected_edges[edge,]))
  }
  
  isolates = unique(unlist(pairwise_relatedness_l[,cols]))[!(unique(unlist(pairwise_relatedness_l[,cols])) %in% unique(edges))]
  
  network_object = graph(edges = edges, isolates = isolates, directed=F)
  
  node_colors = c(unique(edges), isolates)
  
  names(node_colors) = node_colors
  
  for(level in 1:length(levels)){
    node_colors[node_colors %in% metadata[metadata[[group_by]] == levels[level],][[sample_id]]] = colors[level]
  }
  
  plot_network = plot.igraph(network_object,
                             vertex.color = node_colors,
                             vertex.size=vertex.size,
                             vertex.label.cex=0,
                             vertex.label.dist=.5, 
                             #vertex.label.degree=-pi/2,
                             vertex.label.color="black",
                             vertex.label=NA,
                             asp=0,
                             xlim=c(-1,1),
                             ylim=c(-1,1))
  
  return(list(network_object = network_object,
              plot_network = plot_network))
  
}



### fastGRMcpp----
#' C++ implementation of a Genomic relationship matrix 'GRM'
#' 
#' @param X Matrix of the type 'MatrixXd' for which the GRM will be calculated.
#' 
#' @return Genomic relationship matrix (GRM).
#' 
#' @importFrom Rdpack reprompt
#' @references https://doi.org/10.1016/j.ajhg.2010.11.011
#' 
#' @examples
#' require(fastGRM)
#' Data = matrix(sample(0:1, 9000, TRUE, c(.9,.1)), 90)
#' X = grm(Data)
#' 
#' @export
#' 

grm = function(X){
  grmCpp(X)
}


#' C++ implementation of a fast singular value decomposition (SVD)
#' 
#' @param X Symmetric matrix of the type 'MatrixXd' for which the SVD will be calculated.
#' @param k Number of first k eigen vectors to return
#' @param q Auxiliary exponent
#' 
#' @return SVD matrix of size .
#' 
#' @importFrom Rdpack reprompt
#' @references https://doi.org/10.48550/arXiv.0909.4061
#' 
#' @examples
#' require(fastGRM)
#' Data = matrix(sample(0:1, 9000, TRUE, c(.9,.1)), 90)
#' X = grm(Data)
#' V = fastSVD(X, 2)
#' 
#' @export

fastSVD = function(X, k, q = 2){
  
  
  
  fastSVDCpp(X, k, q)
}

#' C++ implementation of a fast GRM function
#' 
#' @param X Matrix of the type 'MatrixXd' for which the fastGRM matrix will be calculated.
#' @param k Number of first k eigen vectors to return
#' @param q Auxiliary exponent
#' 
#' @return SVD matrix of size .
#' 
#' @importFrom Rdpack reprompt
#' @references https://doi.org/10.1016/j.ajhg.2010.11.011
#' @references https://doi.org/10.48550/arXiv.0909.4061
#' 
#' @examples
#' require(fastGRM)
#' Data = matrix(sample(0:1, 9000, TRUE, c(.9,.1)), 90)
#' X = fastGRM(Data, 2)
#' 
#' @export
#' 

GRM_evectors = function(gt, k = nrow(gt), metadata, Pop = 'Population', q = 2){
  
  X = t(gsub('_.+$','',gt))
  
  X = matrix(as.numeric(X), ncol = ncol(X),
             nrow = nrow(X), 
             dimnames = list(
               rownames(X),
               colnames(X)
             ))
  
  X[is.na(X)] = 0
  
  evector = fastGRMCpp(X, k, q)
  
  #### Add metadata to the PCA
  Pop_col = merge(data.frame(Sample_id = gsub('_C[1,2]$','',colnames(X)),
                             order = 1:ncol(X)), metadata[,c('Sample_id', Pop)], by = 'Sample_id', all.x = T)
  
  Pop_col = Pop_col[order(Pop_col$order),]
  
  evector = data.frame(Pop_col, evector)
  names(evector) = c(colnames(Pop_col), paste0(rep('PC', k), 1:k))
  
  return(evector)
  
}

### IBD_evectors----
IBD_evectors = function(ampseq_object, relatedness_table, k = NULL, Pop = 'Population', q = 2){
  
  metadata = ampseq_object@metadata
  
  
  
  pairwise_relatedness_matrix = matrix(data = NA,
                                       ncol = nrow(ampseq_object@metadata),
                                       nrow = nrow(ampseq_object@metadata),
                                       dimnames = list(ampseq_object@metadata$Sample_id,
                                                       ampseq_object@metadata$Sample_id))
  
  #### Add metadata to the PCA
  Pop_col = merge(data.frame(Sample_id = gsub('_C[1,2]$','',colnames(pairwise_relatedness_matrix)),
                             order = 1:ncol(pairwise_relatedness_matrix)), metadata[,c('Sample_id', Pop)], by = 'Sample_id', all.y = T)
  
  Pop_col = Pop_col[order(Pop_col$order),]
  
  
  relatedness_table %<>% filter(Yi %in% metadata$Sample_id,
                                Yj %in% metadata$Sample_id)
  
  
  for(pair in 1:nrow(relatedness_table)){
    
    pairwise_relatedness_matrix[relatedness_table[pair,][['Yi']],
                                relatedness_table[pair,][['Yj']]] = 
      relatedness_table[pair,][['rhat']]
    
    pairwise_relatedness_matrix[relatedness_table[pair,][['Yj']],
                                relatedness_table[pair,][['Yi']]] = 
      relatedness_table[pair,][['rhat']]
  }
  
  pairwise_relatedness_matrix[is.na(pairwise_relatedness_matrix)] = 1
  
  pairwise_dist_matrix = 1 - pairwise_relatedness_matrix
  
  ## Using_fastSVDCpp
  
  # evector = fastSVDCpp(pairwise_relatedness_matrix, k, q)
  # 
  # for(i in 1:k){
  #   evector[,i] = sign(evector[1,i])*evector[,i]
  # }
  # 
  # evalues = NULL
  # 
  # for(i in 1:k){
  #   evalues = c(evalues, unlist((pairwise_dist_matrix %*% evector[,i])/evector[,i])[1])
  # }
  # 
  # contrib = 100*(evalues)^2/sum((evalues)^2)
  # 
  # evector = data.frame(Pop_col, evector)
  # names(evector) = c(colnames(Pop_col), paste0(rep('PC', k), 1:k))
  # 
  # ibd_pca_fast = list(eigenvector=evector, eigenvalues = evalues, contrib = contrib)
  
  ## Using princomp from R
  
  ibd_pca = princomp(pairwise_dist_matrix, cor = T)
  
  ibd_evector = ibd_pca$scores
  ibd_evalues = ibd_pca$sdev
  
  ibd_contrib = 100*(ibd_evalues)^2/sum((ibd_evalues)^2)
  
  ibd_evector = data.frame(Pop_col, ibd_evector)
  names(ibd_evector) = c(colnames(Pop_col), paste0(rep('PC', k), 1:k))
  
  ibd_pca = list(eigenvector=ibd_evector, eigenvalues = ibd_evalues, contrib = ibd_contrib)
  
  return(ibd_pca)
  
}

## Functions for genetic diversity----

### get_locus_diversity----
get_locus_diversity = function(locus, variance = TRUE){
  
  alleles = levels(as.factor(unlist(strsplit(gsub(":[0-9]+", "", locus), "_"))))
  
  if(length(alleles) > 0){
    
    n.all = length(alleles)
    
    n = sum(!is.na(unlist(strsplit(locus, "_"))))
    
    freq = sort(sapply(alleles, function(allele){
      sum(grepl(paste("(^|_)", allele, ":", sep = ""), locus))/n
    }), decreasing = T)
    
    sp2 <- sum(freq^2)
    H <- n * (1 - sp2)/(n - 1)
    na.e <- 1/(1-H)
    if (variance) {
      sp3 <- sum(freq^3)
      var.H <- 2 * (2 * (n - 2) * (sp3 - sp2^2) + sp2 - sp2^2) / (n * (n - 1))
      Nei.var.H <- (2*(n-1)/n^3)*((3-2*n)*(sp2^2)+2*(n-2)*sp3+sp2)
      return(c(n.all = n.all, na.e = na.e, Hexp = H, var.Hexp = var.H, Nei.var.Hexp = Nei.var.H))
    }
    else return(c(n.all = n.all, na.e = na.e, Hexp = H))
    
  }else{
    
    if (variance) {
      return(c(n.all = NA, na.e = NA, Hexp = NA, var.Hexp = NA, Nei.var.Hexp = NA))
    }
    else return(c(n.all = NA, na.e = NA, Hexp = NA))
    
  }
  
}


### get_loci_diversity----
get_loci_diversity = function(loci_abd_table = NULL, variance = TRUE){
  
  #source("fx_locus_diversity.R")
  
  loci_diversity = data.frame(t(sapply(colnames(loci_abd_table), function(locus) get_locus_diversity(loci_abd_table[,locus], variance = variance))))
  
  return(loci_diversity)
  
  
}

### get_pop_diversity----
get_pop_diversity = function(ampseq_object, strata){
  
  # source("fx_loci_diversity.R")
  
  ampseq_loci_abd_table = ampseq_object@gt
  metadata = ampseq_object@metadata
  
  mlg = as.factor(apply(gsub(":[0-9]+", "",ampseq_loci_abd_table[metadata[["max_nAlleles"]] == 1, ]), 1, function(sample) paste(sample, collapse = "_")))
  
  n_mlg = nlevels(mlg) # Richness: number of different variants (species, genus, families, ASVs)
  p = summary(mlg)/length(mlg) # Frequency of each variant
  sp2 = sum(p^2)
  D = length(mlg) * (1 - sp2)/(length(mlg) - 1) # Simpson: probability of choosing two different variants
  S.e = 1/(1-D) # effective richness: The number of variants that explain the diversity of Simpson
  H = -sum(p*log(p)) # Shannon
  E = H/log(n_mlg) # Evenness: from 0 to 1
  
  loci_diversity = get_loci_diversity(ampseq_loci_abd_table, variance = F)
  
  pop_diversity = data.frame(pop = "Total",
                             Richness = n_mlg,
                             Effective.richness = S.e,
                             Simpson = D,
                             Shannon = H,
                             Evenness = E,
                             n.all = mean(loci_diversity$n.all),
                             na.e = mean(loci_diversity$na.e),
                             Hexp = mean(loci_diversity$Hexp)
  )
  
  for(pop in levels(as.factor(metadata[[strata]]))){
    
    if(is.null(nrow(ampseq_loci_abd_table[metadata[[strata]] == pop & metadata[["max_nAlleles"]] == 1, ]))){
      
      pop_diversity = rbind(pop_diversity,data.frame(pop = pop,
                                                     Richness = 1,
                                                     Effective.richness = 1,
                                                     Simpson = 0,
                                                     Shannon = -sum(1*log(1)),
                                                     Evenness = -sum(1*log(1))/log(1),
                                                     n.all = 1,
                                                     na.e = 1,
                                                     Hexp = 0))
    }else{
      mlg = as.factor(apply(gsub(":[0-9]+", "",ampseq_loci_abd_table[metadata[[strata]] == pop & metadata[["max_nAlleles"]] == 1, ]), 1, function(sample) paste(sample, collapse = "_")))
      n_mlg = nlevels(mlg) # Richness: number of different variants (species, genus, families, ASVs)
      p = summary(mlg)/length(mlg)# Frequency of each variant
      sp2 = sum(p^2)
      D = length(mlg) * (1 - sp2)/(length(mlg) - 1) # Simpson: probability of choosing two different variants
      S.e = 1/(1-D) # effective richness: The number of variants that explain the diversity of Simpson
      H = -sum(p*log(p)) # Shannon
      E = H/log(n_mlg) # Evenness: from 0 to 1
      
      loci_diversity = get_loci_diversity(loci_abd_table = ampseq_loci_abd_table[metadata[[strata]] == pop, ], variance = F)
      
      pop_diversity = rbind(pop_diversity,data.frame(pop = pop,
                                                     Richness = n_mlg,
                                                     Effective.richness = S.e,
                                                     Simpson = D,
                                                     Shannon = H,
                                                     Evenness = E,
                                                     n.all = mean(loci_diversity$n.all, na.rm = T),
                                                     na.e = mean(loci_diversity$na.e, na.rm = T),
                                                     Hexp = mean(loci_diversity$Hexp, na.rm = T)))
    }
    
  }
  
  
  return(pop_diversity)
  
}


# Other functions----
## log_scale_histogram----

log_scale_histogram = function(data, var, binwidth, group_by, levels, x_label, fill_color, y_breaks, ncol,
                               na.rm = FALSE,
                               filters = NULL){
  
  data = data[, c(var, group_by)]
  
  names(data) = c("variable", "groups")
  
  if(na.rm){
    data = data[!is.na(data[['groups']]) & !grepl('NA', data[['groups']]),]
  }else if(length(data[is.na(data[['groups']]) | grepl('NA', data[['groups']]),][['groups']])>0){
    data[is.na(data[['groups']]) | grepl('NA', data[['groups']]),][['groups']] = 'missing data'
  }
  
  if(!is.null(filters)){
    data = data[data[['groups']] %in% filters,]
  }
  
  bins = seq(min(data[["variable"]]), max(data[["variable"]]) + binwidth, binwidth)
  
  bins = data.frame(start = bins[-length(bins)], end = bins[-1])
  
  count_table = NULL
  
  for(group in levels(as.factor(data[["groups"]]))){
    temp_data = data[data[["groups"]]==group,]
    
    counts = sapply(1:nrow(temp_data), function(x) {
      which(bins[,1] <= temp_data[x,][["variable"]] & bins[,2] > temp_data[x,][["variable"]] )
    })
    
    for(bin in 1:nrow(bins)){
      count_table = rbind(count_table,
                          data.frame(group = group,
                                     bin = bins[bin,1],
                                     Count = length(counts[counts == bin])))
    }
  }
  
  
  plot_histogram = count_table %>% ggplot(aes(x = bin, y = log10(Count + 1), fill = factor(group, levels = levels)))+
    geom_col(alpha = .7)+
    facet_wrap(.~factor(group,
                        levels = levels), ncol = ncol)+
    scale_y_continuous(breaks = log10(y_breaks + 1), labels = y_breaks)+
    labs(y = "Count",
         x = x_label)+
    scale_fill_manual(values = fill_color)+
    theme_bw()+
    theme(legend.position = "none")
  
  return(plot_histogram)
  
}


## gadm_loadtCountries----
## Function : gadm_loadtCountries (constructor)
## Description : load a file from local system or from GADM repository
##               You just have to specify the countries (ISO3 CODE) of the
##               file name, like "ARG" for Argentina.
##               Optionally you can specify which level you want to have and
##               simplify or not theshapefile (a value less or equal 0.01
##               is recommended)
## Return : This function creates a gadm_sp that contains a
##          SpatialPolygonsDataFrame object that contains all maps you
##          specify in "fileNames".

gadm_sp_loadCountries <- function (fileNames,
                                   level = 0,
                                   basefile=GADM_BASE,
                                   baseurl=GADM_URL,
                                   simplify=NULL)
{
  loadNamespace("sp")
  
  # Load file and change Prefix 
  loadChangePrefix <- function (fileName, level = 0) {
    FILENAME = sprintf("%s_adm%d.rds", fileName,level)
    LOCAL_FILE = sprintf("%s%s", basefile, FILENAME)
    if (!file.exists(LOCAL_FILE)) {
      .OS <- toupper(Sys.info()["sysname"])
      REMOTEFILE = sprintf("gadm36_%s_%d_sp.rds", fileName,level)
      REMOTE_LINK <- sprintf("%s%s", baseurl, REMOTEFILE)
      if (.OS == "WINDOWS") {
        download.file(REMOTE_LINK, LOCAL_FILE, method="wininet",mode="wb")
      } else {
        download.file(REMOTE_LINK, LOCAL_FILE, method = 'auto')
      }
    }
    gadm <- readRDS(LOCAL_FILE)
    if (!is.null(gadm)) {
      theFile <- sp::spChFIDs(gadm, paste(fileName, row.names(gadm), sep = "_"))
      theFile
    }
  }
  
  polygon <- sapply(fileNames, loadChangePrefix, level)
  polyMap <- do.call("rbind", polygon)
  # ---- Simplify polygones if requested by user
  if (!is.null(simplify)) {
    S <- gSimplify(polyMap, simplify, topologyPreserve = TRUE)
    polyMap@polygons <- S@polygons
  }
  
  # ---- Create gadm_sp object
  structure(list("basename"=basefile,
                 "spdf"=polyMap,
                 "level"=level,
                 "holes" = c(),
                 "L360" = FALSE,
                 "stripped" = FALSE,
                 "hasBGND"  = FALSE),
            class = "gadm_sp")
}


## nthroot----
nthroot = function(x,n) {
  (abs(x)^(1/n))*sign(x)
}



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
        
        if(sum(grepl('(/|-|:)', colnames(cigar_run))) > 0){
          
          colnames(cigar_run) = gsub('(/|-|:)', '_', colnames(cigar_run))
          print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the cigar table.')
          
        }
        
      }else if(file.exists(file.path(paths, Run, "dada2/run_dada2/CIGARVariants_Bfilter.out.tsv"))){
        cigar_run = read.table(file.path(paths, Run, "dada2/run_dada2/CIGARVariants_Bfilter.out.tsv"), header = T, check.names = FALSE)
        
        if(sum(grepl('(/|-|:)', colnames(cigar_run))) > 0){
          
          colnames(cigar_run) = gsub('(/|-|:)', '_', colnames(cigar_run))
          print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the cigar table.')
          
        }
        
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
        
        
        if(sum(grepl('(/|-|:)', asv_table_run[['Amplicon']])) > 0){
          
          asv_table_run[['Amplicon']] = gsub('(/|-|:)', '_', asv_table_run[['Amplicon']])
          print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the column Amplicon of the ASV table.')
          
        }
        
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
        print(paste0('ZeroReadSamples file ', file.path(paths, Run, "dada2/run_dada2/zeroReadSamples.txt"), ' not found'))
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
    print("Iterating through cigar_files...")
    
    
    for(file in 1:length(cigar_files)){
      cigar_run = read.table(cigar_files[file], header = T, check.names = FALSE
      )
      
      if(sum(grepl('(/|-|:)', colnames(cigar_run))) > 0){
        
        colnames(cigar_run) = gsub('(/|-|:)', '_', colnames(cigar_run))
        print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the cigar table.')
        
      }
      
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
          
          if(sum(grepl('(/|-|:)', asv_table_run[['Amplicon']])) > 0){
            
            asv_table_run[['Amplicon']] = gsub('(/|-|:)', '_', asv_table_run[['Amplicon']])
            print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the column Amplicon of the ASV table.')
            
          }
          
          # Saving asv_table in cigar_tables object
          cigar_tables[[cigar_files[file]]][['asv_table']] = asv_table_run
          
          # Uploading asv_seqs
          asv_seqs = Biostrings::readDNAStringSet(asv_seq_files[file])
          
          # Saving asv_seqs in cigar_tables object
          cigar_tables[[cigar_files[file]]][['asv_seqs']] = asv_seqs
          
          # adding ZeroReadSamples
          
          ZeroReadSamples = read.table(zero_read_sample_list[file], header = T)
          
          if(nrow(ZeroReadSamples) > 0){
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
          }
          
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
      
      # ASVs in the new data set
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
        unconsitent_cigar_strings = data.frame(temp_asv_table2[temp_asv_table2$CIGAR != temp_asv_table3$CIGAR,], 
                                               CIGAR2 = temp_asv_table3[temp_asv_table3$CIGAR != temp_asv_table2$CIGAR,][['CIGAR']])
        
        for(allele in unconsitent_cigar_strings$CIGAR){
          
          cigar_table[
            cigar_table$alleles == 
              paste(unconsitent_cigar_strings[unconsitent_cigar_strings$CIGAR == allele, c('Amplicon', 'CIGAR')], collapse = ','),][['alleles']] = 
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
      
      asv_table = asv_table[order(as.integer(gsub('ASV','',asv_table$hapid))),]
      
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
  asv_seqs_masked = 'ANY',
  vcf_like = 'ANY',
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
                         asv_seqs_masked = NULL,
                         vcf_like = NULL,
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
  obj@asv_seqs_masked = asv_seqs_masked
  obj@vcf_like = vcf_like
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
    
    if(sum(grepl('(/|-|:)', markers[["amplicon"]])) > 0){
      
      markers[["amplicon"]] = gsub('(/|-|:)', '_', markers[["amplicon"]])
      print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the column amplicon of the markers table.')
      
    }
    
    ampseq_loci_vector = markers[["amplicon"]]
    
  }else if(!is.null(markers_pattern)){
    all_alleles = gsub(markers_pattern, '', colnames(cigar_table))
    all_alleles = paste0(paste(c(all_alleles[all_alleles != '.'], '\\.'), collapse = '|'), '$')
    ampseq_loci_vector = unique(gsub(all_alleles, '',colnames(cigar_table)))
    
  }else if(is.null(markers) & is.null(markers_pattern)){
    ampseq_loci_vector = unique(sapply(strsplit(colnames(cigar_table), ","), function(x) x[1]))
    markers = data.frame(amplicon = ampseq_loci_vector, length=NA)
  }

  # Check - If cigar_table and ampseq_loci_vector are not equal, recorrect ampseq_loci_vector to prevent errors propagated by different markersTable information
  # if(!setequal(unique(colnames(cigar_table)), unique(ampseq_loci_vector))){
  #   ampseq_loci_vector = unique(sapply(strsplit(colnames(cigar_table), ","), function(x) x[1]))
  # }
  
  ampseq_loci_abd_table = matrix(NA, nrow = nrow(cigar_table), ncol = length(ampseq_loci_vector), dimnames = list(rownames(cigar_table), ampseq_loci_vector))
  
  for(sample in rownames(ampseq_loci_abd_table)){
    for(locus in colnames(ampseq_loci_abd_table)){
      alleles = cigar_table[sample, grepl(paste0("^",locus,',|\\.)'), colnames(cigar_table))]
      
      if(length(alleles) == 1){
        names(alleles) <- colnames(cigar_table)[grepl(paste0("^",locus,',|\\.)'), colnames(cigar_table))]
      }
      
      if(length(alleles[which(alleles >= min_abd)]) == 1){
        
        if(!is.null(markers_pattern)){
          ampseq_loci_abd_table[sample, locus] = paste(gsub(locus,'',names(alleles[which(alleles >= min_abd)]), ","), alleles[which(alleles >= min_abd)], sep = ":")
        }else{
          ampseq_loci_abd_table[sample, locus] = paste(gsub(paste0(locus, '.'), '', names(alleles[which(alleles >= min_abd)]), ","), alleles[which(alleles >= min_abd)], sep = ":")
        }
        
      }else if(length(alleles[which(alleles >= min_abd)]) > 1){
        allele_names = names(alleles)
        alleles = as.integer(alleles)
        names(alleles) = allele_names
        alleles = sort(alleles, decreasing = T)
        
        if(!is.null(markers_pattern)){
          ampseq_loci_abd_table[sample, locus] = gsub(locus, "", paste(paste(names(alleles[alleles/max(alleles) >= min_ratio]), alleles[alleles/max(alleles) >= min_ratio], sep = ":"), collapse = "_"))  
        }else{
          ampseq_loci_abd_table[sample, locus] = gsub(paste0(locus, "."), "", paste(paste(names(alleles[alleles/max(alleles) >= min_ratio]), alleles[alleles/max(alleles) >= min_ratio], sep = ":"), collapse = "_"))  
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
                                  asv_seqs_masked = asv_seqs,
                                  metadata = metadata,
                                  controls = list(gt = controls_ampseq_loci_abd_table,
                                                  metadata = controls_metadata),
                                  markers = markers,
                                  loci_performance = NULL,
                                  pop_summary = NULL)
    
    cigars_in_gt = get_cigar_alleles(ampseq_object, 
                                           from = 'gt',
                                           as = 'vector')
    
    cigars_in_asvtab = get_cigar_alleles(ampseq_object, 
                                               from = 'asv_table',
                                               as = 'vector',
                                               col = 'CIGAR_masked')
    
    if(sum(!(cigars_in_asvtab %in% cigars_in_gt)) > 0){
      
      removed_alleles = cigars_in_asvtab[!(cigars_in_asvtab %in% cigars_in_gt)]
      
      for(allele in removed_alleles){
        
        Amplicon = gsub(';.+$', '', allele)
        CIGAR_masked = gsub('^.+;', '', allele)
        
        # ampseq_object@asv_table[ampseq_object@asv_table[['Amplicon']] == Amplicon &
        #                           ampseq_object@asv_table[['CIGAR_masked']] == CIGAR_masked &
        #                           !is.na(ampseq_object@asv_table[['CIGAR_masked']]),
        # ][['CIGAR_masked']] = NA
        
        ampseq_object@asv_table = ampseq_object@asv_table[!(ampseq_object@asv_table[['Amplicon']] == Amplicon &
                                                              ampseq_object@asv_table[['CIGAR_masked']] == CIGAR_masked),]
        
        ampseq_object@asv_seqs = ampseq_object@asv_seqs[names(ampseq_object@asv_seqs) %in%
                                                          ampseq_object@asv_table$hapid]
        ampseq_object@asv_seqs_masked = ampseq_object@asv_seqs_masked[names(ampseq_object@asv_seqs_masked) %in%
                                                          ampseq_object@asv_table$hapid]
        
        ampseq_object@asv_table$hapid = paste0('ASV', 1:nrow(ampseq_object@asv_table))
        names(ampseq_object@asv_seqs) = ampseq_object@asv_table$hapid
        names(ampseq_object@asv_seqs_masked) = ampseq_object@asv_table$hapid
        
      }
      
    }
    
    return(ampseq_object)
    
  }else{
    
    metadata = cigar_object@metadata
    
    ampseq_object = create_ampseq(gt = ampseq_loci_abd_table,
                                  asv_table = asv_table,
                                  asv_seqs = asv_seqs,
                                  asv_seqs_masked = asv_seqs,
                                  metadata = metadata,
                                  markers = markers,
                                  loci_performance = NULL,
                                  pop_summary = NULL)
    
    cigars_in_gt = get_cigar_alleles(ampseq_object, 
                                     from = 'gt',
                                     as = 'vector')
    
    cigars_in_asvtab = get_cigar_alleles(ampseq_object, 
                                         from = 'asv_table',
                                         as = 'vector',
                                         col = 'CIGAR_masked')
    
    if(sum(!(cigars_in_asvtab %in% cigars_in_gt)) > 0){
      
      removed_alleles = cigars_in_asvtab[!(cigars_in_asvtab %in% cigars_in_gt)]
      
      for(allele in removed_alleles){
        
        Amplicon = gsub(';.+$', '', allele)
        CIGAR_masked = gsub('^.+;', '', allele)
        
        # ampseq_object@asv_table[ampseq_object@asv_table[['Amplicon']] == Amplicon &
        #                           ampseq_object@asv_table[['CIGAR_masked']] == CIGAR_masked &
        #                           !is.na(ampseq_object@asv_table[['CIGAR_masked']]),
        # ][['CIGAR_masked']] = NA
        
        ampseq_object@asv_table = ampseq_object@asv_table[!(ampseq_object@asv_table[['Amplicon']] == Amplicon &
                                  ampseq_object@asv_table[['CIGAR_masked']] == CIGAR_masked),]
        
        ampseq_object@asv_seqs = ampseq_object@asv_seqs[names(ampseq_object@asv_seqs) %in%
                                                          ampseq_object@asv_table$hapid]
        
        ampseq_object@asv_seqs_masked = ampseq_object@asv_seqs_masked[names(ampseq_object@asv_seqs_masked) %in%
                                                          ampseq_object@asv_table$hapid]
        
      }
      
    }
    
    return(ampseq_object)
    
  }
  
}

## ampseq2cigar----

ampseq2cigar = function(ampseq_object){
  
  cigar_object = create_cigar()
  
  gt = ampseq_object@gt
  
  long_cigar = NULL
  
  for(mhap in 1:ncol(gt)){
    
    temp_cigar_long = unlist(strsplit(gt[,mhap], '_'))
    
    nclones = sapply(gt[,mhap], function(Sample){
      sum(nchar(unlist(str_extract_all(Sample, '_'))))
    }) + 1
    
    nclones[is.na(nclones)] = 1
    
    temp_cigar_long = data.frame(Sample_id = rep(names(nclones), nclones),
                                 Amplicon = colnames(gt)[mhap],
                                 Cigar_allele = gsub(':\\d+$', '', temp_cigar_long),
                                 Read_depth = as.integer(gsub('^.+:', '', temp_cigar_long))
                                 )
    
    temp_cigar_long %<>% filter(!is.na(Cigar_allele))
    
    long_cigar = rbind(long_cigar,
                       temp_cigar_long)
     
  }
  
  long_cigar %<>% mutate(Cigar_allele = paste(Amplicon, Cigar_allele, sep = ','))
  
  long_cigar$Amplicon = NULL
  
  wide_cigar = pivot_wider(long_cigar,
                           names_from = 'Cigar_allele',
                           values_from = 'Read_depth')
  
  cigar_object@cigar_table = wide_cigar
  cigar_object@metadata = ampseq_object@metadata
  cigar_object@asv_table = ampseq_object@asv_table
  cigar_object@asv_seqs = ampseq_object@asv_seqs
  
  return(cigar_object)
  
}

## join_ampseq----

# Use only if the objects to merge haven't been masked

join_ampseq = function(ampseq_obj_list = NULL, remove_replicates = TRUE){
  
  gt = NULL
  asv_table = NULL
  asv_seqs = NULL
  metadata = NULL
  markers = NULL
  
  for(obj in 1:length(ampseq_obj_list)){
    
    #if(obj == 2){stop()}
    
    obj = ampseq_obj_list[[obj]]
    
    
    if(is.null(gt)){ # for the first object
      
      gt = obj@gt
      asv_table = obj@asv_table
      asv_seqs = obj@asv_seqs
      metadata = obj@metadata
      markers = obj@markers
      
      asv_table$hapid = paste0('ASV', 1:nrow(asv_table))
      names(asv_seqs) = asv_table$hapid
      
    }else{# for subsequent objects
      
      temp_gt1 = obj@gt
      temp_asv_table1 = obj@asv_table
      temp_asv_seqs1 = obj@asv_seqs
      temp_metadata1 = obj@metadata
      temp_markers1 = obj@markers
      
      temp_asv_table1$hapid = paste0('ASV', 1:nrow(temp_asv_table1))
      names(temp_asv_seqs1) = temp_asv_table1$hapid
      
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

        unconsitent_cigar_strings = data.frame(temp_asv_table2[temp_asv_table2$CIGAR != temp_asv_table3$CIGAR,], # Information from the new dataset
                                               Amplicon2 = temp_asv_table3[temp_asv_table3$CIGAR != temp_asv_table2$CIGAR,][['Amplicon']], # Name of the amplicon in the previous  dataset
                                               CIGAR2 = temp_asv_table3[temp_asv_table3$CIGAR != temp_asv_table2$CIGAR,][['CIGAR']], # Original cigar string in the previous dataset
                                               CIGAR_masked2 = temp_asv_table3[temp_asv_table3$CIGAR != temp_asv_table2$CIGAR,][['CIGAR_masked']]) # Masked cigar string in the previous dataset
        
        
        
        if(nrow(unconsitent_cigar_strings) > 0){
          
          for(pos in 1:nrow(unconsitent_cigar_strings)){
            
            cigar_string_masked_replacment = unconsitent_cigar_strings[pos, ][['CIGAR_masked2']]
            cigar_string_masked_pattern = unconsitent_cigar_strings[pos, ][['CIGAR_masked']]
            
            #if(cigar_string_masked_replacment == "14G18A19G21C26T"){stop("Wrong_pattern")}
            
            mhap = unconsitent_cigar_strings[pos, ][['Amplicon']]
            
            
            if(mhap != unconsitent_cigar_strings[pos, ][['Amplicon2']]){stop()}
            #if(mhap == 'pvpi3k_2' & cigar_string_masked_pattern == '26C30T84T'){stop()}
            
            for(samp in 1:nrow(temp_gt1)){
              
              cigar_string_masked_pattern_in_samp = 
                str_extract(temp_gt1[samp, mhap], 
                            paste0('(^|_)',
                                   ifelse(cigar_string_masked_pattern == '.', '\\.', cigar_string_masked_pattern), ':'))
              
              if(!is.na(cigar_string_masked_pattern_in_samp)){
                
                if(cigar_string_masked_pattern_in_samp == '.:'){
                  cigar_string_masked_pattern_in_samp = '\\.:'
                }
                
                if(!grepl('^_',cigar_string_masked_pattern_in_samp) & cigar_string_masked_pattern_in_samp != '\\.:'){
                  
                  cigar_string_masked_pattern_in_samp = paste0('^', cigar_string_masked_pattern_in_samp)
                  
                }
                
                cigar_string_masked_replacment_in_samp = gsub(ifelse(cigar_string_masked_pattern == '.', '\\.', cigar_string_masked_pattern),
                                                              cigar_string_masked_replacment,
                                                              ifelse(cigar_string_masked_pattern_in_samp == "\\.:", '.:', cigar_string_masked_pattern_in_samp))
                
                cigar_sample_replacement = gsub('\\^', '',gsub(cigar_string_masked_pattern_in_samp,
                                            cigar_string_masked_replacment_in_samp, 
                                            temp_gt1[samp, mhap]))
                
                new_alleles = unlist(strsplit(gsub(':\\d+', '', cigar_sample_replacement), '_'))
                
                if(sum(duplicated(new_alleles))){
                  
                  new_allele_reads = as.integer(gsub(':', '', unlist(str_extract_all(cigar_sample_replacement, ':\\d+'))))
                  
                  duplicated_new_alleles = new_alleles[duplicated(new_alleles)]
                  
                  unique_new_alleles = new_alleles[!duplicated(new_alleles)]
                  unique_new_allele_reads = new_allele_reads[!duplicated(new_alleles)]
                  
                  for(duplicated_new_allele in duplicated_new_alleles){
                    
                    unique_new_allele_reads[which(unique_new_alleles %in% duplicated_new_allele)] = sum(new_allele_reads[which(new_alleles %in% duplicated_new_allele)])
                  }
                  
                  cigar_sample_replacement = paste(paste(unique_new_alleles, unique_new_allele_reads, sep = ':'), collapse = '_')
                  
                }
                
                temp_gt1[samp, mhap] = cigar_sample_replacement
                
                rm(cigar_string_masked_replacment_in_samp)
                
              }
              rm(list = c('samp', 'cigar_string_masked_pattern_in_samp'))
            }
            
            rm(list = c('pos', 
                        'mhap', 
                        'cigar_string_masked_replacment', 
                        'cigar_string_masked_pattern'
            ))
            
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
      asv_table = asv_table[order(as.integer(gsub('ASV','',asv_table$hapid))),]
      
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
        
        shared_loci = markers$amplicon[!(markers$amplicon %in% unshared_loci)]
        
        unshared_loci_with_respc_tablen = markers$amplicon[(markers$amplicon %in% unshared_loci)]
        
        unshared_loci_in_tablen = temp_markers1$amplicon[(temp_markers1$amplicon %in% unshared_loci)]
        
        shared_loci_df = markers[markers$amplicon%in% shared_loci,]
        
        unshared_loci_with_tablen_df = markers[markers$amplicon %in% unshared_loci_with_respc_tablen, ]
        
        unshared_loci_in_tablen_df = temp_markers1[temp_markers1$amplicon %in% unshared_loci_in_tablen, names(markers)]
        
        markers = rbind(shared_loci_df, unshared_loci_with_tablen_df, unshared_loci_in_tablen_df)
        markers %<>% arrange(chromosome, start)
        
        rm(list = c('shared_loci_df', 
                    'unshared_loci_with_tablen_df', 
                    'unshared_loci_in_tablen_df',
                    'shared_loci',
                    'unshared_loci_with_respc_tablen',
                    'unshared_loci_in_tablen'
                    ))
        
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
  
  #sum(grepl('26C30T84T', gt[,'pvpi3k_2']))
  
  #asv_table[asv_table$Amplicon == 'pvpi3k_2',]
  
  # Remove duplicates
  if(remove_replicates){
    
    print("Duplicated samples are being removed...")
    
    duplicated_samples = metadata[duplicated(metadata$Sample_id),][['Sample_id']]
    
    gt_replicates = gt[duplicated(metadata$Sample_id),]
    gt = gt[!duplicated(metadata$Sample_id),]
    
    metadata_replicates = metadata[duplicated(metadata$Sample_id),]
    metadata = metadata[!duplicated(metadata$Sample_id),]
    
    rownames(metadata) = metadata[['Sample_id']]
    rownames(gt) = metadata[['Sample_id']]
    
    rownames(metadata_replicates) = metadata_replicates[['Sample_id']]
    rownames(gt_replicates) = metadata_replicates[['Sample_id']]
    
    for(sample in duplicated_samples){
      
      temp_data = rbind(gt[sample,],
                        gt_replicates[sample,])
      
      temp_metadata = rbind(metadata[metadata$Sample_id == sample, ],
                            metadata_replicates[metadata_replicates$Sample_id == sample, ])
      
      temp_test = rowSums(!is.na(temp_data), na.rm = T)
      
      gt[sample,] = temp_data[which.max(temp_test),]
      metadata[metadata$Sample_id == sample, ] = temp_metadata[which.max(temp_test),]
      
    }
    
  }
## ADD VERIFICATION STEP OF "gt", "asv_table" AND "asv_seqs slots"
  
  ampseq_object = create_ampseq(gt = gt,
                                asv_table = asv_table,
                                asv_seqs = asv_seqs,
                                asv_seqs_masked = asv_seqs,
                                metadata = metadata,
                                markers = markers)
  
  return(ampseq_object)
  
}

## write_ampseq----

write_ampseq = function(ampseq_object, format = c('excel', 'csv', 'json', 'tsv'), name = 'wb.xlsx', sep = '\t'){
  
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
                       'asv_seqs_masked',
                       'vcf_like',
                       'discarded_loci',
                       'discarded_samples',
                       'controls')){
      
      if(temp_slot == 'gt'){
        print("Printing gt slot...")

        if(!is.null(slot(ampseq_object, temp_slot))){
          temp_sheet = data.frame(Sample_id = rownames(slot(ampseq_object, temp_slot)),
                                as.data.frame(slot(ampseq_object, temp_slot)))
        } else {
          temp_sheet = NULL
        } 
      }else if(temp_slot %in% c('asv_seqs', 'asv_seqs_masked')){
        print("Printing asv_seqs slot...")
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_sheet = data.frame(asv_id = names(slot(ampseq_object, temp_slot)),
                                  asv_seq = as.character(slot(ampseq_object, temp_slot)))
          
        }else{
          temp_sheet = NULL
        }
        
      }else if(temp_slot == 'markers'){
        print("Printing markers slot...")
        if(!is.null(slot(ampseq_object, temp_slot))) {
          temp_sheet = as.data.frame(slot(ampseq_object, temp_slot))
          
          if(sum(is.infinite(temp_sheet[['distance']])) > 0){
            temp_sheet[is.infinite(temp_sheet[['distance']]),][['distance']] = NA
          }
        } else {
          temp_sheet = NULL
        }
        
      }else if(temp_slot == 'discarded_loci'){
        print("Printing discarded_loci slots...")
        
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
          
          # Before setting values to NA, check if any infinite values exist
          if(sum(is.infinite(temp_discarded_loci_markers[['distance']])) > 0){
            temp_discarded_loci_markers[is.infinite(temp_discarded_loci_markers[['distance']]),][['distance']] = NA
          }
          
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
        print("Printing discarded_samples slots...")
        temp_sheet = NULL
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_discarded_samples = slot(ampseq_object, temp_slot)
          
          # write gt
          
           if (!is.null(temp_discarded_samples[['gt']])) {
            temp_discarded_samples_gt <- data.frame(Sample_id = rownames(temp_discarded_samples[['gt']]), 
                                                    as.data.frame(temp_discarded_samples[['gt']]))
          } else {
            temp_discarded_samples_gt <- data.frame(Sample_id = NULL)
          }
          
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
        print("Printing controls slots...")
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
                       'asv_seqs_masked',
                       'vcf_like',
                       'discarded_loci',
                       'discarded_samples',
                       'controls')){
      
      if(temp_slot == 'gt'){
        
        temp_sheet = data.frame(Sample_id = rownames(slot(ampseq_object, temp_slot)),
                                as.data.frame(slot(ampseq_object, temp_slot)))
        
      }else if(temp_slot %in% c('asv_seqs', 'asv_seqs_masked')){
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_sheet = data.frame(asv_id = names(slot(ampseq_object, temp_slot)),
                                  asv_seq = as.character(slot(ampseq_object, temp_slot)))
          
        }else{
          temp_sheet = NULL
        }
        
      }else if(temp_slot == 'markers'){
        
        temp_sheet = as.data.frame(slot(ampseq_object, temp_slot))
        
        if(sum(is.infinite(temp_sheet[['distance']])) > 0){
          temp_sheet[is.infinite(temp_sheet[['distance']]),][['distance']] = NA
        }
        
        
      }else if(temp_slot == 'discarded_loci'){
        
        temp_sheet = NULL
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_discarded_loci = slot(ampseq_object, temp_slot)
          
          # write gt
          
          if (!is.null(temp_discarded_loci[['gt']])) {
            temp_discarded_loci_gt = data.frame(Sample_id = rownames(temp_discarded_loci[['gt']]),
                                              as.data.frame(temp_discarded_loci[['gt']]))
          } else {
            temp_discarded_loci_gt = data.frame(Sample_id = NULL)
          }
          
          write.csv(temp_discarded_loci_gt,
                    file.path(name, 'discarded_loci_gt.csv'), 
                    quote = F, row.names = F)
          
          # write markers
          
          temp_discarded_loci_markers = as.data.frame(temp_discarded_loci[['markers']])
          
          if(sum(is.infinite(temp_discarded_loci_markers[['distance']])) > 0){
            temp_discarded_loci_markers[is.infinite(temp_discarded_loci_markers[['distance']]),][['distance']] = NA  
          }
          
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
    
    
  }else if(format == 'tsv'){
    
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
                       'asv_seqs_masked',
                       'vcf_like',
                       'discarded_loci',
                       'discarded_samples',
                       'controls')){
      
      if(temp_slot == 'gt'){
        
        temp_sheet = data.frame(Sample_id = rownames(slot(ampseq_object, temp_slot)),
                                as.data.frame(slot(ampseq_object, temp_slot)))
        
      }else if(temp_slot %in% c('asv_seqs', 'asv_seqs_masked')){
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_sheet = data.frame(asv_id = names(slot(ampseq_object, temp_slot)),
                                  asv_seq = as.character(slot(ampseq_object, temp_slot)))
          
        }else{
          temp_sheet = NULL
        }
        
      }else if(temp_slot == 'markers'){
        
        temp_sheet = as.data.frame(slot(ampseq_object, temp_slot))
        
        if(sum(is.infinite(temp_sheet[['distance']])) > 0){
          temp_sheet[is.infinite(temp_sheet[['distance']]),][['distance']] = NA
        }
        
        
      }else if(temp_slot == 'discarded_loci'){
        
        temp_sheet = NULL
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_discarded_loci = slot(ampseq_object, temp_slot)
          
          # write gt
          
          if(!is.null(nrow(temp_discarded_loci[['gt']]))){
            temp_discarded_loci_gt = data.frame(Sample_id = rownames(temp_discarded_loci[['gt']]),
                                                as.data.frame(temp_discarded_loci[['gt']]))
          }else{
            
            temp_discarded_loci_gt = data.frame(Sample_id = names(temp_discarded_loci[['gt']]),
                                                as.data.frame(temp_discarded_loci[['gt']]))
            
          }
          
          
          
          write.table(temp_discarded_loci_gt,
                    file.path(name, 'discarded_loci_gt.tsv'), 
                    quote = F, row.names = F, sep = sep)
          
          # write markers
          
          temp_discarded_loci_markers = as.data.frame(temp_discarded_loci[['markers']])
          
          if(sum(is.infinite(temp_discarded_loci_markers[['distance']])) > 0){
            temp_discarded_loci_markers[is.infinite(temp_discarded_loci_markers[['distance']]),][['distance']] = NA  
          }
          
          write.table(temp_discarded_loci_markers,
                    file.path(name, 'discarded_loci_markers.tsv'), 
                    quote = F, row.names = F, sep = sep)
          
          # write loci_performance
          
          temp_discarded_loci_loci_performance = as.data.frame(temp_discarded_loci[['loci_performance']])
          
          write.table(temp_discarded_loci_loci_performance,
                    file.path(name, 'discarded_loci_loci_performance.tsv'), 
                    quote = F, row.names = F, sep = sep)
          
        }
        
        
      }else if(temp_slot == 'discarded_samples'){
        
        temp_sheet = NULL
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_discarded_samples = slot(ampseq_object, temp_slot)
          
          # write gt
          
          temp_discarded_samples_gt = data.frame(Sample_id = rownames(temp_discarded_samples[['gt']]),
                                              as.data.frame(temp_discarded_samples[['gt']]))
          
          write.table(temp_discarded_samples_gt,
                    file.path(name, 'discarded_samples_gt.tsv'), 
                    quote = F, row.names = F, sep = sep)
          
          # write metadata
          
          temp_discarded_samples_metadata = as.data.frame(temp_discarded_samples[['metadata']])
          
          write.table(temp_discarded_samples_metadata,
                    file.path(name, 'discarded_samples_metadata.tsv'),
                    quote = F, row.names = F, sep = sep)
        }
        
      }else if(temp_slot == 'controls'){
        
        temp_sheet = NULL
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          
          temp_controls = slot(ampseq_object, temp_slot)
          
          # write gt
          
          temp_controls_gt = data.frame(Sample_id = rownames(temp_controls[['gt']]),
                                        as.data.frame(temp_controls[['gt']]))
          
          write.table(temp_controls_gt,
                    file.path(name, 'controls_gt.csv'), 
                    quote = F, row.names = F, sep = sep)
          
          # write metadata
          
          temp_controls_metadata = as.data.frame(temp_controls[['metadata']])
          
          write.table(temp_controls_metadata,
                    file.path(name, 'controls_metadata.csv'),
                    quote = F, row.names = F, sep = sep)
        }
        
      }else{
        
        if(!is.null(slot(ampseq_object, temp_slot))){
          temp_sheet = as.data.frame(slot(ampseq_object, temp_slot))
        }else{
          temp_sheet = NULL
        }
        
        
      }
      
      if(!is.null(temp_sheet)){
        
        write.table(temp_sheet, paste0(file.path(name, temp_slot), '.tsv'), quote = F, row.names = F, sep = sep)
        
      }
      
      
    }
    
    
  }else if(format == 'json'){
    # In development
    
    # library(jsonlite)
    # 
    # ampseq_object2 = NULL
    # 
    # for(temp_slot in slotNames(ampseq_object)){
    #   
    #   ampseq_object2[[temp_slot]] = slot(ampseq_object, temp_slot)
    #   
    # }
    # 
    # ampseq_json =  toJSON(ampseq_object2$asv_table)
  }
  
}

## write_cigar----

write_cigar = function(cigar_object, format = c('excel', 'csv', 'json', 'tsv'), name = 'wb.xlsx', sep = '\t'){
  
  if(format == 'excel'){
    
    if(file.exists(name)){
      system(paste0('rm ', name))
    }
    
    excel_wb = loadWorkbook(name, create = T)
    
    for(temp_slot in c('cigar_table', 
                       'metadata', 
                       'asv_table',
                       'asv_seqs')){
      
      if(temp_slot == 'cigar_table'){
        
        temp_sheet = data.frame(Sample_id = rownames(slot(cigar_object, temp_slot)),
                                as.data.frame(slot(cigar_object, temp_slot)))
        
      }else if(temp_slot == 'asv_seqs'){
        
        if(!is.null(slot(cigar_object, temp_slot))){
          
          temp_sheet = data.frame(asv_id = names(slot(cigar_object, temp_slot)),
                                  asv_seq = as.character(slot(cigar_object, temp_slot)))
          
        }else{
          temp_sheet = NULL
        }
        
      }else{
        
        if(!is.null(slot(cigar_object, temp_slot))){
          temp_sheet = as.data.frame(slot(cigar_object, temp_slot))
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
    
    for(temp_slot in c('cigar_table', 
                       'metadata',
                       'asv_table', 
                       'asv_seqs')){
      
      if(temp_slot == 'cigar_table'){
        
        temp_sheet = data.frame(Sample_id = rownames(slot(cigar_object, temp_slot)),
                                as.data.frame(slot(cigar_object, temp_slot)))
        
      }else if(temp_slot == 'asv_seqs'){
        
        if(!is.null(slot(cigar_object, temp_slot))){
          
          temp_sheet = data.frame(asv_id = names(slot(cigar_object, temp_slot)),
                                  asv_seq = as.character(slot(cigar_object, temp_slot)))
          
        }else{
          temp_sheet = NULL
        }
        
      }else{
        
        if(!is.null(slot(cigar_object, temp_slot))){
          temp_sheet = as.data.frame(slot(cigar_object, temp_slot))
        }else{
          temp_sheet = NULL
        }
        
        
      }
      
      if(!is.null(temp_sheet)){
        
        write.csv(temp_sheet, paste0(file.path(name, temp_slot), '.csv'), quote = F, row.names = F)
        
      }
      
      
    }
    
    
  }else if(format == 'tsv'){
    
    if(file.exists(name)){
      system(paste0('rm -r ', name))
    }
    
    system(paste0('mkdir ', name))
    
    for(temp_slot in c('cigar_table', 
                       'metadata',
                       'asv_table', 
                       'asv_seqs')){
      
      if(temp_slot == 'cigar_table'){
        
        temp_sheet = data.frame(Sample_id = rownames(slot(cigar_object, temp_slot)),
                                as.data.frame(slot(cigar_object, temp_slot)))
        
      }else if(temp_slot == 'asv_seqs'){
        
        if(!is.null(slot(cigar_object, temp_slot))){
          
          temp_sheet = data.frame(asv_id = names(slot(cigar_object, temp_slot)),
                                  asv_seq = as.character(slot(cigar_object, temp_slot)))
          
        }else{
          temp_sheet = NULL
        }
        
      }else{
        
        if(!is.null(slot(cigar_object, temp_slot))){
          temp_sheet = as.data.frame(slot(cigar_object, temp_slot))
        }else{
          temp_sheet = NULL
        }
        
        
      }
      
      if(!is.null(temp_sheet)){
        
        write.table(temp_sheet, paste0(file.path(name, temp_slot), '.tsv'), quote = F, row.names = F, sep = sep)
        
      }
      
      
    }
    
    
  }else if(format == 'json'){
    # In development
    
    # library(jsonlite)
    # 
    # cigar_object2 = NULL
    # 
    # for(temp_slot in slotNames(cigar_object)){
    #   
    #   ampseq_object2[[temp_slot]] = slot(cigar_object, temp_slot)
    #   
    # }
    # 
    # ampseq_json =  toJSON(cigar_object2$asv_table)
  }
  
}

## read_ampseq----

read_ampseq = function(file = NULL, format = 'excel', sep = '\t'){
  
  ampseq_object = create_ampseq()
  
  if(format == 'excel'){
    
    temp_wb = loadWorkbook(file)
    
    for(sheet in getSheets(temp_wb)){
      print(paste0("Loading sheet ", sheet))
      if(sheet == 'gt'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = as.matrix(temp_sheet[,-1])
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(ampseq_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet %in% c('asv_seqs', 'asv_seqs_masked')){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_names = temp_sheet[[1]]
        temp_sheet = DNAStringSet(temp_sheet[[2]])
        names(temp_sheet) = temp_sheet_names
        
        slot(ampseq_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet %in% c('metadata', 'loci_performance', 'vcf_like')){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(ampseq_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet == 'markers'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        if (any(is.na(temp_sheet[['distance']]))){
          temp_sheet[['distance']][is.na(temp_sheet[['distance']])] = Inf 
        }
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
          if (any(is.na(temp_sheet[['distance']]))){
            temp_sheet[['distance']][is.na(temp_sheet[['distance']])] = Inf 
          }
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
        
      }else if(sheet %in% c('asv_seqs.csv', 'asv_seqs_masked.csv')){
        
        temp_sheet = read.csv(file.path(file, sheet))
        temp_sheet_names = temp_sheet[[1]]
        temp_sheet = DNAStringSet(temp_sheet[[2]])
        names(temp_sheet) = temp_sheet_names
        
        slot(ampseq_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet %in% c('metadata.csv', 'loci_performance.csv', 'vcf_like.csv')){
        
        temp_sheet = read.csv(file.path(file, sheet))
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(ampseq_object, gsub('.csv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet == 'markers.csv'){
        
        temp_sheet = read.csv(file.path(file, sheet))
        if(sum(is.na(temp_sheet[['distance']]))>0){
          temp_sheet[is.na(temp_sheet[['distance']]),][['distance']] = Inf  
        }
        
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
          if(sum(is.na(temp_sheet[['distance']]))>0){
            temp_sheet[is.na(temp_sheet[['distance']]),][['distance']] = Inf
          }
          
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
    
  }else if(format == 'tsv'){
    
    for(sheet in list.files(file)){
      if(sheet == 'gt.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), sep = sep, header = T)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = as.matrix(temp_sheet[,-1])
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(ampseq_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet %in% c('asv_seqs.tsv', 'asv_seqs_masked.tsv')){
        
        temp_sheet = read.table(file.path(file, sheet), sep = sep, header = T)
        temp_sheet_names = temp_sheet[[1]]
        temp_sheet = DNAStringSet(temp_sheet[[2]])
        names(temp_sheet) = temp_sheet_names
        
        slot(ampseq_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet %in% c('metadata.tsv', 'loci_performance.tsv', 'vcf_like.tsv')){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(ampseq_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet == 'markers.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        
        if(sum(is.na(temp_sheet[['distance']])) > 0){
          temp_sheet[is.na(temp_sheet[['distance']]),][['distance']] = Inf
        }
        slot(ampseq_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(grepl('discarded_loci', sheet)){
        
        if(sheet == 'discarded_loci_gt.tsv'){
          temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
          temp_sheet_rownames = temp_sheet[,1]
          temp_sheet = as.matrix(temp_sheet[,-1])
          rownames(temp_sheet) = temp_sheet_rownames
          
          ampseq_object@discarded_loci[['gt']] = temp_sheet
          
        }else if(sheet == 'discarded_loci_markers.tsv'){
          
          temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
          
          if(sum(is.na(temp_sheet[['distance']])) > 0){
            temp_sheet[is.na(temp_sheet[['distance']]),][['distance']] = Inf
          }
          
          ampseq_object@discarded_loci[['markers']] = temp_sheet
          
        }else if(sheet == 'discarded_loci_loci_performance.tsv'){
          
          temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
          temp_sheet_rownames = temp_sheet[,1]
          rownames(temp_sheet) = temp_sheet_rownames
          ampseq_object@discarded_loci[['loci_performance']] = temp_sheet
          
        }
        
      }else if(grepl('discarded_samples', sheet)){
        
        if(sheet == 'discarded_samples_gt.tsv'){
          temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
          temp_sheet_rownames = temp_sheet[,1]
          temp_sheet = as.matrix(temp_sheet[,-1])
          rownames(temp_sheet) = temp_sheet_rownames
          
          ampseq_object@discarded_samples[['gt']] = temp_sheet
          
        }else if(sheet == 'discarded_samples_metadata.tsv'){
          
          temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
          temp_sheet_rownames = temp_sheet[,1]
          rownames(temp_sheet) = temp_sheet_rownames
          ampseq_object@discarded_samples[['metadata']] = temp_sheet
          
        }
        
      }else if(grepl('controls', sheet)){
        
        if(sheet == 'controls_gt.tsv'){
          temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
          temp_sheet_rownames = temp_sheet[,1]
          temp_sheet = as.matrix(temp_sheet[,-1])
          rownames(temp_sheet) = temp_sheet_rownames
          
          ampseq_object@controls[['gt']] = temp_sheet
          
        }else if(sheet == 'controls_metadata.tsv'){
          
          temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
          temp_sheet_rownames = temp_sheet[,1]
          rownames(temp_sheet) = temp_sheet_rownames
          ampseq_object@controls[['metadata']] = temp_sheet
          
        }
        
      }else{
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        slot(ampseq_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }
    }
    
  }else if(format == 'json'){
    # In development
  }
  
  return(ampseq_object)
  
}



## loci S4class and create_loci----

## loci S4 class

setClass('loci', slots = c(
  loci_table = "ANY",
  metadata = "ANY",
  freq_table = 'ANY',
  ampseq_alleles_per_locus = 'ANY',
  markers = 'ANY'
  
))

## loci constructor
create_loci = function(loci_table = NULL,
                        metadata = NULL,
                        freq_table = NULL,
                        ampseq_alleles_per_locus = NULL,
                        markers = NULL
                        ){
  
  obj = new('loci')
  obj@loci_table = loci_table
  obj@metadata = metadata
  obj@freq_table = freq_table
  obj@ampseq_alleles_per_locus = ampseq_alleles_per_locus
  obj@markers = markers
  
  return(obj)
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
  
  loci_object = create_loci(
    loci_table = loci_table,
    metadata = ampseq_object@metadata,
    freq_table = freq_table,
    ampseq_alleles_per_locus = ampseq_allele_freq,
    markers = ampseq_object@markers)
  
  return(loci_object)
  
}


## ampseq2vcf----

ampseq2vcf = function(ampseq_object, monoclonals = NULL, polyclonals = NULL, ref_fasta){
  
  markers = ampseq_object@markers
  asv_table = ampseq_object@asv_table
  gt_ampseq = ampseq_object@gt
  
  reference_genome = readDNAStringSet(ref_fasta)
  
  names(reference_genome) = gsub(' \\| .+', '', names(reference_genome))
  
  loci_table = NULL
  gt = NULL
  
  if(is.null(monoclonals) & is.null(polyclonals)){
    monoclonals = ampseq_object@metadata$Sample_id
  }
  
  # amplicon = markers$amplicon[1]
  # 
  # amplicon = "PVP01_0109700_2"

  # 
  # which(markers$amplicon == amplicon)
  # 
  #amplicon = markers$amplicon[1]
  for(amplicon in markers$amplicon){
    
    print(amplicon)
    
    #View(asv_table[asv_table$Amplicon == amplicon,])
    
    amplicon_cigars = asv_table[asv_table$Amplicon == amplicon,][['CIGAR_masked']]
    
    amplicon_cigars = unique(amplicon_cigars[!is.na(amplicon_cigars) &
                                               amplicon_cigars != '.'])
    
    if(length(amplicon_cigars) != 0){
      
      cigar_position_variants = NULL
      
      for(cigar in amplicon_cigars){
        cigar_position_variants = 
          rbind(cigar_position_variants,
                data.frame(cigar = cigar,
                           CHROM = markers[markers$amplicon == amplicon,][['chromosome']],
                           amplicon_POS = unlist(str_extract_all(cigar, '\\d+')),
                           POS = as.integer(unlist(str_extract_all(cigar, '\\d+'))) + 
                             markers[markers$amplicon == amplicon,][['start']] - 1,
                           ALT = unlist(str_extract_all(cigar, '[DI]?=?[ATGC]+')),
                           ALT_label = NA
                ))
      }
      
      amplicon_loci_table = NULL
      
      
      for(position in as.character(sort(as.integer(unique(cigar_position_variants$amplicon_POS))))){
        ALT_alleles = unique(cigar_position_variants[cigar_position_variants$amplicon_POS == position,][['ALT']])
        
        for(ALT_allele in 1:length(ALT_alleles)){
          
          cigar_position_variants[
            cigar_position_variants$amplicon_POS == position &
              cigar_position_variants$ALT == ALT_alleles[ALT_allele],][['ALT_label']] = ALT_allele
          
        }
        
        if(sum(grepl('D', ALT_alleles)) > 0 & sum(grepl('I', ALT_alleles)) == 0){
          
          REF = gsub('D=', '', ALT_alleles[which.max(nchar(ALT_alleles))])
          
          ALT_alleles[which.max(nchar(ALT_alleles))] = '*'
          
          ALT_alleles = gsub('D=', '', ALT_alleles)
          
        }else if(sum(grepl('I', ALT_alleles)) > 0 & sum(grepl('D', ALT_alleles)) == 0){
          
          REF = as.character(Biostrings::substr(reference_genome[[unique(cigar_position_variants$CHROM)]],
                                                start = unique(cigar_position_variants[cigar_position_variants$amplicon_POS == position,][['POS']]),
                                                stop = unique(cigar_position_variants[cigar_position_variants$amplicon_POS == position,][['POS']])))
          
          ALT_alleles[which(grepl('I=', ALT_alleles))] = paste0(REF, gsub('I=', '', ALT_alleles[which(grepl('I=', ALT_alleles))]))
          
        }else if(sum(grepl('I', ALT_alleles)) == 0 & sum(grepl('D', ALT_alleles)) == 0){
          
          REF = Biostrings::substr(reference_genome[[unique(cigar_position_variants$CHROM)]],
                                   start = unique(cigar_position_variants[cigar_position_variants$amplicon_POS == position,][['POS']]),
                                   stop = unique(cigar_position_variants[cigar_position_variants$amplicon_POS == position,][['POS']]))
          
        }else if(sum(grepl('I', ALT_alleles)) > 0 & sum(grepl('D', ALT_alleles)) > 0){
          
          ALT_del = ALT_alleles[grepl('D=', ALT_alleles)]
          
          REF = gsub('D=', '', ALT_del[which.max(nchar(ALT_del))])
          
          ALT_alleles[ALT_alleles == paste0('D=', REF)] = '*'
          
          ALT_alleles = gsub('D=', '', ALT_alleles)
          
          ALT_alleles[which(grepl('I=', ALT_alleles))] = paste0(REF, gsub('I=', '', ALT_alleles[which(grepl('I=', ALT_alleles))]))
          
        }
        
        amplicon_loci_table = rbind(amplicon_loci_table,
                                    data.frame(CHROM = unique(cigar_position_variants$CHROM),
                                               POS = unique(cigar_position_variants[cigar_position_variants$amplicon_POS == position,][['POS']]),
                                               ID = NA,
                                               REF = REF,
                                               ALT = paste(ALT_alleles, collapse = ','),
                                               QUAL = NA,
                                               FILTER = NA,
                                               INFO = NA, 
                                               FORMAT = NA))
        
      }
      
      amplicon_gt = NULL
      
      # Sample_id = rownames(gt_ampseq)[3]
      
      for(Sample_id in rownames(gt_ampseq)){
        
        Sample_amplicon_gt = gt_ampseq[Sample_id, amplicon]
        
        if(Sample_id %in% monoclonals){
          
          Sample_amplicon_gt = gsub('_.+', '', Sample_amplicon_gt)
          
          if(is.na(Sample_amplicon_gt)){
            
            sample_vcf_gt = NULL
            
            for(position in amplicon_loci_table$POS){
              sample_vcf_gt = c(sample_vcf_gt,
                                paste0('./.:', 
                                       paste(rep(0,
                                                 length(unlist(str_split(amplicon_loci_table[amplicon_loci_table$POS == position,][['ALT']], ','))) + 1
                                       ), collapse = ','), ':0'))
            }
            
          }else if(grepl('\\.', Sample_amplicon_gt)){
            
            sample_vcf_gt = NULL
            
            for(position in amplicon_loci_table$POS){
              sample_vcf_gt = c(sample_vcf_gt,
                                paste0('0/0:', 
                                       paste(c(gsub('^\\.:', '', Sample_amplicon_gt),
                                               rep(0,
                                                   length(unlist(str_split(amplicon_loci_table[amplicon_loci_table$POS == position,][['ALT']], ',')))
                                               )), collapse = ','), ':', gsub('^\\.:', '', Sample_amplicon_gt)))
            }
            
          }else{
            
            sample_positions = data.frame(amplicon_POS = unlist(str_extract_all(gsub(':\\d+', '',Sample_amplicon_gt), '\\d+')),
                                          POS = as.integer(unlist(str_extract_all(gsub(':\\d+', '',Sample_amplicon_gt), '\\d+'))) + 
                                            markers[markers$amplicon == amplicon,][['start']] - 1,
                                          ALT = unlist(str_extract_all(gsub(':\\d+', '',Sample_amplicon_gt), '[DI]?=?[ATGC]+')),
                                          read_depth = gsub('^.+:', '',Sample_amplicon_gt)
            )
            
            sample_vcf_gt = NULL
            
            for(position in amplicon_loci_table$POS){
              
              if(position %in% sample_positions$POS){
                
                GT = cigar_position_variants[
                  cigar_position_variants[['cigar']] == gsub(':\\d+', '', Sample_amplicon_gt) &
                    cigar_position_variants[['POS']] == sample_positions[sample_positions$POS == position,][['POS']] &
                    cigar_position_variants[['ALT']] == sample_positions[sample_positions$POS == position,][['ALT']],][['ALT_label']]
                
                if(length(GT) == 0){
                  GT = '.'
                  print('Cigar string not in Cigar table')
                }
                
                AD = rep(0,
                         length(unlist(str_split(amplicon_loci_table[amplicon_loci_table$POS == position,][['ALT']], ','))) + 1)
                
                AD[as.integer(GT) + 1] = gsub('^.+:', '', Sample_amplicon_gt)
                
                
                sample_vcf_gt = c(sample_vcf_gt,
                                  paste0(GT, '/', GT,':', 
                                         paste(AD, collapse = ','), ':', gsub('^.+:', '', Sample_amplicon_gt)))
                
              }else{
                
                sample_vcf_gt = c(sample_vcf_gt,
                                  paste0('0/0:', 
                                         paste(c(gsub('^.+:', '', Sample_amplicon_gt),
                                                 rep(0,
                                                     length(unlist(str_split(amplicon_loci_table[amplicon_loci_table$POS == position,][['ALT']], ',')))
                                                 )), collapse = ','), ':', gsub('^.+:', '', Sample_amplicon_gt)))
                
              }
              
            }
            
          }
          
          sample_vcf_gt = matrix(sample_vcf_gt,
                                 nrow = nrow(amplicon_loci_table),
                                 ncol = 1,
                                 dimnames = list(paste(amplicon_loci_table$CHROM, amplicon_loci_table$POS, sep = '_'),
                                                 Sample_id))
          
        }else if(Sample_id %in% polyclonals){
          
          nclones = ampseq_object@metadata[ampseq_object@metadata$Sample_id == Sample_id,][['max_nAlleles']]
          
          clones_amplicon_gt = unlist(strsplit(Sample_amplicon_gt, '_'))
          
          if(length(clones_amplicon_gt) < nclones){
            
            clones_amplicon_gt = c(clones_amplicon_gt, rep(clones_amplicon_gt[1], nclones - length(clones_amplicon_gt)))
            
          }
          
          sample_vcf_gt = NULL
          
          for(clone in 1:nclones){
            
            Clone_amplicon_gt = clones_amplicon_gt[clone]
            
            if(grepl('\\.', Clone_amplicon_gt)){
              
              clone_vcf_gt = NULL
              
              for(position in amplicon_loci_table$POS){
                clone_vcf_gt = c(clone_vcf_gt,
                                 paste0('0:', 
                                        paste(c(gsub('^\\d+[DI]?=?[ATGC]+:', '', Clone_amplicon_gt),
                                                rep(0,
                                                    length(unlist(str_split(amplicon_loci_table[amplicon_loci_table$POS == position,][['ALT']], ',')))
                                                )), collapse = ','), ':', gsub('^\\d+[DI]?=?[ATGC]+:', '', Clone_amplicon_gt)))
              }
              
            }else{
              
              sample_positions = data.frame(amplicon_POS = unlist(str_extract_all(gsub(':\\d+', '',Clone_amplicon_gt), '\\d+')),
                                            POS = as.integer(unlist(str_extract_all(gsub(':\\d+', '',Clone_amplicon_gt), '\\d+'))) + 
                                              markers[markers$amplicon == amplicon,][['start']] - 1,
                                            ALT = unlist(str_extract_all(gsub(':\\d+', '',Clone_amplicon_gt), '[DI]?=?[ATGC]+')),
                                            read_depth = gsub('^.+:', '',Clone_amplicon_gt)
              )
              
              clone_vcf_gt = NULL
              
              for(position in amplicon_loci_table$POS){
                
                if(position %in% sample_positions$POS){
                  
                  GT = cigar_position_variants[
                    cigar_position_variants[['cigar']] == gsub(':\\d+', '', Clone_amplicon_gt) &
                      cigar_position_variants[['POS']] == sample_positions[sample_positions$POS == position,][['POS']] &
                      cigar_position_variants[['ALT']] == sample_positions[sample_positions$POS == position,][['ALT']],][['ALT_label']]
                  
                  AD = rep(0,
                           length(unlist(str_split(amplicon_loci_table[amplicon_loci_table$POS == position,][['ALT']], ','))) + 1)
                  
                  AD[as.integer(GT) + 1] = gsub('^\\d+[DI]?=?[ATGC]+:', '', Clone_amplicon_gt)
                  
                  
                  clone_vcf_gt = c(clone_vcf_gt,
                                   paste0(GT,':', 
                                          paste(AD, collapse = ','), ':', gsub('^\\d+[DI]?=?[ATGC]+:', '', Clone_amplicon_gt)))
                  
                }else{
                  
                  clone_vcf_gt = c(clone_vcf_gt,
                                   paste0('0:', 
                                          paste(c(gsub('^\\d+[DI]?=?[ATGC]+:', '', Clone_amplicon_gt),
                                                  rep(0,
                                                      length(unlist(str_split(amplicon_loci_table[amplicon_loci_table$POS == position,][['ALT']], ',')))
                                                  )), collapse = ','), ':', gsub('^\\d+[DI]?=?[ATGC]+:', '', Clone_amplicon_gt)))
                  
                }
                
              }
              
            }
            
            clone_vcf_gt = matrix(clone_vcf_gt,
                                  nrow = nrow(amplicon_loci_table),
                                  ncol = 1,
                                  dimnames = list(paste(amplicon_loci_table$CHROM, amplicon_loci_table$POS, sep = '_'),
                                                  paste0(Sample_id, '_C', clone)))
            
            sample_vcf_gt = cbind(sample_vcf_gt, clone_vcf_gt)
            
          }
          
        }else{
          sample_vcf_gt = NULL
          
          print(paste0('Sample_ID ', Sample_id, ' not in monoclonal or polyclonal list,\nthis sample will be omited in the output'))
          
        }
        
        amplicon_gt = cbind(amplicon_gt, sample_vcf_gt)
        
        if(sum(duplicated(colnames(amplicon_gt))) > 0){
          stop('Duplicated ID')
        }
        
      }
      
      loci_table = rbind(loci_table, amplicon_loci_table)
      gt = rbind(gt, amplicon_gt)
      
    }
    
  }
  
  loci_table[['ID']] = '.'
  loci_table[['INFO']] = '.'
  loci_table[['QUAL']] = '.'
  loci_table[['FILTER']] = '.'
  loci_table[['FORMAT']] = 'GT:AD:DP'
  
  VCF_object = cbind(loci_table, gt)
  
  return(VCF_object)
  
}

# Functions to filter or mask alleles, loci, or samples on the ampseq object----

## filter_samples----

filter_samples = function(obj, v, update_cigars = TRUE){
  
  if(class(obj) == 'ampseq'){
    obj2 = obj
    obj2@gt = obj@gt[v,]
    obj2@metadata = obj@metadata[v,]
    
    if(update_cigars){
      
      cigars_gt = get_cigar_alleles(obj2, 
                                    from = 'gt',
                                    as = 'vector')
      
      cigars_asvtab = get_cigar_alleles(obj2, 
                                        from = 'asv_table',
                                        as = 'vector',
                                        col = 'CIGAR_masked')
      
      if(sum(!(cigars_gt %in% cigars_asvtab)) > 0 & 
         sum(!(cigars_asvtab %in% cigars_gt)) > 0){
        
        cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                   paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                   '\n\n',
                   'Cigar strings in asv_table that are not present in gt matrix are:\n',
                   paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                   '\n'
        ))
        
        stop('There are cigar strings that are different in the gt and the asv_table')

        #cigar_strings_to_remove = cigars_asvtab[!(cigars_asvtab %in% cigars_gt)]
                
        #for(cigar_string_to_remove in cigar_strings_to_remove){

        #  Amplicon = gsub(';.+$', '', cigar_string_to_remove)
        #  CIGAR = gsub('^.+;', '', cigar_string_to_remove)

        #  obj2@asv_table = 
        #    obj2@asv_table[!(obj2@asv_table[['Amplicon']] == Amplicon &
        #                       obj2@asv_table[['CIGAR_masked']] == CIGAR),]

        #  obj2@asv_seqs = obj2@asv_seqs[
        #    names(obj2@asv_seqs) %in% obj2@asv_table$hapid]

        #  obj2@asv_table$hapid = paste0('ASV', 1:nrow(obj2@asv_table))
        #  names(obj2@asv_seqs) = obj2@asv_table$hapid

        #}
        
      }else if(sum(!(cigars_gt %in% cigars_asvtab)) > 0){
        
        cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                   paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                   '\n'))
        
        stop('There are cigar strings in the gt that are not present in the asv_table')
        
      }else if(sum(!(cigars_asvtab %in% cigars_gt)) > 0){
        
        cat(paste0('The following cigar strings in asv_table are not longer present in gt matrix:\n',
                   paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                   '\n'))
        
        cigar_strings_to_remove = cigars_asvtab[!(cigars_asvtab %in% cigars_gt)]
        
        for(cigar_string_to_remove in unique(cigar_strings_to_remove)){
          
          Amplicon = gsub(';.+$', '', cigar_string_to_remove)
          CIGAR = gsub('^.+;', '', cigar_string_to_remove)
          
          # obj2@asv_table[(obj2@asv_table[['Amplicon']] == Amplicon &
          #                    obj2@asv_table[['CIGAR_masked']] == CIGAR),]
          
          obj2@asv_table = 
            obj2@asv_table[!(obj2@asv_table[['Amplicon']] == Amplicon &
                               obj2@asv_table[['CIGAR_masked']] == CIGAR),]
          
          obj2@asv_seqs = obj2@asv_seqs[
            names(obj2@asv_seqs) %in% obj2@asv_table$hapid]
          
          obj2@asv_seqs_masked = obj2@asv_seqs_masked[
            names(obj2@asv_seqs_masked) %in% obj2@asv_table$hapid]
          
          obj2@asv_table$hapid = paste0('ASV', 1:nrow(obj2@asv_table))
          names(obj2@asv_seqs) = obj2@asv_table$hapid
          names(obj2@asv_seqs_masked) = obj2@asv_table$hapid
          
        }
        
        
      }else{
        print('cigar strings are consistent between gt and asv_table')
      }
      
      consistency_between_gt_and_asvtab(obj2)
      
    }
    
    
    
    
  }else if(class(obj) == 'rGenome'){
    
    if(is.logical(v)){
      if(sum(v) > 1){
        
        gt = (obj@gt[,v])
        
      }else if(sum(v) == 1){
        
        gt = matrix(obj@gt[,v], 
                    ncol = 1,
                    nrow = nrow(obj@gt),
                    dimnames = list(rownames(obj@gt),
                                    obj@metadata$Sample_id[v]))
        
      }else if(sum(v) == 0){
        
        stop('Any sample meet the filtering criteria')
        
      }
    }
    
    
    obj2 = rGenome(gt = gt,
                   loci_table = obj@loci_table,
                   metadata = obj@metadata[v,])
    
  }else(
    stop('Object should be of class ampseq or rGenome')
  )
  
  return(obj2)
}

## filter_loci----

filter_loci = function(obj, v, update_cigars = TRUE){
  
  if(class(obj) == 'ampseq'){
    
    obj2 = obj
    
    if(is.logical(v)){
      obj2@gt = obj@gt[, v]
      obj2@markers = obj@markers[v, ]
      
      if(!is.null(obj2@loci_performance)){
        obj2@loci_performance = obj@loci_performance[v, ]
      }
      
    }else if(is.character(v)){
      obj2@gt = obj@gt[, colnames(obj@gt) %in% v]
      obj2@markers = obj@markers[obj@markers$amplicon %in% v, ]
      
      if(!is.null(obj2@loci_performance)){
        obj2@loci_performance = obj@loci_performance[obj@loci_performance$loci %in% v, ]
      }
      
    }
    
    obj2@markers[["distance"]] = Inf
    
    for(chromosome in levels(as.factor(obj2@markers[["chromosome"]]))){
      for(amplicon in 1:(nrow(obj2@markers[obj2@markers[["chromosome"]] == chromosome,])-1)){
        obj2@markers[
          obj2@markers[["chromosome"]] == chromosome,
        ][amplicon, "distance"] = 
          obj2@markers[obj2@markers[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - 
          obj2@markers[obj2@markers[["chromosome"]] == chromosome,][amplicon, "pos"]
      }
    }
    
    if(update_cigars){
      
      cigars_gt = get_cigar_alleles(obj2, 
                                    from = 'gt',
                                    as = 'vector')
      
      cigars_asvtab = get_cigar_alleles(obj2, 
                                        from = 'asv_table',
                                        as = 'vector',
                                        col = 'CIGAR_masked')
      
      if(sum(!(cigars_gt %in% cigars_asvtab)) > 0 & 
         sum(!(cigars_asvtab %in% cigars_gt)) > 0){
        
        cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                   paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                   '\n\n',
                   'Cigar strings in asv_table that are not present in gt matrix are:\n',
                   paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                   '\n'
        ))
        
        stop('There are cigar strings that are different in the gt and the asv_table')
        
      }else if(sum(!(cigars_gt %in% cigars_asvtab)) > 0){
        
        cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                   paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                   '\n'))
        
        stop('There are cigar strings in the gt that are not present in the asv_table')
        
      }else if(sum(!(cigars_asvtab %in% cigars_gt)) > 0){
        
        cat(paste0('The following cigar strings in asv_table are not longer present in gt matrix:\n',
                   paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                   '\n'))
        
        cigar_strings_to_remove = cigars_asvtab[!(cigars_asvtab %in% cigars_gt)]
        
        for(cigar_string_to_remove in unique(cigar_strings_to_remove)){
          
          Amplicon = gsub(';.+$', '', cigar_string_to_remove)
          CIGAR = gsub('^.+;', '', cigar_string_to_remove)
          
          obj2@asv_table = 
            obj2@asv_table[!(obj2@asv_table[['Amplicon']] == Amplicon &
                               obj2@asv_table[['CIGAR_masked']] == CIGAR),]
          
          obj2@asv_seqs = obj2@asv_seqs[
            names(obj2@asv_seqs) %in% obj2@asv_table$hapid]
          
          obj2@asv_seqs_masked = obj2@asv_seqs_masked[
            names(obj2@asv_seqs_masked) %in% obj2@asv_table$hapid]
          
          obj2@asv_table$hapid = paste0('ASV', 1:nrow(obj2@asv_table))
          names(obj2@asv_seqs) = obj2@asv_table$hapid
          names(obj2@asv_seqs_masked) = obj2@asv_table$hapid
          
        }
        
      }else{
        print('cigar strings are consistent between gt and asv_table')
      }
      
      consistency_between_gt_and_asvtab(obj2)
      
    }
    
    
  }else if(class(obj) == 'rGenome'){
    
    obj2 = obj
    obj2@gt = obj2@gt[v,]
    obj2@loci_table = obj2@loci_table[v,]
    
    if(is.null(nrow(obj2@gt)) & ncol(obj@gt) > 1){
      
      obj2@gt = matrix(obj2@gt, nrow = 1, ncol = length(obj2@gt),
                       dimnames = list(
                         rownames(obj2@loci_table),
                         names(obj2@gt)
                       ))
      
    }else if(is.null(nrow(obj2@gt)) & ncol(obj@gt) == 1){
      
      obj2@gt = matrix(obj2@gt, ncol = 1, nrow = length(obj2@gt),
                       dimnames = list(
                         rownames(obj2@loci_table),
                         colnames(obj@gt)
                       ))
      
    }
    
  }else{
    stop('Object should be of calss ampseq or rGenome')
  }
  
  return(obj2)
}

## Mask alternative alleles----
setGeneric("mask_alt_alleles", function(obj = NULL, ref_fasta = NULL, mask_formula = "dVSITES_ij > 0.3", homopolymer_length = 5) standardGeneric("mask_alt_alleles"))

setMethod("mask_alt_alleles", signature(obj = "ampseq"),
          
          function(obj = NULL,  ref_fasta = NULL, mask_formula = "dVSITES_ij > 0.3", homopolymer_length = 5){
            
            gt = obj@gt
            mhaps = obj@markers
            asv_table = obj@asv_table
            
            mhaps$homopolymer_regions = NA
            mhaps$homopolymer_anchors = NA
            
            mhaps$STR_regions = NA
            
            ref_sequences = readDNAStringSet(ref_fasta)
            
            if(sum(grepl('(/|-|:)', names(ref_sequences))) > 0){
              
              names(ref_sequences) = gsub('(/|-|:)', '_', names(ref_sequences))
              print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the names of the fasta sequences.')
              
            }
            
            homopolymer_pattern = '(A{length,}|T{length,}|G{length,}|C{length,})'
            
            homopolymer_pattern = gsub('length', homopolymer_length, homopolymer_pattern)
            
            for(mhap in mhaps$amplicon){
              
              # Identification of Homopolymers
              
              homopolymers = unlist(str_extract_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern))
              
              if(length(homopolymers) > 0){
                
                homopolymers_location = str_locate_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern)
                
                mhaps[mhaps$amplicon == mhap, ][['homopolymer_regions']] = paste(paste(homopolymers,
                                                                                       paste(homopolymers_location[[1]][,'start'],
                                                                                             homopolymers_location[[1]][,'end'], sep = '-'), sep = ':'),
                                                                                 collapse = ',')
                
                homopolymer_starts = homopolymers_location[[1]][,'start']
                
                anchor_position = NULL
                for(homopolymer in 1:length(homopolymers_location[[1]][,'start'])){
                  
                  homopolymer_start = homopolymers_location[[1]][,'start'][homopolymer]
                  
                  if(homopolymer_start > 1){
                    
                    anchor_position = c(anchor_position, 
                                        as.character(Biostrings::substr(ref_sequences[[mhap]], 
                                                                        start = homopolymer_start - 1, 
                                                                        stop = homopolymer_start - 1)))
                  }else{
                    anchor_position = c(anchor_position, '')
                  }
                  
                }
                
                if(!is.null(anchor_position)){
                  homopolymer_anchors = paste0(homopolymers_location[[1]][,'start'] - 1, anchor_position)
                  mhaps[mhaps$amplicon == mhap, ][['homopolymer_anchors']] = paste(homopolymer_anchors, collapse = ',')  
                }
                
              }
              
              # Identification of STRs in flanking regions
              alleles = unique(unlist(strsplit(gsub(':\\d+', '',gt[,mhap]), '_')))
              alleles = alleles[!is.na(alleles)]
              
              if(sum(grepl('=', alleles)) > 0){
                
                insertions = unique(unlist(str_extract_all(alleles, '\\d+I=[ATCG]+')))
                deletions = unique(unlist(str_extract_all(alleles, '\\d+D=[ATCG]+')))
                
                if(length(insertions) > 0){
                  insertions_at_3prime = 
                    as.integer(gsub('I=[ATCG]+', '', insertions)) == mhaps[mhaps$amplicon == mhap, ][['length']]
                }else{
                  
                  insertions_at_3prime = NULL
                  
                }
                
                if(length(deletions) > 0){
                  
                  deletions_at_3prime = sapply(deletions, function(deletion){
                    grepl(paste0(gsub('\\d+D=', '', deletion), '$'),
                          as.character(ref_sequences[[mhap]]))},
                    simplify = T)
                  
                }else{
                  
                  deletions_at_3prime = NULL
                  
                }
                
                STR_at_3prime_start = NULL
                STR_at_3prime = NULL
                
                if(sum(insertions_at_3prime) > 0 | sum(deletions_at_3prime) > 0){
                  
                  insertions_at_3prime = insertions[insertions_at_3prime]
                  
                  insertions_at_3prime = gsub('\\d+I=.', '',insertions_at_3prime)
                  
                  deletions_at_3prime = deletions[deletions_at_3prime]
                  
                  deletions_at_3prime = gsub('\\d+D=', '',deletions_at_3prime)
                  
                  indels_at_3prime = unique(c(insertions_at_3prime, deletions_at_3prime))
                  
                  for(indel_at_3prime in indels_at_3prime){
                    
                    if(grepl(paste0(indel_at_3prime, '$'), as.character(ref_sequences[[mhap]]))){
                      
                      temp_ref_sequence = as.character(ref_sequences[[mhap]])
                      
                      while(grepl(paste0(indel_at_3prime, '$'), temp_ref_sequence)){
                        
                        temp_ref_sequence = gsub(paste0(indel_at_3prime, '$'), '', temp_ref_sequence)
                        
                      }
                      
                      STR_at_3prime_start = c(STR_at_3prime_start, nchar(temp_ref_sequence) + 1)
                      STR_at_3prime = c(STR_at_3prime, indel_at_3prime)
                      
                    }
                    
                  }
                  
                  if(!is.null(STR_at_3prime_start)){
                    STR_at_3prime = STR_at_3prime[which.min(STR_at_3prime_start)]
                    STR_at_3prime_start = min(STR_at_3prime_start)
                    
                    if(is.infinite(STR_at_3prime_start)){stop()}
                    
                    STR_at_3prime_end = mhaps[mhaps$amplicon == mhap, ][['length']]
                  }
                  
                }
                
                
                
                if(length(insertions) > 0){
                  insertions_at_5prime = 
                    as.integer(gsub('I=[ATCG]+', '', insertions)) == 0
                }else{
                  
                  insertions_at_5prime = NULL
                  
                }
                
                if(length(deletions) > 0){
                  
                  deletions_at_5primet = as.integer(gsub('D=[ATCG]+', '', deletions)) == 1
                  
                }else{
                  
                  deletions_at_5primet = NULL
                  
                }
                
                
                STR_at_5prime_end = NULL
                STR_at_5prime = NULL
                
                if(sum(insertions_at_5prime) > 0 | sum(deletions_at_5primet) > 0){
                  
                  insertions_at_5prime = insertions[insertions_at_5prime]
                  
                  insertions_at_5prime = gsub('\\d+I=', '',insertions_at_5prime)
                  
                  deletions_at_5primet = deletions[deletions_at_5primet]
                  
                  deletions_at_5primet = gsub('\\d+D=', '',deletions_at_5primet)
                  
                  indels_at_5prime = unique(c(insertions_at_5prime, deletions_at_5primet))
                  
                  for(indel_at_5prime in indels_at_5prime){
                    
                    if(grepl(paste0('^', indel_at_5prime), as.character(ref_sequences[[mhap]]))){
                      
                      temp_ref_sequence = as.character(ref_sequences[[mhap]])
                      
                      while(grepl(paste0('^', indel_at_5prime), temp_ref_sequence)){
                        
                        temp_ref_sequence = gsub(paste0('^', indel_at_5prime), '', temp_ref_sequence)
                        
                      }
                      
                      STR_at_5prime_end = c(STR_at_5prime_end, 
                                            nchar(as.character(ref_sequences[[mhap]])) -
                                              nchar(temp_ref_sequence))
                      STR_at_5prime = c(STR_at_5prime, indel_at_5prime)
                      
                    }
                    
                  }
                  
                  if(!is.null(STR_at_5prime_end)){
                    
                    STR_at_5prime = STR_at_5prime[which.max(STR_at_5prime_end)]
                    STR_at_5prime_end = max(STR_at_5prime_end)
                    STR_at_5prime_start = 0
                  }
                  
                }
                
                if(!is.null(STR_at_5prime_end) & !is.null(STR_at_3prime_start)){
                  
                  mhaps[mhaps$amplicon == mhap, ][['STR_regions']] = 
                    paste(c(paste0(STR_at_5prime, ':', STR_at_5prime_start, '-', STR_at_5prime_end),
                            paste0(STR_at_3prime, ':', STR_at_3prime_start, '-', STR_at_3prime_end)), 
                          collapse = ',')
                  
                }else if(!is.null(STR_at_5prime_end) & is.null(STR_at_3prime_start)){
                  
                  mhaps[mhaps$amplicon == mhap, ][['STR_regions']] = paste0(STR_at_5prime, ':', STR_at_5prime_start, '-', STR_at_5prime_end)
                  
                }else if(is.null(STR_at_5prime_end) & !is.null(STR_at_3prime_start)){
                  
                  mhaps[mhaps$amplicon == mhap, ][['STR_regions']] = paste0(STR_at_3prime, ':', STR_at_3prime_start, '-', STR_at_3prime_end)
                  
                }
                
              }
              
              
            }
            
            alt = sapply(colnames(gt), function(mhap){
              alt = unique(gsub(':\\d+', '',unlist(strsplit(gt[,mhap], '_'))))
              
              alt = paste(alt[!is.na(alt) & alt != '.'], collapse = ',')
            })
            
            gt_masked = gt
            
            gt = gsub(':\\d+', '',gt)
            
            # Heterozygous positions
            HetPos = matrix(grepl('_', gt), ncol = ncol(gt), nrow = nrow(gt))
            
            # Check if formula is correct
            
            if(grepl("(h_ij|h_ijminor|p_ij|P_ij|H_ij|H_ijminor|nVSITES_ij|dVSITES_ij|nSNVs_ij|dSNVs_ij|nINDELs_ij|dINDELs_ij|flanking_INDEL|SNV_in_homopolymer|INDEL_in_homopolymer|Variant_in_STR|bimera)(<|>|!|=)+", mask_formula)){
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
            
            if(grepl("Variant_in_STR ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "Variant_in_STR (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "Variant_in_STR (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
                mask_formula = gsub("Variant_in_STR ", "ASVs_attributes_table_temp[['Variant_in_STR']] ", mask_formula)
              }else{
                stop("Filter Variant_in_STR is been called but there are spelling issues in this part of the formula")
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
              stop("Execution halted, revise mask_filter argument.\nPossible filters are:\nh_ij, h_ijminor, p_ij, P_ij, H_ij, H_ijminor, nVSITES_ij, dVSITES_ij, nSNVs_ij, dSNVs_ij, nINDELs_ij, dINDELs_ij, flanking_INDEL, SNV_in_homopolymer, INDEL_in_homopolymer, Variant_in_STR, bimera")
            }
            
            #mhap = which(colnames(gt) == 'CCR4_3')
            
            for(mhap in 1:ncol(gt)){
              
              #if(colnames(gt)[mhap] == 'pvmdr1_4'){stop('1')}
              
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
                                      'SNV'}else if(grepl('I', polymorphism)){
                                        'Insertion'
                                      }else if(grepl('D', polymorphism)){
                                        'Deletion'
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
                                  
                                  # CODE NEEDED FOR MASKING SNVs THAT ARE INSIDE OF INDELs, THIS PART IS STILL IN DEVELOPMENT
                                  # if('Insertion' %in% type_polymorphism & !is.na(mhaps[mhap,][['homopolymer_regions']])){
                                  #   snv_in_insertion_positions = positions[type_polymorphism == 'Insertion']
                                  #   snv_in_insertion_polymorphisms = polymorphisms[type_polymorphism == 'Insertion']
                                  #   snv_in_insertion_polymorphisms = substr(snv_in_insertion_polymorphisms, 3, 3)
                                  #   
                                  #   if(sum(snv_in_insertion_positions %in% masked_region) > 0){
                                  #     
                                  #     snv_in_insertion_polymorphisms_in_homopolymers = snv_in_insertion_polymorphisms[snv_in_insertion_positions %in% masked_region]
                                  #     snv_in_insertion_positions_in_homopolymers = snv_in_insertion_positions[snv_in_insertion_positions %in% masked_region]
                                  #     
                                  #     start = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '\\d+-'))))
                                  #     end = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '-\\d+'))))
                                  #     
                                  #     homopolymer_repeat = unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '[ATGC]+'))
                                  #     homopolymer_repeat = sapply(homopolymer_repeat, function(repeated_nucleotide){
                                  #       unique(unlist(strsplit(repeated_nucleotide, '')))
                                  #     })
                                  #     
                                  #     homopolymer_repeat = data.frame(
                                  #       start,
                                  #       end,
                                  #       homopolymer_repeat)
                                  #     
                                  #     homopolymer_repeat_motive_in_asv = NULL
                                  #     
                                  #     for(snv_in_insertion_position in snv_in_insertion_positions_in_homopolymers){
                                  #       homopolymer_repeat_motive_in_asv = c(homopolymer_repeat_motive_in_asv,
                                  #                                         homopolymer_repeat[homopolymer_repeat[["start"]] <= snv_in_insertion_position &
                                  #                                         homopolymer_repeat[["end"]] >= snv_in_insertion_position,][["start"]])
                                  #     }
                                  #     
                                  #     homopolymer_repeat_motive_in_asv = homopolymer_repeat[homopolymer_repeat$start %in% homopolymer_repeat_motive_in_asv,]
                                  #     
                                  #     SNV_in_Insertion_in_homopolymer = sum(homopolymer_repeat_motive_in_asv$homopolymer_repeat != snv_in_insertion_polymorphisms_in_homopolymers)
                                  #     
                                  #   }
                                  #   
                                  # }
                                  
                                  
                                  # SNV_in_homopolymer pattern and replacement
                                  
                                  if(SNV_in_homopolymer == 1){
                                    
                                    SNV_in_homopolymer_pattern = allele
                                    
                                    snv_positions_remined = snv_positions[!(snv_positions %in% masked_region)]
                                    snv_polymorphisms_remined = snv_polymorphisms[!(snv_positions %in% masked_region)]
                                    
                                    indel_postions_to_add = positions[type_polymorphism %in% c('Insertion', 'Deletion')]
                                    indel_polymorphisms_to_add = polymorphisms[type_polymorphism %in% c('Insertion', 'Deletion')]
                                    
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
                                  
                                  if(sum(c('Insertion', 'Deletion') %in% type_polymorphism & !is.na(mhaps[mhap,][['homopolymer_regions']])) > 0){
                                    
                                    deleted_positions = positions[type_polymorphism %in% c('Deletion')]
                                    inserted_positions = positions[type_polymorphism %in% c('Insertion')]
                                    
                                    deleted_polymorphisms = polymorphisms[type_polymorphism %in% c('Deletion')]
                                    inserted_polymorphisms = polymorphisms[type_polymorphism %in% c('Insertion')]
                                    
                                    
                                    start = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '\\d+-'))))
                                    end = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '-\\d+'))))
                                    
                                    masked_region_for_deletions = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos], end[pos], 1)})))
                                    
                                    masked_region_for_insertions = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos] - 1, end[pos], 1)})))
                                    
                                    
                                    Insertions_in_homopolymer = sum(inserted_positions %in% masked_region_for_insertions) > 0
                                    Deletions_in_homopolymer = sum(deleted_positions %in% masked_region_for_insertions) > 0
                                    
                                  }else{
                                    Insertions_in_homopolymer = FALSE
                                    Deletions_in_homopolymer = FALSE
                                  }
                                  
                                  
                                  # INDEL_in_homopolymer pattern and replacement
                                  
                                  INDEL_in_homopolymer = Insertions_in_homopolymer == 1 | Deletions_in_homopolymer == 1
                                  
                                  if(INDEL_in_homopolymer == 1){
                                    
                                    INDEL_in_homopolymer_pattern = allele
                                    
                                    inserted_positions_remined = inserted_positions[!(inserted_positions %in% masked_region_for_insertions)]
                                    inserted_polymorphisms_remined = inserted_polymorphisms[!(inserted_positions %in% masked_region_for_insertions)]
                                    
                                    deleted_positions_remined = deleted_positions[!(deleted_positions %in% masked_region_for_deletions)]
                                    deleted_polymorphisms_remined = deleted_polymorphisms[!(deleted_positions %in% masked_region_for_deletions)]
                                    
                                    snv_postions_to_add = positions[type_polymorphism == 'SNV']
                                    snv_polymorphisms_to_add = polymorphisms[type_polymorphism == 'SNV']
                                    
                                    INDEL_in_homopolymer_replacement = data.frame(positions = c(inserted_positions_remined,
                                                                                                deleted_positions_remined,
                                                                                                snv_postions_to_add),
                                                                                  polymorphisms = c(inserted_polymorphisms_remined,
                                                                                                    deleted_polymorphisms_remined,
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
                                  
                                  
                                  # Variant_in_STR
                                  
                                  if(!is.na(mhaps[mhap,][['STR_regions']])){
                                    
                                    SNV_DEL_positions = positions[type_polymorphism %in% c('Deletion', 'SNV')]
                                    INS_positions = positions[type_polymorphism %in% c('Insertion')]
                                    
                                    SNV_DEL_polymorphisms = polymorphisms[type_polymorphism %in% c('Deletion', 'SNV')]
                                    INS_polymorphisms = polymorphisms[type_polymorphism %in% c('Insertion')]
                                    
                                    start = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['STR_regions']], '\\d+-'))))
                                    end = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['STR_regions']], '-\\d+'))))
                                    
                                    masked_region_for_SNV_DEL = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos], end[pos], 1)})))
                                    
                                    masked_region_for_INS = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos] - 1, end[pos], 1)})))
                                    
                                    Insertion_in_STR = sum(positions %in% INS_positions & positions %in% masked_region_for_INS) > 0
                                    SNV_or_Deletion_in_homopolymer = sum(positions %in% masked_region_for_SNV_DEL) > 0
                                    
                                  }else{
                                    Insertion_in_STR = FALSE
                                    SNV_or_Deletion_in_homopolymer = FALSE
                                  }
                                  
                                  # Variant_in_STR pattern and replacement
                                  
                                  Variant_in_STR = Insertion_in_STR == 1 | SNV_or_Deletion_in_homopolymer == 1
                                  
                                  if(Variant_in_STR == 1){
                                    
                                    Variant_in_STR_pattern = allele
                                    
                                    INS_positions_remined = INS_positions[!(INS_positions %in% masked_region_for_INS)]
                                    INS_polymorphisms_remined = INS_polymorphisms[!(INS_positions %in% masked_region_for_INS)]
                                    
                                    SNV_DEL_positions_remined = SNV_DEL_positions[!(SNV_DEL_positions %in% masked_region_for_SNV_DEL)]
                                    SNV_DEL_polymorphisms_remined = SNV_DEL_polymorphisms[!(SNV_DEL_positions %in% masked_region_for_SNV_DEL)]
                                    
                                    Variant_in_STR_replacement = data.frame(positions = c(INS_positions_remined,
                                                                                          SNV_DEL_positions_remined),
                                                                            polymorphisms = c(INS_polymorphisms_remined,
                                                                                              SNV_DEL_polymorphisms_remined))
                                    
                                    Variant_in_STR_replacement = Variant_in_STR_replacement[order(Variant_in_STR_replacement$positions),]
                                    
                                    Variant_in_STR_replacement = gsub(' ',
                                                                      '',
                                                                      paste(apply(Variant_in_STR_replacement, 
                                                                                  1, 
                                                                                  function(pos){paste(pos, collapse = '')}), 
                                                                            collapse = ''))
                                    
                                    if(Variant_in_STR_replacement == ''){
                                      Variant_in_STR_replacement = '.'
                                    }
                                    
                                  }else{
                                    
                                    Variant_in_STR_pattern = NA
                                    Variant_in_STR_replacement = NA
                                    
                                  }
                                  
                                  
                                  # INDELs in flanking regions
                                  
                                  flanking_Insertion = as.integer(grepl(paste0(paste('(^0', mhaps[mhap,][['length']], sep = '|'),')I(=|\\.)[ATGC]+') ,allele))
                                  flanking_Deletion = as.integer(grepl(paste0(paste('(^1', mhaps[mhap,][['length']], sep = '|'),')D(=|\\.)[ATGC]+') ,allele)) # ERROR: this recognition pattern only takes into account deletions that occurs at the last position and not deletions that includes the last position
                                  
                                  flanking_INDEL = flanking_Insertion == 1 | flanking_Deletion == 1
                                  
                                  if(flanking_INDEL == 1){
                                    temp_flanking_Insertion_pattern = str_extract_all(allele, paste0(paste('(^0', paste0('[ATGC]?', mhaps[mhap,][['length']]), sep = '|'),')I(=|\\.)[ATGC]+\\d?'))[[1]]
                                    temp_flanking_Insertion_replacement = gsub(paste0(paste('(^0', mhaps[mhap,][['length']], sep = '|'),')I(=|\\.)[ATGC]+'), '', temp_flanking_Insertion_pattern)
                                    
                                    temp_flanking_Deletion_pattern = str_extract_all(allele, paste0(paste('(^1', paste0('[ATGC]?', mhaps[mhap,][['length']]), sep = '|'),')D(=|\\.)[ATGC]+\\d?'))[[1]]
                                    temp_flanking_Deletion_replacement = gsub(paste0(paste('(^1', mhaps[mhap,][['length']], sep = '|'),')D(=|\\.)[ATGC]+'), '', temp_flanking_Deletion_pattern)
                                    
                                    
                                    temp_flanking_INDEL_pattern = c(temp_flanking_Insertion_pattern,
                                                                    temp_flanking_Deletion_pattern
                                    )
                                    
                                    temp_flanking_INDEL_replacement = c(temp_flanking_Insertion_replacement,
                                                                        temp_flanking_Deletion_replacement)
                                    
                                    if('' %in% temp_flanking_INDEL_replacement){
                                      
                                      flanking_INDEL_replacement = '.'
                                      
                                    }else{
                                      
                                      flanking_INDEL_replacement = allele
                                      
                                      for(replacement in 1:length(temp_flanking_INDEL_pattern)){
                                        
                                        flanking_INDEL_replacement = gsub(temp_flanking_INDEL_pattern[replacement], temp_flanking_INDEL_replacement[replacement], flanking_INDEL_replacement)
                                      }
                                      
                                    }
                                    
                                    flanking_INDEL_pattern = allele
                                    
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
                                    
                                    mhaps[mhap,][['homopolymer_anchors']],
                                    mhaps[mhap,][['homopolymer_regions']],
                                    
                                    mhaps[mhap,][['STR_regions']],
                                    
                                    SNV_in_homopolymer,
                                    SNV_in_homopolymer_pattern,
                                    SNV_in_homopolymer_replacement,
                                    
                                    INDEL_in_homopolymer,
                                    INDEL_in_homopolymer_pattern,
                                    INDEL_in_homopolymer_replacement,
                                    
                                    Variant_in_STR,
                                    Variant_in_STR_pattern,
                                    Variant_in_STR_replacement,
                                    
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
                                                      'homopolymer_anchors',
                                                      'homopolymer_regions',
                                                      'STR_regions',
                                                      'SNV_in_homopolymer',
                                                      'SNV_in_homopolymer_pattern',
                                                      'SNV_in_homopolymer_replacement',
                                                      'INDEL_in_homopolymer',
                                                      'INDEL_in_homopolymer_pattern',
                                                      'INDEL_in_homopolymer_replacement',
                                                      'Variant_in_STR',
                                                      'Variant_in_STR_pattern',
                                                      'Variant_in_STR_replacement',
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
                
                ASVs_attributes_table_temp[['Variant_in_STR']] = as.logical(ASVs_attributes_table_temp[['Variant_in_STR']])
                
                
                #ASVs_attributes_table_temp[['flanking_INDEL']] = as.logical(as.integer(ASVs_attributes_table_temp[['flanking_INDEL']]))
                
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
                    
                    # mask_formula2 = str_extract(mask_formula,
                    #                             "ASVs_attributes_table_temp\\[\\['SNV_in_homopolymer'\\]\\] (=|!)+ (TRUE|FALSE)")
                    # 
                    # replaced_alleles = ASVs_attributes_table_temp[
                    #   eval(parse(text = mask_formula2)),][['SNV_in_homopolymer_pattern']]
                    # 
                    # 
                    # replacement_alleles = ASVs_attributes_table_temp[
                    #   eval(parse(text = mask_formula2)),][['SNV_in_homopolymer_replacement']]
                    
                    replacement_alleles = ASVs_attributes_table_temp[
                      eval(parse(text = mask_formula)),][['SNV_in_homopolymer_replacement']]
                    
                    
                    for(replaced_allele in 1:length(replaced_alleles)){
                      
                      temp_replaced_allele = replaced_alleles[replaced_allele]
                      temp_replacement_allele = replacement_alleles[replaced_allele]
                      
                      asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                  !is.na(asv_table[['Amplicon']]) &
                                  replace_na(asv_table[['CIGAR_masked']], 'NAs') == temp_replaced_allele
                                ,][['CIGAR_masked']] = temp_replacement_allele
                      
                      for(sample in 1:nrow(gt_masked)){
                        
                        # which(grepl(paste0('_',temp_replaced_allele, ':'), gt_masked[,mhap]))
                        # 
                        # if(grepl(paste0('(^|_)',temp_replaced_allele, ':'), gt_masked[sample,mhap])){
                        #   stop('Pattern in sample')
                        # }
                        
                        replaced_pattern = str_extract(gt_masked[sample,mhap], paste0('(^|_)',temp_replaced_allele, ':'))
                        
                        if(!is.na(replaced_pattern)){
                          #stop('Sample with pattern')
                          replacement_pattern = gsub(temp_replaced_allele, temp_replacement_allele, replaced_pattern)
                          
                          # if(grepl('30C62A', gsub(replaced_pattern, replacement_pattern, gt_masked[sample,mhap]))){
                          #   stop('Sample with pattern')
                          # }
                          
                          gt_masked[sample,mhap] = gsub(ifelse(grepl('_',replaced_pattern), replaced_pattern, paste0('^', replaced_pattern)), 
                                                        replacement_pattern, 
                                                        gt_masked[sample,mhap])
                          
                          
                          # gt_masked[,mhap][grepl('40I=TGGGGG\\.',gt_masked[,mhap])]
                          
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
                    
                    # sum(grepl(paste0('(^|_)',temp_replaced_allele, ':'), gt_masked[,mhap]))
                    # gt_masked[,mhap][grepl(paste0('(^|_)',temp_replaced_allele, ':'), gt_masked[,mhap])]
                    
                  }
                  
                  
                }else if(grepl('INDEL_in_homopolymer',mask_formula)){
                  
                  replaced_alleles = ASVs_attributes_table_temp[
                    eval(parse(text = mask_formula)),][['INDEL_in_homopolymer_pattern']]
                  
                  if(length(replaced_alleles) > 0){
                    
                    # mask_formula2 = str_extract(mask_formula,
                    #                             "ASVs_attributes_table_temp\\[\\['INDEL_in_homopolymer'\\]\\] (=|!)+ (TRUE|FALSE)")
                    # 
                    # replaced_alleles = ASVs_attributes_table_temp[
                    #   eval(parse(text = mask_formula2)),][['INDEL_in_homopolymer_pattern']]
                    # 
                    # 
                    # replacement_alleles = ASVs_attributes_table_temp[
                    #   eval(parse(text = mask_formula2)),][['INDEL_in_homopolymer_replacement']]
                    
                    replacement_alleles = ASVs_attributes_table_temp[
                      eval(parse(text = mask_formula)),][['INDEL_in_homopolymer_replacement']]
                    
                    
                    for(replaced_allele in 1:length(replaced_alleles)){
                      
                      temp_replaced_allele = replaced_alleles[replaced_allele]
                      temp_replacement_allele = replacement_alleles[replaced_allele]
                      
                      asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                  !is.na(asv_table[['Amplicon']]) &
                                  replace_na(asv_table[['CIGAR_masked']], 'NAs') == temp_replaced_allele
                                ,][['CIGAR_masked']] = temp_replacement_allele
                      
                      for(sample in 1:nrow(gt_masked)){
                        
                        replaced_pattern = str_extract(gt_masked[sample,mhap], paste0('(^|_)',temp_replaced_allele, ':'))
                        
                        if(!is.na(replaced_pattern)){
                          
                          replacement_pattern = gsub(temp_replaced_allele, temp_replacement_allele, replaced_pattern)
                          gt_masked[sample,mhap] = gsub(ifelse(grepl('_',replaced_pattern), replaced_pattern, paste0('^', replaced_pattern)), 
                                                        replacement_pattern, 
                                                        gt_masked[sample,mhap])
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
                  
                  
                }else if(grepl('Variant_in_STR',mask_formula)){
                  
                  replaced_alleles = ASVs_attributes_table_temp[
                    eval(parse(text = mask_formula)),][['Variant_in_STR_pattern']]
                  
                  if(length(replaced_alleles) > 0){
                    
                    replacement_alleles = ASVs_attributes_table_temp[
                      eval(parse(text = mask_formula)),][['Variant_in_STR_replacement']]
                    
                    for(replaced_allele in 1:length(replaced_alleles)){
                      
                      temp_replaced_allele = replaced_alleles[replaced_allele]
                      temp_replacement_allele = replacement_alleles[replaced_allele]
                      
                      asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                  !is.na(asv_table[['Amplicon']]) &
                                  replace_na(asv_table[['CIGAR_masked']], 'NAs') == temp_replaced_allele
                                ,][['CIGAR_masked']] = temp_replacement_allele
                      
                      for(sample in 1:nrow(gt_masked)){
                        
                        replaced_pattern = str_extract(gt_masked[sample,mhap], paste0('(^|_)',temp_replaced_allele, ':'))
                        
                        if(!is.na(replaced_pattern)){
                          
                          replacement_pattern = gsub(temp_replaced_allele, temp_replacement_allele, replaced_pattern)
                          gt_masked[sample,mhap] = gsub(ifelse(grepl('_',replaced_pattern), replaced_pattern, paste0('^', replaced_pattern)), 
                                                        replacement_pattern, 
                                                        gt_masked[sample,mhap])
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
                    
                    # mask_formula2 = str_extract(mask_formula,
                    #                             "ASVs_attributes_table_temp\\[\\['flanking_INDEL'\\]\\] (=|!)+ (TRUE|FALSE)")
                    # 
                    # replaced_alleles = ASVs_attributes_table_temp[
                    #   eval(parse(text = mask_formula2)),][['flanking_INDEL_pattern']]
                    # 
                    # 
                    # replacement_alleles = ASVs_attributes_table_temp[
                    #   eval(parse(text = mask_formula2)),][['flanking_INDEL_replacement']]
                    
                    replacement_alleles = ASVs_attributes_table_temp[
                      eval(parse(text = mask_formula)),][['flanking_INDEL_replacement']]
                    
                    replacement_alleles = replacement_alleles[!duplicated(replaced_alleles)]
                    replaced_alleles = replaced_alleles[!duplicated(replaced_alleles)]
                    
                    for(replaced_allele in 1:length(replaced_alleles)){
                      
                      temp_replaced_allele = unlist(strsplit(replaced_alleles[replaced_allele], '\\|\\|'))
                      temp_replacement_allele = unlist(strsplit(replacement_alleles[replaced_allele], '\\|\\|'))
                      
                      
                      for(i in 1:length(temp_replaced_allele)){
                        
                        replaced_pattern = temp_replaced_allele[i]
                        
                        asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                    !is.na(asv_table[['Amplicon']]) &
                                    grepl(paste0('^',replaced_pattern, '$'), asv_table[['CIGAR_masked']])
                                  ,][['CIGAR_masked']] = gsub(paste0('^',replaced_pattern, '$'),
                                                              temp_replacement_allele, 
                                                              asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                                                          !is.na(asv_table[['Amplicon']]) &
                                                                          grepl(paste0('^',replaced_pattern, '$'), asv_table[['CIGAR_masked']])
                                                                        ,][['CIGAR_masked']])
                        
                        for(sample in 1:nrow(gt_masked)){
                          
                          
                          if(grepl(paste0('(^|_)',replaced_pattern, ':'), gt_masked[sample,mhap])){
                            
                            sample_replaced_pattern = str_extract(gt_masked[sample,mhap], paste0('(^|_)',replaced_pattern, ':'))
                            sample_replacement_pattern = gsub(replaced_pattern, temp_replacement_allele[i], sample_replaced_pattern)
                            
                            gt_masked[sample,mhap] = gsub(sample_replaced_pattern, sample_replacement_pattern, gt_masked[sample,mhap])
                            
                          }
                          
                          
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
                    
                    gt_ad = sapply(gt_masked[,mhap], function(Sample_id){
                      sum(as.integer(gsub(':', '',unlist(str_extract_all(Sample_id, ':\\d+')))), na.rm = T)
                    }, simplify = T)
                    
                    for(removed_pattern in paste('(^|_)', removed_alleles, ':\\d+_?', sep = '')){
                      
                      # Mask alleles below threshold
                      gt_masked[,mhap] = gsub(removed_pattern, 'R', gt_masked[,mhap])
                      
                      gt_masked[,mhap] = gsub('^R', '', gt_masked[,mhap])
                      
                      gt_masked[,mhap] = gsub('R$', '', gt_masked[,mhap])
                      
                      gt_masked[,mhap] = gsub('R', '_', gt_masked[,mhap])
                      
                      
                      if(length(gt_masked[gt_masked[,mhap] == '' &
                                          !is.na(gt_masked[,mhap])
                                          , mhap]) > 0){
                        
                        # THINK IN POTENTIAL BIAS OF THIS STRATEGY
                        if('.' %in% asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                              !is.na(asv_table[['Amplicon']])
                                              ,][['CIGAR_masked']]){
                          
                          gt_masked[gt_masked[,mhap] == '' &
                                      !is.na(gt_masked[,mhap])
                                    , mhap] = paste0('.:', gt_ad[gt_masked[,mhap] == '' &
                                                                   !is.na(gt_masked[,mhap])])
                          
                        }else{
                          
                          gt_masked[gt_masked[,mhap] == '' &
                                      !is.na(gt_masked[,mhap])
                                    , mhap] = NA
                          
                        }
                        
                        
                      }
                      
                    }
                    
                    if(length(asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                        !is.na(asv_table[['Amplicon']]) &
                                        asv_table[['CIGAR_masked']] %in% removed_alleles
                                        ,][['CIGAR_masked']]) > 0){
                      
                      # asv_table[asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                      #             !is.na(asv_table[['Amplicon']]) &
                      #             asv_table[['CIGAR_masked']] %in% removed_alleles
                      #           ,][['CIGAR_masked']] = NA
                      
                      asv_table =
                        asv_table[!(asv_table[['Amplicon']] == unique(ASVs_attributes_table_temp[['MHap']]) &
                                      !is.na(asv_table[['Amplicon']]) &
                                      asv_table[['CIGAR_masked']] %in% removed_alleles)
                                  ,]
                      
                      obj@asv_seqs = obj@asv_seqs[names(obj@asv_seqs) %in% asv_table[['hapid']]]
                      obj@asv_seqs_masked = obj@asv_seqs_masked[names(obj@asv_seqs_masked) %in% asv_table[['hapid']]]
                      
                    }
                    
                    
                    
                  }
                  
                }
                
                
              }
              
              
              # if(sum(is.na(asv_table$hapid)) > 0){stop("NA's introduced in asv_table")}
              # 
              # if(grepl("14G18A19G21C26T:\\d+_14G18A19G21C26T:\\d+", gt_masked[,mhap])){stop("Duplicated cigar string not masked")}
              
            }
            
            obj@asv_seqs_masked = cigar_strings2fasta(obj = ampseq_object@asv_table, ref_fasta = ref_fasta, cigar_string_col = 'CIGAR_masked', amplicon_col = 'Amplicon', format = 'DNAStringSet')
            
            names(obj@asv_seqs_masked) = names(obj@asv_seqs)
            
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

get_ReadDepth_coverage = function(ampseq_object, variable = NULL, plot = TRUE){
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
  
  coverage_by_sample = coverage %>% 
    summarise(Total_Read_depth = sum(Read_depth, na.rm = T), .by = c(var, Sample_id))
  
  if(plot == TRUE){
    
    sample_order = coverage %>%
      group_by(Sample_id) %>%
      summarise(Read_depth = sum(Read_depth, na.rm = T)) %>%
      arrange(Read_depth) %>%
      dplyr::select(Sample_id) %>% unlist
    
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
      theme_bw() +
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
    
  }else{
    
    if(!is.null(variable)){
      names(coverage) = c(names(coverage)[1], variable, names(coverage)[c(-1, -2)])
      names(coverage_by_sample) = c(variable, names(coverage_by_sample)[-1])
      }else{
      
      names(coverage) = c(names(coverage)[1], 'Population', names(coverage)[c(-1, -2)])
      names(coverage_by_sample) = c('Population', names(coverage_by_sample)[-1])
      
    }
    
    return(list(coverage_by_marker_by_sample = coverage,
                coverage_by_sample = coverage_by_sample))
    
  }
  
}


## locus_amplification_rate----

locus_amplification_rate = function(ampseq_object, threshold = .65, update_loci = TRUE, strata = NULL, based_on_strata = FALSE, chr_lengths = NULL){
  
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
  
  print("Calculating loci performance...")
  loci_performance = data.frame(loci = colnames(ampseq_loci_abd_table),
                                loci_ampl_rate_Total = apply(ampseq_loci_abd_table, 2, function(x) 1 - sum(is.na(x))/length(x)))
  print(loci_performance)
  print("Finished calculating loci performance!") 
  
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
  
  
  if(!is.null(chr_lengths)){
    amplification_rate_per_locus = ggplot()+
      geom_segment(data = data.frame(x = 0, y = 1:14, xend =chr_lengths, yend = 1:14),
                   aes(x=x, y=y, xend=xend, yend=yend, size = 2), alpha = .5)+
      geom_point(data = data.frame(Position = ampseq_object@markers$pos,
                                   Amplification_rate =
                                     loci_performance$loci_ampl_rate,
                                   Chr = as.integer(gsub("(Pf3D7_|_v3|PvP01_|_v1)", "", ampseq_object@markers$chromosome))),
                 aes(x = Position, y = Chr, color = Amplification_rate),
                 pch = "|", size = 5)+
      theme_bw()+
      labs(title = "Amplification rate by locus",
           x = "Chromosome position",
           y = "Chromosome",
           color = NULL)+
      theme(legend.position = 'none',
            text = element_text(size = 18))+
      scale_y_continuous(breaks = 1:14)+
      scale_color_continuous(type = "viridis")
  }else{
    amplification_rate_per_locus = NULL
  }
  print("Updating loci based on locus amplification rate...")
  
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
        matrix(ampseq_loci_abd_table[, colnames(ampseq_loci_abd_table) %in% discarded_loci],
               nrow = nrow(ampseq_loci_abd_table),
               ncol = sum(colnames(ampseq_loci_abd_table) %in% discarded_loci),
               dimnames = list(rownames(ampseq_loci_abd_table),
                               colnames(ampseq_loci_abd_table)[colnames(ampseq_loci_abd_table) %in% discarded_loci])
               )
      
      ampseq_loci_abd_table = 
        matrix(ampseq_loci_abd_table[, colnames(ampseq_loci_abd_table) %in% keeped_loci],
               nrow = nrow(ampseq_loci_abd_table),
               ncol = sum(colnames(ampseq_loci_abd_table) %in% keeped_loci),
               dimnames = list(rownames(ampseq_loci_abd_table),
                               colnames(ampseq_loci_abd_table)[colnames(ampseq_loci_abd_table) %in% keeped_loci]))
      
      markers = ampseq_object@markers
      
      discarded_markers = markers[markers[['amplicon']] %in% discarded_loci,]
      markers = markers[markers[['amplicon']] %in% keeped_loci,]
      
      markers[["distance"]] = Inf
      
      for(chromosome in levels(as.factor(markers[["chromosome"]]))){
        for(amplicon in 1:(nrow(markers[markers[["chromosome"]] == chromosome,])-1)){
          markers[
            markers[["chromosome"]] == chromosome,][amplicon, "distance"] = 
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
      
      
      cigars_gt = get_cigar_alleles(ampseq_object, 
                                    from = 'gt',
                                    as = 'vector')
      
      cigars_asvtab = get_cigar_alleles(ampseq_object, 
                                        from = 'asv_table',
                                        as = 'vector',
                                        col = 'CIGAR_masked')
      
      if(sum(!(cigars_gt %in% cigars_asvtab)) > 0 & 
          sum(!(cigars_asvtab %in% cigars_gt)) > 0){
        
        cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                   paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                   '\n\n',
                   'Cigar strings in asv_table that are not present in gt matrix are:\n',
                   paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                   '\n'
        ))
        
        stop('There are cigar strings that are different in the gt and the asv_table')
        
      }else if(sum(!(cigars_gt %in% cigars_asvtab)) > 0){
        
        cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                   paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                   '\n'))
        
        stop('There are cigar strings in the gt that are not present in the asv_table')
        
      }else if(sum(!(cigars_asvtab %in% cigars_gt)) > 0){
        
        cat(paste0('The following cigar strings in asv_table are not longer present in gt matrix:\n',
                   paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                   '\n'))
        
        cigar_strings_to_remove = cigars_asvtab[!(cigars_asvtab %in% cigars_gt)]
        
        for(cigar_string_to_remove in unique(cigar_strings_to_remove)){
          
          Amplicon = gsub(';.+$', '', cigar_string_to_remove)
          CIGAR = gsub('^.+;', '', cigar_string_to_remove)
          
          ampseq_object@asv_table = 
            ampseq_object@asv_table[!(ampseq_object@asv_table[['Amplicon']] == Amplicon &
                               ampseq_object@asv_table[['CIGAR_masked']] == CIGAR),]
          
          ampseq_object@asv_seqs = ampseq_object@asv_seqs[
            names(ampseq_object@asv_seqs) %in% ampseq_object@asv_table$hapid]
          
          ampseq_object@asv_seqs_masked = ampseq_object@asv_seqs_masked[
            names(ampseq_object@asv_seqs_masked) %in% ampseq_object@asv_table$hapid]
          
          ampseq_object@asv_table$hapid = paste0('ASV', 1:nrow(ampseq_object@asv_table))
          names(ampseq_object@asv_seqs) = ampseq_object@asv_table$hapid
          names(ampseq_object@asv_seqs_masked) = ampseq_object@asv_table$hapid
          
        }
        
        
      }else{
        print('cigar strings are consistent between gt and asv_table')
      }
      
      consistency_between_gt_and_asvtab(ampseq_object)
      
      
      return(ampseq_object)
      
    }else if((is.null(strata) | !based_on_strata)){
      print("Deciding which loci to keep/discard...")

      discarded_loci = loci_performance[loci_performance[["loci_ampl_rate_Total"]] <= threshold,][["loci"]]
      
      if(length(discarded_loci) > 0){
        
        keeped_loci = loci_performance[loci_performance[["loci_ampl_rate_Total"]] > threshold,][["loci"]]
        
        # Ensure that these stay as dataframes

        ampseq_loci_abd_table_discarded_loci =
          matrix(ampseq_loci_abd_table[, colnames(ampseq_loci_abd_table) %in% discarded_loci],
                 nrow = nrow(ampseq_loci_abd_table),
                 ncol = sum(colnames(ampseq_loci_abd_table) %in% discarded_loci),
                 dimnames = list(rownames(ampseq_loci_abd_table),
                                 colnames(ampseq_loci_abd_table)[colnames(ampseq_loci_abd_table) %in% discarded_loci])
          )
        
        ampseq_loci_abd_table = 
          matrix(ampseq_loci_abd_table[, colnames(ampseq_loci_abd_table) %in% keeped_loci],
                 nrow = nrow(ampseq_loci_abd_table),
                 ncol = sum(colnames(ampseq_loci_abd_table) %in% keeped_loci),
                 dimnames = list(rownames(ampseq_loci_abd_table),
                                 colnames(ampseq_loci_abd_table)[colnames(ampseq_loci_abd_table) %in% keeped_loci]))
        
        print("Finished deciding loci to keep/discard!")

        markers = ampseq_object@markers
        
        discarded_markers = markers[markers[['amplicon']] %in% discarded_loci,]
        markers = markers[markers[['amplicon']] %in% keeped_loci,]
        
        markers[["distance"]] = Inf
        
        for(chromosome in levels(as.factor(markers[["chromosome"]]))){
          for(amplicon in 1:(nrow(markers[markers[["chromosome"]] == chromosome,])-1)){
            markers[markers[["chromosome"]] == chromosome,][amplicon, "distance"] = 
              markers[markers[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - 
              markers[markers[["chromosome"]] == chromosome,][amplicon, "pos"]
          }
        }

        print("Finished calculating marker distance!")

        loci_performance_complete = loci_performance
        loci_performance = loci_performance[keeped_loci,]
        
        ampseq_object@gt = ampseq_loci_abd_table
        ampseq_object@markers = markers
        ampseq_object@loci_performance = loci_performance
        ampseq_object@discarded_loci = list(gt = ampseq_loci_abd_table_discarded_loci,
                                            loci_performance = loci_performance_complete,
                                            markers = discarded_markers)
        
      }
      
      ampseq_object@plots[["all_loci_amplification_rate"]] = all_loci_performance_plot
      #ampseq_object@plots[["amplification_rate_per_locus"]] = amplification_rate_per_locus
      
      
      
      cigars_gt = get_cigar_alleles(ampseq_object, 
                                    from = 'gt',
                                    as = 'vector')
      
      cigars_asvtab = get_cigar_alleles(ampseq_object, 
                                        from = 'asv_table',
                                        as = 'vector',
                                        col = 'CIGAR_masked')
      
      if(sum(!(cigars_gt %in% cigars_asvtab)) > 0 & 
         sum(!(cigars_asvtab %in% cigars_gt)) > 0){
        
        cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                   paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                   '\n\n',
                   'Cigar strings in asv_table that are not present in gt matrix are:\n',
                   paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                   '\n'
        ))
        
        stop('There are cigar strings that are different in the gt and the asv_table')
        
      }else if(sum(!(cigars_gt %in% cigars_asvtab)) > 0){
        
        cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                   paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                   '\n'))
        
        stop('There are cigar strings in the gt that are not present in the asv_table')
        
      }else if(sum(!(cigars_asvtab %in% cigars_gt)) > 0){
        
        cat(paste0('The following cigar strings in asv_table are not longer present in gt matrix:\n',
                   paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                   '\n'))
        
        cigar_strings_to_remove = cigars_asvtab[!(cigars_asvtab %in% cigars_gt)]
        
        for(cigar_string_to_remove in unique(cigar_strings_to_remove)){
          
          Amplicon = gsub(';.+$', '', cigar_string_to_remove)
          CIGAR = gsub('^.+;', '', cigar_string_to_remove)
          
          ampseq_object@asv_table = 
            ampseq_object@asv_table[!(ampseq_object@asv_table[['Amplicon']] == Amplicon &
                                        ampseq_object@asv_table[['CIGAR_masked']] == CIGAR),]
          
          ampseq_object@asv_seqs = ampseq_object@asv_seqs[
            names(ampseq_object@asv_seqs) %in% ampseq_object@asv_table$hapid]
          
          ampseq_object@asv_seqs_masked = ampseq_object@asv_seqs_masked[
            names(ampseq_object@asv_seqs_masked) %in% ampseq_object@asv_table$hapid]
          
          ampseq_object@asv_table$hapid = paste0('ASV', 1:nrow(ampseq_object@asv_table))
          names(ampseq_object@asv_seqs) = ampseq_object@asv_table$hapid
          names(ampseq_object@asv_seqs_masked) = ampseq_object@asv_table$hapid
          
        }
        
        
      }else{
        print('cigar strings are consistent between gt and asv_table')
      }
      
      consistency_between_gt_and_asvtab(ampseq_object)
      
      
      
      return(ampseq_object)
      
    }
    
    
  }else{
    
    return(list(all_loci_performance_plot = all_loci_performance_plot,
                amplification_rate_per_locus = amplification_rate_per_locus))
    
  }
  
  
}

## sample_amplification_rate----

sample_amplification_rate = function(ampseq_object, threshold = .8, update_samples = TRUE, strata = NULL){
  
  metadata = ampseq_object@metadata
  ampseq_loci_abd_table = ampseq_object@gt
  loci_performance = ampseq_object@loci_performance

   # Check if ampseq_object is a vector with one locus due to locus performance dropout
  if(is.vector(ampseq_object@gt) && nrow(loci_performance) == 1) {
    locus_name <- rownames(loci_performance)
    ampseq_loci_abd_table = matrix(ampseq_object@gt, nrow=length(ampseq_object@gt), ncol=1, dimnames=list(names(ampseq_object@gt), locus_name))
  }
  
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
  
  if(update_samples & sum(metadata[["sample_ampl_rate"]] <= threshold) > 0){
    
    name_of_discarded_samples = rownames(ampseq_loci_abd_table)[rownames(ampseq_loci_abd_table) %in% metadata[metadata[["sample_ampl_rate"]] <= threshold ,][["Sample_id"]]]
    number_of_discarded_samples = length(name_of_discarded_samples)
    
    name_of_kept_samples = rownames(ampseq_loci_abd_table)[rownames(ampseq_loci_abd_table) %in% metadata[metadata[["sample_ampl_rate"]] > threshold ,][["Sample_id"]]]
    number_of_kept_samples = length(name_of_kept_samples)
    
    ampseq_loci_abd_table_discarded_samples = matrix(ampseq_loci_abd_table[name_of_discarded_samples,], 
                                                     nrow = number_of_discarded_samples, 
                                                     ncol = ncol(ampseq_loci_abd_table),
                                                     dimnames = list(name_of_discarded_samples,
                                                                     colnames(ampseq_loci_abd_table)))
    print("Printing ampseq_loci_abd_table_discarded_samples...")
    print(name_of_discarded_samples)
    ampseq_loci_abd_table = matrix(ampseq_loci_abd_table[name_of_kept_samples,],
                                   nrow = number_of_kept_samples,
                                   ncol = ncol(ampseq_loci_abd_table),
                                   dimnames = list(name_of_kept_samples,
                                                   colnames(ampseq_loci_abd_table)))
    print("Printing ampseq_loci_abd_table...")
    print(rownames(ampseq_loci_abd_table))
    metadata_complete = metadata
    metadata = metadata[metadata[["Sample_id"]] %in% name_of_kept_samples ,]
    metadata_discarded = metadata_complete[metadata_complete[["Sample_id"]] %in% name_of_discarded_samples ,]
    
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
                                           metadata = metadata_discarded)
    
    # ampseq_object@plots[["loci_amplification_rate"]] = loci_performance_plot
    ampseq_object@plots[["samples_amplification_rate"]] = all_samples_performance_plot
    
    
    cigars_gt = get_cigar_alleles(ampseq_object, 
                                  from = 'gt',
                                  as = 'vector')
    
    cigars_asvtab = get_cigar_alleles(ampseq_object, 
                                      from = 'asv_table',
                                      as = 'vector',
                                      col = 'CIGAR_masked')
    
    if(sum(!(cigars_gt %in% cigars_asvtab)) > 0 & 
       sum(!(cigars_asvtab %in% cigars_gt)) > 0){
      
      cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                 paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                 '\n\n',
                 'Cigar strings in asv_table that are not present in gt matrix are:\n',
                 paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                 '\n'
      ))
      
      stop('There are cigar strings that are different in the gt and the asv_table')
      
    }else if(sum(!(cigars_gt %in% cigars_asvtab)) > 0){
      
      cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                 paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                 '\n'))
      
      stop('There are cigar strings in the gt that are not present in the asv_table')
      
    }else if(sum(!(cigars_asvtab %in% cigars_gt)) > 0){
      
      cat(paste0('The following cigar strings in asv_table are not longer present in gt matrix:\n',
                 paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                 '\n'))
      
      cigar_strings_to_remove = cigars_asvtab[!(cigars_asvtab %in% cigars_gt)]
      
      for(cigar_string_to_remove in unique(cigar_strings_to_remove)){
        
        Amplicon = gsub(';.+$', '', cigar_string_to_remove)
        CIGAR = gsub('^.+;', '', cigar_string_to_remove)
        
        ampseq_object@asv_table = 
          ampseq_object@asv_table[!(ampseq_object@asv_table[['Amplicon']] == Amplicon &
                                      ampseq_object@asv_table[['CIGAR_masked']] == CIGAR),]
        
        ampseq_object@asv_seqs = ampseq_object@asv_seqs[
          names(ampseq_object@asv_seqs) %in% ampseq_object@asv_table$hapid]
        
        ampseq_object@asv_seqs_masked = ampseq_object@asv_seqs_masked[
          names(ampseq_object@asv_seqs_masked) %in% ampseq_object@asv_table$hapid]
        
        ampseq_object@asv_table$hapid = paste0('ASV', 1:nrow(ampseq_object@asv_table))
        names(ampseq_object@asv_seqs) = ampseq_object@asv_table$hapid
        names(ampseq_object@asv_seqs_masked) = ampseq_object@asv_table$hapid
        
      }
      
      
    }else{
      print('cigar strings are consistent between gt and asv_table')
    }
    
    consistency_between_gt_and_asvtab(ampseq_object)
    
    return(ampseq_object)
    
  }else if(update_samples & sum(metadata[["sample_ampl_rate"]] > threshold) > 0){
    
    
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
            mhaps$homopolymer_anchors = NA
            
            mhaps$STR_regions = NA
            
            ref_sequences = readDNAStringSet(ref_fasta)
            
            if(sum(grepl('(/|-|:)', names(ref_sequences))) > 0){
              
              names(ref_sequences) = gsub('(/|-|:)', '_', names(ref_sequences))
              print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the names of the fasta sequences.')
              
            }
            
            homopolymer_pattern = '(A{length,}|T{length,}|G{length,}|C{length,})'
            
            homopolymer_pattern = gsub('length', homopolymer_length, homopolymer_pattern)
            
            
            for(mhap in mhaps$amplicon){
              
              # Identification of Homopolymers
              
              homopolymers = unlist(str_extract_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern))
              
              if(length(homopolymers) > 0){
                
                homopolymers_location = str_locate_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern)
                
                mhaps[mhaps$amplicon == mhap, ][['homopolymer_regions']] = paste(paste(homopolymers,
                                                                                       paste(homopolymers_location[[1]][,'start'],
                                                                                             homopolymers_location[[1]][,'end'], sep = '-'), sep = ':'),
                                                                                 collapse = ',')
                
                homopolymer_starts = homopolymers_location[[1]][,'start']
                
                anchor_position = NULL
                for(homopolymer in 1:length(homopolymers_location[[1]][,'start'])){
                  
                  homopolymer_start = homopolymers_location[[1]][,'start'][homopolymer]
                  
                  if(homopolymer_start > 1){
                    
                    anchor_position = c(anchor_position, 
                                        as.character(Biostrings::substr(ref_sequences[[mhap]], 
                                                                        start = homopolymer_start - 1, 
                                                                        stop = homopolymer_start - 1)))
                  }else{
                    anchor_position = c(anchor_position, '')
                  }
                  
                }
                
                if(!is.null(anchor_position)){
                  homopolymer_anchors = paste0(homopolymers_location[[1]][,'start'] - 1, anchor_position)
                  mhaps[mhaps$amplicon == mhap, ][['homopolymer_anchors']] = paste(homopolymer_anchors, collapse = ',')  
                }
                
              }
              
              # Identification of STRs in flanking regions
              alleles = unique(unlist(strsplit(gsub(':\\d+', '',gt[,mhap]), '_')))
              alleles = alleles[!is.na(alleles)]
              
              if(sum(grepl('=', alleles)) > 0){
                
                insertions = unique(unlist(str_extract_all(alleles, '\\d+I=[ATCG]+')))
                deletions = unique(unlist(str_extract_all(alleles, '\\d+D=[ATCG]+')))
                
                if(length(insertions) > 0){
                  insertions_at_3prime = 
                    as.integer(gsub('I=[ATCG]+', '', insertions)) == mhaps[mhaps$amplicon == mhap, ][['length']]
                }else{
                  
                  insertions_at_3prime = NULL
                  
                }
                
                if(length(deletions) > 0){
                  
                  deletions_at_3prime = sapply(deletions, function(deletion){
                    grepl(paste0(gsub('\\d+D=', '', deletion), '$'),
                          as.character(ref_sequences[[mhap]]))},
                    simplify = T)
                  
                }else{
                  
                  deletions_at_3prime = NULL
                  
                }
                
                STR_at_3prime_start = NULL
                STR_at_3prime = NULL
                
                if(sum(insertions_at_3prime) > 0 | sum(deletions_at_3prime) > 0){
                  
                  insertions_at_3prime = insertions[insertions_at_3prime]
                  
                  insertions_at_3prime = gsub('\\d+I=.', '',insertions_at_3prime)
                  
                  deletions_at_3prime = deletions[deletions_at_3prime]
                  
                  deletions_at_3prime = gsub('\\d+D=', '',deletions_at_3prime)
                  
                  indels_at_3prime = unique(c(insertions_at_3prime, deletions_at_3prime))
                  
                  for(indel_at_3prime in indels_at_3prime){
                    
                    if(grepl(paste0(indel_at_3prime, '$'), as.character(ref_sequences[[mhap]]))){
                      
                      temp_ref_sequence = as.character(ref_sequences[[mhap]])
                      
                      while(grepl(paste0(indel_at_3prime, '$'), temp_ref_sequence)){
                        
                        temp_ref_sequence = gsub(paste0(indel_at_3prime, '$'), '', temp_ref_sequence)
                        
                      }
                      
                      STR_at_3prime_start = c(STR_at_3prime_start, nchar(temp_ref_sequence) + 1)
                      STR_at_3prime = c(STR_at_3prime, indel_at_3prime)
                      
                    }
                    
                  }
                  
                  if(!is.null(STR_at_3prime_start)){
                    STR_at_3prime = STR_at_3prime[which.min(STR_at_3prime_start)]
                    STR_at_3prime_start = min(STR_at_3prime_start)
                    
                    if(is.infinite(STR_at_3prime_start)){stop()}
                    
                    STR_at_3prime_end = mhaps[mhaps$amplicon == mhap, ][['length']]
                  }
                  
                }
                
                
                
                if(length(insertions) > 0){
                  insertions_at_5prime = 
                    as.integer(gsub('I=[ATCG]+', '', insertions)) == 0
                }else{
                  
                  insertions_at_5prime = NULL
                  
                }
                
                if(length(deletions) > 0){
                  
                  deletions_at_5primet = as.integer(gsub('D=[ATCG]+', '', deletions)) == 1
                  
                }else{
                  
                  deletions_at_5primet = NULL
                  
                }
                
                
                STR_at_5prime_end = NULL
                STR_at_5prime = NULL
                
                if(sum(insertions_at_5prime) > 0 | sum(deletions_at_5primet) > 0){
                  
                  insertions_at_5prime = insertions[insertions_at_5prime]
                  
                  insertions_at_5prime = gsub('\\d+I=', '',insertions_at_5prime)
                  
                  deletions_at_5primet = deletions[deletions_at_5primet]
                  
                  deletions_at_5primet = gsub('\\d+D=', '',deletions_at_5primet)
                  
                  indels_at_5prime = unique(c(insertions_at_5prime, deletions_at_5primet))
                  
                  for(indel_at_5prime in indels_at_5prime){
                    
                    if(grepl(paste0('^', indel_at_5prime), as.character(ref_sequences[[mhap]]))){
                      
                      temp_ref_sequence = as.character(ref_sequences[[mhap]])
                      
                      while(grepl(paste0('^', indel_at_5prime), temp_ref_sequence)){
                        
                        temp_ref_sequence = gsub(paste0('^', indel_at_5prime), '', temp_ref_sequence)
                        
                      }
                      
                      STR_at_5prime_end = c(STR_at_5prime_end, 
                                            nchar(as.character(ref_sequences[[mhap]])) -
                                              nchar(temp_ref_sequence))
                      STR_at_5prime = c(STR_at_5prime, indel_at_5prime)
                      
                    }
                    
                  }
                  
                  if(!is.null(STR_at_5prime_end)){
                    
                    STR_at_5prime = STR_at_5prime[which.max(STR_at_5prime_end)]
                    STR_at_5prime_end = max(STR_at_5prime_end)
                    STR_at_5prime_start = 0
                  }
                  
                }
                
                if(!is.null(STR_at_5prime_end) & !is.null(STR_at_3prime_start)){
                  
                  mhaps[mhaps$amplicon == mhap, ][['STR_regions']] = 
                    paste(c(paste0(STR_at_5prime, ':', STR_at_5prime_start, '-', STR_at_5prime_end),
                            paste0(STR_at_3prime, ':', STR_at_3prime_start, '-', STR_at_3prime_end)), 
                          collapse = ',')
                  
                }else if(!is.null(STR_at_5prime_end) & is.null(STR_at_3prime_start)){
                  
                  mhaps[mhaps$amplicon == mhap, ][['STR_regions']] = paste0(STR_at_5prime, ':', STR_at_5prime_start, '-', STR_at_5prime_end)
                  
                }else if(is.null(STR_at_5prime_end) & !is.null(STR_at_3prime_start)){
                  
                  mhaps[mhaps$amplicon == mhap, ][['STR_regions']] = paste0(STR_at_3prime, ':', STR_at_3prime_start, '-', STR_at_3prime_end)
                  
                }
                
              }
              
              
            }
            
            
            
            
            alt = sapply(colnames(gt), function(mhap){
              alt = unique(unlist(strsplit(gsub(':\\d+', '',gt[,mhap]), '_')))
              
              alt = paste(alt[!is.na(alt) & alt != '.'], collapse = ',')
            })
            
            
            gt = gsub(':\\d+', '',gt)
            
            # Heterozygous positions
            HetPos = matrix(grepl('_', gt), ncol = ncol(gt), nrow = nrow(gt))
            
            ASVs_attributes_table = NULL
            
            #mhap = which(colnames(gt) == 'SOD2_1')
            
            for(mhap in 1:ncol(gt)){
              temp_gts = gt[,mhap] # genotypes observed in that site
              alleles = strsplit(alt[mhap], ',')[[1]] # alternative alleles observed in that site
              
              if(length(alleles) != 0){
                
                # Vector of presence or absence of each alternative allele
                h_ij = t(sapply(alleles,
                                function(allele){
                                  #allele = alleles[4]
                                  
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
                                      'SNV'}else if(grepl('I', polymorphism)){
                                        'Insertion'
                                      }else if(grepl('D', polymorphism)){
                                        'Deletion'
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
                                  
                                  # CODE NEEDED FOR MASKING SNVs THAT ARE INSIDE OF INDELs
                                  # if('Insertion' %in% type_polymorphism & !is.na(mhaps[mhap,][['homopolymer_regions']])){
                                  #   snv_in_insertion_positions = positions[type_polymorphism == 'Insertion']
                                  #   snv_in_insertion_polymorphisms = polymorphisms[type_polymorphism == 'Insertion']
                                  #   snv_in_insertion_polymorphisms = substr(snv_in_insertion_polymorphisms, 3, 3)
                                  #   
                                  #   if(sum(snv_in_insertion_positions %in% masked_region) > 0){
                                  #     
                                  #     snv_in_insertion_polymorphisms_in_homopolymers = snv_in_insertion_polymorphisms[snv_in_insertion_positions %in% masked_region]
                                  #     snv_in_insertion_positions_in_homopolymers = snv_in_insertion_positions[snv_in_insertion_positions %in% masked_region]
                                  #     
                                  #     start = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '\\d+-'))))
                                  #     end = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '-\\d+'))))
                                  #     
                                  #     homopolymer_repeat = unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '[ATGC]+'))
                                  #     homopolymer_repeat = sapply(homopolymer_repeat, function(repeated_nucleotide){
                                  #       unique(unlist(strsplit(repeated_nucleotide, '')))
                                  #     })
                                  #     
                                  #     homopolymer_repeat = data.frame(
                                  #       start,
                                  #       end,
                                  #       homopolymer_repeat)
                                  #     
                                  #     homopolymer_repeat_motive_in_asv = NULL
                                  #     
                                  #     for(snv_in_insertion_position in snv_in_insertion_positions_in_homopolymers){
                                  #       homopolymer_repeat_motive_in_asv = c(homopolymer_repeat_motive_in_asv,
                                  #                                         homopolymer_repeat[homopolymer_repeat[["start"]] <= snv_in_insertion_position &
                                  #                                         homopolymer_repeat[["end"]] >= snv_in_insertion_position,][["start"]])
                                  #     }
                                  #     
                                  #     homopolymer_repeat_motive_in_asv = homopolymer_repeat[homopolymer_repeat$start %in% homopolymer_repeat_motive_in_asv,]
                                  #     
                                  #     SNV_in_Insertion_in_homopolymer = sum(homopolymer_repeat_motive_in_asv$homopolymer_repeat != snv_in_insertion_polymorphisms_in_homopolymers)
                                  #     
                                  #   }
                                  #   
                                  # }
                                  
                                  
                                  # SNV_in_homopolymer pattern and replacement
                                  
                                  if(SNV_in_homopolymer == 1){
                                    
                                    SNV_in_homopolymer_pattern = allele
                                    
                                    snv_positions_remined = snv_positions[!(snv_positions %in% masked_region)]
                                    snv_polymorphisms_remined = snv_polymorphisms[!(snv_positions %in% masked_region)]
                                    
                                    indel_postions_to_add = positions[type_polymorphism %in% c('Insertion', 'Deletion')]
                                    indel_polymorphisms_to_add = polymorphisms[type_polymorphism %in% c('Insertion', 'Deletion')]
                                    
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
                                  
                                  if(sum(c('Insertion', 'Deletion') %in% type_polymorphism & !is.na(mhaps[mhap,][['homopolymer_regions']])) > 0){
                                    
                                    deleted_positions = positions[type_polymorphism %in% c('Deletion')]
                                    inserted_positions = positions[type_polymorphism %in% c('Insertion')]
                                    
                                    deleted_polymorphisms = polymorphisms[type_polymorphism %in% c('Deletion')]
                                    inserted_polymorphisms = polymorphisms[type_polymorphism %in% c('Insertion')]
                                    
                                    start = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '\\d+-'))))
                                    end = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['homopolymer_regions']], '-\\d+'))))
                                    
                                    masked_region_for_deletions = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos], end[pos], 1)})))
                                    
                                    masked_region_for_insertions = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos] - 1, end[pos], 1)})))
                                    
                                    Insertions_in_homopolymer = sum(inserted_positions %in% masked_region_for_insertions) > 0
                                    Deletions_in_homopolymer = sum(deleted_positions %in% masked_region_for_deletions) > 0
                                    
                                  }else{
                                    Insertions_in_homopolymer = FALSE
                                    Deletions_in_homopolymer = FALSE
                                  }
                                  
                                  
                                  # INDEL_in_homopolymer pattern and replacement
                                  
                                  INDEL_in_homopolymer = Insertions_in_homopolymer == 1 | Deletions_in_homopolymer == 1
                                  
                                  if(INDEL_in_homopolymer == 1){
                                    
                                    INDEL_in_homopolymer_pattern = allele
                                    
                                    inserted_positions_remined = inserted_positions[!(inserted_positions %in% masked_region_for_insertions)]
                                    inserted_polymorphisms_remined = inserted_polymorphisms[!(inserted_positions %in% masked_region_for_insertions)]
                                    
                                    deleted_positions_remined = deleted_positions[!(deleted_positions %in% masked_region_for_deletions)]
                                    deleted_polymorphisms_remined = deleted_polymorphisms[!(deleted_positions %in% masked_region_for_deletions)]
                                    
                                    snv_postions_to_add = positions[type_polymorphism == 'SNV']
                                    snv_polymorphisms_to_add = polymorphisms[type_polymorphism == 'SNV']
                                    
                                    INDEL_in_homopolymer_replacement = data.frame(positions = c(inserted_positions_remined,
                                                                                                deleted_positions_remined,
                                                                                                snv_postions_to_add),
                                                                                  polymorphisms = c(inserted_polymorphisms_remined,
                                                                                                    deleted_polymorphisms_remined,
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
                                  
                                  # Variant_in_STR
                                  
                                  if(!is.na(mhaps[mhap,][['STR_regions']])){
                                    
                                    SNV_DEL_positions = positions[type_polymorphism %in% c('Deletion', 'SNV')]
                                    INS_positions = positions[type_polymorphism %in% c('Insertion')]
                                    
                                    SNV_DEL_polymorphisms = polymorphisms[type_polymorphism %in% c('Deletion', 'SNV')]
                                    INS_polymorphisms = polymorphisms[type_polymorphism %in% c('Insertion')]
                                    
                                    start = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['STR_regions']], '\\d+-'))))
                                    end = as.integer(gsub('-', '', unlist(str_extract_all(mhaps[mhap,][['STR_regions']], '-\\d+'))))
                                    
                                    masked_region_for_SNV_DEL = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos], end[pos], 1)})))
                                    
                                    masked_region_for_INS = c(unlist(sapply(1:length(start), function(pos){
                                      seq(start[pos] - 1, end[pos], 1)})))
                                    
                                    Insertion_in_STR = sum(positions %in% INS_positions & positions %in% masked_region_for_INS) > 0
                                    SNV_or_Deletion_in_homopolymer = sum(positions %in% masked_region_for_SNV_DEL) > 0
                                    
                                  }else{
                                    Insertion_in_STR = FALSE
                                    SNV_or_Deletion_in_homopolymer = FALSE
                                  }
                                  
                                  # Variant_in_STR pattern and replacement
                                  
                                  Variant_in_STR = Insertion_in_STR == 1 | SNV_or_Deletion_in_homopolymer == 1
                                  
                                  if(Variant_in_STR == 1){
                                    
                                    Variant_in_STR_pattern = allele
                                    
                                    INS_positions_remined = INS_positions[!(INS_positions %in% masked_region_for_INS)]
                                    INS_polymorphisms_remined = INS_polymorphisms[!(INS_positions %in% masked_region_for_INS)]
                                    
                                    SNV_DEL_positions_remined = SNV_DEL_positions[!(SNV_DEL_positions %in% masked_region_for_SNV_DEL)]
                                    SNV_DEL_polymorphisms_remined = SNV_DEL_polymorphisms[!(SNV_DEL_positions %in% masked_region_for_SNV_DEL)]
                                    
                                    Variant_in_STR_replacement = data.frame(positions = c(INS_positions_remined,
                                                                                                SNV_DEL_positions_remined),
                                                                                  polymorphisms = c(INS_polymorphisms_remined,
                                                                                                    SNV_DEL_polymorphisms_remined))
                                    
                                    Variant_in_STR_replacement = Variant_in_STR_replacement[order(Variant_in_STR_replacement$positions),]
                                    
                                    Variant_in_STR_replacement = gsub(' ',
                                                                      '',
                                                                      paste(apply(Variant_in_STR_replacement, 
                                                                                  1, 
                                                                                  function(pos){paste(pos, collapse = '')}), 
                                                                            collapse = ''))
                                    
                                    if(Variant_in_STR_replacement == ''){
                                      Variant_in_STR_replacement = '.'
                                    }
                                    
                                  }else{
                                    
                                    Variant_in_STR_pattern = NA
                                    Variant_in_STR_replacement = NA
                                    
                                  }
                                  
                                  
                                  # INDELs in flanking regions
                                  
                                  flanking_Insertion = as.integer(grepl(paste0(paste('(^0', mhaps[mhap,][['length']], sep = '|'),')I(=|\\.)[ATGC]+') ,allele))
                                  flanking_Deletion = as.integer(grepl(paste0(paste('(^1', mhaps[mhap,][['length']], sep = '|'),')D(=|\\.)[ATGC]+') ,allele)) # ERROR: this recognition pattern only takes into account deletions that occurs at the last position and not deletions that includes the last position
                                  
                                  flanking_INDEL = flanking_Insertion == 1 | flanking_Deletion == 1
                                  
                                  #print(flanking_INDEL)
                                  
                                  if(flanking_INDEL == 1){
                                    temp_flanking_Insertion_pattern = str_extract_all(allele, paste0(paste('(^0', paste0('[ATGC]?', mhaps[mhap,][['length']]), sep = '|'),')I(=|\\.)[ATGC]+\\d?'))[[1]]
                                    temp_flanking_Insertion_replacement = gsub(paste0(paste('(^0', mhaps[mhap,][['length']], sep = '|'),')I(=|\\.)[ATGC]+'), '', temp_flanking_Insertion_pattern)
                                    
                                    temp_flanking_Deletion_pattern = str_extract_all(allele, paste0(paste('(^1', paste0('[ATGC]?', mhaps[mhap,][['length']]), sep = '|'),')D(=|\\.)[ATGC]+\\d?'))[[1]]
                                    temp_flanking_Deletion_replacement = gsub(paste0(paste('(^1', mhaps[mhap,][['length']], sep = '|'),')D(=|\\.)[ATGC]+'), '', temp_flanking_Deletion_pattern)
                                    
                                    
                                    temp_flanking_INDEL_pattern = c(temp_flanking_Insertion_pattern,
                                                                    temp_flanking_Deletion_pattern
                                                                    )
                                    
                                    temp_flanking_INDEL_replacement = c(temp_flanking_Insertion_replacement,
                                                                        temp_flanking_Deletion_replacement)
                                    
                                    if('' %in% temp_flanking_INDEL_replacement){
                                      
                                      flanking_INDEL_replacement = '.'
                                      
                                    }else{
                                      
                                      flanking_INDEL_replacement = allele
                                      
                                      for(replacement in 1:length(temp_flanking_INDEL_pattern)){
                                        
                                        flanking_INDEL_replacement = gsub(temp_flanking_INDEL_pattern[replacement], temp_flanking_INDEL_replacement[replacement], flanking_INDEL_replacement)
                                      }
                                      
                                    }
                                    
                                    flanking_INDEL_pattern = allele
                                    
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
                                    
                                    mhaps[mhap,][['homopolymer_anchors']],
                                    mhaps[mhap,][['homopolymer_regions']],
                                    
                                    mhaps[mhap,][['STR_regions']],
                                    
                                    SNV_in_homopolymer,
                                    SNV_in_homopolymer_pattern,
                                    SNV_in_homopolymer_replacement,
                                    
                                    INDEL_in_homopolymer,
                                    INDEL_in_homopolymer_pattern,
                                    INDEL_in_homopolymer_replacement,
                                    
                                    Variant_in_STR,
                                    Variant_in_STR_pattern,
                                    Variant_in_STR_replacement,
                                    
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
                                                      'homopolymer_anchors',
                                                      'homopolymer_regions',
                                                      'STR_regions',
                                                      'SNV_in_homopolymer',
                                                      'SNV_in_homopolymer_pattern',
                                                      'SNV_in_homopolymer_replacement',
                                                      'INDEL_in_homopolymer',
                                                      'INDEL_in_homopolymer_pattern',
                                                      'INDEL_in_homopolymer_replacement',
                                                      'Variant_in_STR',
                                                      'Variant_in_STR_pattern',
                                                      'Variant_in_STR_replacement',
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
                
                ASVs_attributes_table_temp[['Variant_in_STR']] = as.logical(ASVs_attributes_table_temp[['Variant_in_STR']])
                
                #print(ASVs_attributes_table_temp[['flanking_INDEL']])
                
                #ASVs_attributes_table_temp[['flanking_INDEL']] = as.logical(as.integer(ASVs_attributes_table_temp[['flanking_INDEL']]))
                
                ASVs_attributes_table_temp[['bimera']] = as.logical(ASVs_attributes_table_temp[['bimera']])
                
                ASVs_attributes_table_temp[['p_ij']] = ASVs_attributes_table_temp[['P_ij']]/nrow(gt) # CORRECT THIS LINE OR ADD AN ARGUMENT TO CHANGE THE DENOMINATOR
                
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
              
              SNV_in_homopolymer_pattern = case_when(
                SNV_in_homopolymer_pattern == 'NA' ~ NA,
                SNV_in_homopolymer_pattern != 'NA' ~ SNV_in_homopolymer_pattern),
              
              Variant_in_STR_replacement = case_when(
                Variant_in_STR_replacement == 'NA' ~ NA,
                Variant_in_STR_replacement != 'NA' ~ Variant_in_STR_replacement),
              
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
                                                     homopolymer_anchors,
                                                     homopolymer_regions,
                                                     STR_regions,
                                                     SNV_in_homopolymer,
                                                     SNV_in_homopolymer_pattern,
                                                     SNV_in_homopolymer_replacement,
                                                     INDEL_in_homopolymer,
                                                     INDEL_in_homopolymer_pattern,
                                                     INDEL_in_homopolymer_replacement,
                                                     Variant_in_STR,
                                                     Variant_in_STR_pattern,
                                                     Variant_in_STR_replacement,
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
                                                          'PfMDR2'),
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
  print("Filtering drug resistant markers...")
  if('gene_id' %in% colnames(markers)){
    markers_of_interest = markers[grep(paste(gene_ids, collapse = "|"), markers$gene_id),]
  }else{
    markers = cbind(markers, get_gene_description(markers, gff = gff_file))
    markers_of_interest = markers[grep(paste(gene_ids, collapse = "|"), markers$gene_id),]
  }
  
  if(is.null(gene_names)){
    gene_names = gene_ids
  }
  
  
  # Calculates the start and end position of each drugR marker on the CDS of each gene---
  
  ## Start and end position in 3D7 CDSs---
  markers_of_interest[['start_cds']] = NA
  markers_of_interest[['end_cds']] = NA
  markers_of_interest[['strand']] = NA
  markers_of_interest[['ref_length']] = NA
  markers_of_interest[['upstream_no_cds_pos']] = NA
  markers_of_interest[['downstream_no_cds_pos']] = NA
  markers_of_interest[['intronic_pos']] = NA
  markers_of_interest[['length']] = markers_of_interest[['end']] - markers_of_interest[['start']] + 1
  
  print("Analyzing each amplicon in drug resistance markers...")
  for(amplicon in markers_of_interest$amplicon){ # for each drugR marker
    
    # Gene where the drugR marker is located
    gene = markers_of_interest[markers_of_interest$amplicon == amplicon,][['gene_id']]
    
    # Filter the genes attributes of the CDS
    temp_gff = reference_gff[grepl(gene, reference_gff$attributes, ignore.case = TRUE) &
                               reference_gff$type == 'CDS',]
    
    # start counting from 0 taking into account the CDSs (exons) where the drugR marker is located
    start_cds_pos = 0
    
    # check how many CDSs does the gene have
    if(nrow(temp_gff) > 1){ # if the gene has 2 or more CDSs then
      
      # Identify on which CDS(s) the start position of the drugR marker is located
      cds_start = which(temp_gff$start < markers_of_interest[markers_of_interest$amplicon== amplicon,]$start &
                          temp_gff$end > markers_of_interest[markers_of_interest$amplicon== amplicon,]$start)
      
      # Identify on which CDS(s) the end position of the drugR marker is located
      cds_end = which(temp_gff$start < markers_of_interest[markers_of_interest$amplicon== amplicon,]$end &
                        temp_gff$end > markers_of_interest[markers_of_interest$amplicon== amplicon,]$end)
      
      if(length(cds_start) == 0 & length(cds_end) > 0){# if 5' for the amplicon is outside of the exon
        
        for(pre_cds in 1:(cds_end - 1)){ # sum the length of all previous CDSs
          start_cds_pos = start_cds_pos + 
            temp_gff[pre_cds,][['end']] - temp_gff[pre_cds,][['start']] + 1 # length of each previous CDS
        }
        
        # Starting position of the amplicon in the CDS
        
        start_cds_pos = start_cds_pos + 1
        
        # upstream nucleotides located in no CDS regions
        
        upstream_no_cds_pos = abs(markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[cds_end,][['start']])
        
        # end position of the amplicon in the CDS
        
        end_cds = start_cds_pos + markers_of_interest[markers_of_interest$amplicon== amplicon,]$end - temp_gff[cds_end,][['start']]
        
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds_pos
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[cds_end,][['strand']])
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['upstream_no_cds_pos']] = upstream_no_cds_pos
        
      }else if(length(cds_start) > 0 & length(cds_end) == 0){# if 3' for the amplicon is outside of the exon
        
        for(pre_cds in 1:(cds_start - 1)){ # sum the length of all previous CDSs
          start_cds_pos = start_cds_pos + 
            temp_gff[pre_cds,][['end']] - temp_gff[pre_cds,][['start']] + 1 # length of each previous CDS
        }
        
        # Starting position of the amplicon in the CDS
        
        start_cds_pos = start_cds_pos + markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[cds_start,][['start']] + 1
        
        # end position of the amplicon in the CDS
        
        end_cds = start_cds_pos + temp_gff[cds_start,][['end']] - markers_of_interest[markers_of_interest$amplicon== amplicon,]$start 
        
        # upstream nucleotides located in no CDS regions
        
        downstream_no_cds_pos = temp_gff[cds_start,][['end']] - markers_of_interest[markers_of_interest$amplicon== amplicon,]$start + 1 + 1
        
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds_pos
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[cds_start,][['strand']])
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['downstream_no_cds_pos']] = downstream_no_cds_pos
        
      }else if(cds_start != cds_end & length(cds_end) != 0 & length(cds_end) != 0){ # if start and end are in different cds
        
        for(pre_cds in 1:(cds_start - 1)){ # sum the length of all previous CDSs
          start_cds_pos = start_cds_pos + 
            temp_gff[pre_cds,][['end']] - temp_gff[pre_cds,][['start']] + 1 # length of each previous CDS
        }
        
        # Starting position of the amplicon in the CDS
        
        start_cds_pos = start_cds_pos + markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[cds_start,][['start']] + 1
        
        # end position of the amplicon in the CDS
        
        end_cds = 0
        
        for(pre_cds in 1:(cds_end - 1)){ # sum the length of all previous CDSs
          end_cds = end_cds + 
            temp_gff[pre_cds,][['end']] - temp_gff[pre_cds,][['start']] + 1 # length of each previous CDS
        }
        
        end_cds = end_cds + markers_of_interest[markers_of_interest$amplicon== amplicon,]$end - temp_gff[cds_end,][['start']] + 1
        
        # end_cds = start_cds_pos +  
        #   temp_gff[cds_start,][['end']] - markers_of_interest[markers_of_interest$amplicon== amplicon,]$start + # Remaining nucleotides in the exon where start position is located
        #   markers_of_interest[markers_of_interest$amplicon== amplicon,]$end - temp_gff[cds_end,][['start']] + 1 # nucleotides in the exon where end position is located
        
        if(cds_end - cds_start > 1){
          print("cds_end - cds_start > 1")
          
          intronic_pos = NULL
          
          for(internal_cds in (cds_start + 1):(cds_end - 1)){
            end_cds = end_cds + 
              temp_gff[internal_cds,][['end']] - temp_gff[internal_cds,][['start']] + 1
            
            first_intronic_pos_in_amplicon = temp_gff[cds_start,][['end']] - markers_of_interest[markers_of_interest$amplicon == amplicon,]$start + 1 + 1
            last_intronic_pos_in_amplicon = first_intronic_pos_in_amplicon + temp_gff[cds_end,][['start']] - temp_gff[cds_start,][['end']] - 1 - 1
            
            temp_intronic_pos = paste(first_intronic_pos_in_amplicon,
                                 last_intronic_pos_in_amplicon,
                                 sep = ',')
            
            if(internal_cds > cds_start + 1){
              intronic_pos = temp_intronic_pos
            }else{
              intronic_pos = paste(intronic_pos, temp_intronic_pos, sep = ';')
            }
            
          }
          
        }else{
          
          first_intronic_pos_in_amplicon = temp_gff[cds_start,][['end']] - markers_of_interest[markers_of_interest$amplicon == amplicon,]$start + 1 + 1
          last_intronic_pos_in_amplicon = first_intronic_pos_in_amplicon + temp_gff[cds_end,][['start']] - temp_gff[cds_start,][['end']] - 1 - 1
          
          intronic_pos = paste(first_intronic_pos_in_amplicon,
                               last_intronic_pos_in_amplicon,
                               sep = ',')
          
        }
        
        # upstream nucleotides located in no CDS regions
        
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds_pos
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[cds_start,][['strand']])
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['intronic_pos']] = intronic_pos
        
      }else if(cds_start == cds_end & cds_start > 1){ # if the drugR marker is located in CDS 2 or above then
        
        for(pre_cds in 1:(cds_start - 1)){ # sum the length of all previous CDSs
          start_cds_pos = start_cds_pos + temp_gff[pre_cds,][['end']] - temp_gff[pre_cds,][['start']] + 1
        }
        
        # sum the position on the located CDS
        start_cds_pos = start_cds_pos + markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[cds_start,][['start']] + 1
        end_cds = start_cds_pos + markers_of_interest[markers_of_interest$amplicon== amplicon,]$length - 1
        
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds_pos
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[cds_start,][['strand']])
        
        
      }else if(cds_start == cds_end & cds_start == 1){ # if the drugR marker is located in CDS 1
        
        # Take the position on the located CDS
        
        start_cds_pos = start_cds_pos + markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[cds_start,][['start']] + 1
        end_cds = start_cds_pos + markers_of_interest[markers_of_interest$amplicon== amplicon,]$length - 1
        
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds_pos
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
        markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[cds_start,][['strand']])
        
      }
      
    }else{# if the gene has 1 CDS then
      
      # Take the position on the CDS 1
      start_cds_pos = start_cds_pos + markers_of_interest[markers_of_interest$amplicon== amplicon,]$start - temp_gff[['start']] + 1
      end_cds = start_cds_pos + markers_of_interest[markers_of_interest$amplicon== amplicon,]$length - 1
      
      markers_of_interest[markers_of_interest$amplicon== amplicon,][['start_cds']] = start_cds_pos
      markers_of_interest[markers_of_interest$amplicon== amplicon,][['end_cds']] = end_cds
      markers_of_interest[markers_of_interest$amplicon== amplicon,][['strand']] = as.character(temp_gff[['strand']])
      
    }
  }
  
  
  # Generate a reference sequence of each gene---
  # Genes located in the negative strand will not be transformed to their reverse complement yet
  
  # Reference sequences from 3D7---
  ref_seqs = NULL
  
  for(gene in unique(markers_of_interest$gene_id)){
    
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
    
    markers_of_interest[markers_of_interest$gene_id == gene,'ref_length'] = nchar(temp_refseq)
    
  }
  
  print("Generating the reference sequences")
  names(ref_seqs) = unique(markers_of_interest$gene_id)
  ref_seqs = DNAStringSet(ref_seqs)
  
  
  # Convert cigar format to standard mutation nomenclature PMC1867422---
  
  ## Filter drugR markers---
  # markers of interest loci abundance table
  print("Filter drugR markers")
  if(length(gene_ids) > 1){
    moi_loci_abd_table = ampseq_object@gt[, markers$gene_id %in% gene_ids]  
  }else{
    moi_loci_abd_table = matrix(ampseq_object@gt[, markers$gene_id %in% gene_ids],
                                ncol = 1,
                                dimnames = list(rownames(ampseq_object@gt),
                                                colnames(ampseq_object@gt)[markers$gene_id %in% gene_ids]))
  }
  
  
  ## Remove read abundace---
  print("Removing read abundance...")
  moi_loci_abd_table = gsub(":[0-9]+", "", moi_loci_abd_table)
  
  moi_loci_dna_table = moi_loci_abd_table
  moi_loci_aa_table = moi_loci_abd_table
  
  for(amplicon in colnames(moi_loci_abd_table)){ # For each amplicon in columns
    for(sample in (1:nrow(moi_loci_abd_table))){ # For each sample in rows
      
      locus = moi_loci_abd_table[sample, amplicon] # Get the genotype in the locus
      
      
      locus = gsub('\\d+(D|I)=[ATGC]+', '', locus) # REMOVE WHEN INDELs DETECTION IS IMPLEMENTED

      if(!is.na(locus)){# REMOVE WHEN INDELs DETECTION IS IMPLEMENTED

        if(locus == ''){
          locus = '.'
        }

        if(grepl('^_', locus)){
          locus = gsub('^_', '._', locus)
        }

        if(grepl('\\._\\.', locus)){
          locus = '.'
        }

        if(grepl('_', locus)){
          locus = gsub('^_', '._', locus)
        }

      }

      
      if(is.na(locus)){ # if the locus is NULL complete the cell with NA
        
        moi_loci_dna_table[sample, amplicon] = NA
        moi_loci_aa_table[sample, amplicon] = NA
        
      }else{
        
        clones = unlist(strsplit(locus, '_')) # get all different clones at that locus
        
        dna_clones = NULL
        aa_clones = NULL
        
        for(clone in clones){ # for each clone get alleles of all different SNPs
          
          alleles = paste0(unlist(str_extract_all(clone, '\\d+')),
                           unlist(strsplit(gsub('^[0-9]+', '', clone), '[0-9]+')))
          
          
          if(!is.na(markers_of_interest[markers_of_interest[['amplicon']] == amplicon, ][['downstream_no_cds_pos']])){
            downstream_no_cds_pos = markers_of_interest[markers_of_interest[['amplicon']] == amplicon, ][['downstream_no_cds_pos']]
            final_pos = markers_of_interest[markers_of_interest[['amplicon']] == amplicon, ][['length']]
            
            intron_alleles = alleles[as.integer(str_extract(alleles, '\\d+')) %in% downstream_no_cds_pos:final_pos]
            
            alleles = alleles[!(as.integer(str_extract(alleles, '\\d+')) %in% downstream_no_cds_pos:final_pos)]
            
            if(length(alleles) == 0){
              alleles = '.'
            }
            
          }
          
          if(!is.na(markers_of_interest[markers_of_interest[['amplicon']] == amplicon, ][['upstream_no_cds_pos']])){
            upstream_no_cds_pos = markers_of_interest[markers_of_interest[['amplicon']] == amplicon, ][['upstream_no_cds_pos']]
            final_pos = markers_of_interest[markers_of_interest[['amplicon']] == amplicon, ][['length']]
            
            intron_alleles = alleles[as.integer(str_extract(alleles, '\\d+')) %in% 1:upstream_no_cds_pos]
            
            alleles = alleles[!(as.integer(str_extract(alleles, '\\d+')) %in% 1:upstream_no_cds_pos)]
            
            if(length(alleles) == 0){
              alleles = '.'
            }else{
              
              temp_alleles = data.frame(position = str_extract(alleles, '\\d+'),
                                        allele = gsub('\\d+', '', alleles))
              
              temp_alleles[["position"]] = as.integer(temp_alleles[["position"]]) - 
                as.integer(upstream_no_cds_pos)
              
              alleles = paste0(temp_alleles[["position"]], temp_alleles[["allele"]])
              
            }
          }
          
          
          if(!is.na(markers_of_interest[markers_of_interest[['amplicon']] == amplicon, ][['intronic_pos']])){
            start_intronic_pos = gsub(',\\d+$', '', markers_of_interest[markers_of_interest[['amplicon']] == amplicon, ][['intronic_pos']])
            end_intronic_pos = gsub('^\\d+,', '', markers_of_interest[markers_of_interest[['amplicon']] == amplicon, ][['intronic_pos']])
            
            intron_alleles = alleles[as.integer(str_extract(alleles, '\\d+')) %in% start_intronic_pos:end_intronic_pos]
            
            alleles = alleles[!(as.integer(str_extract(alleles, '\\d+')) %in% start_intronic_pos:end_intronic_pos)]
            
            if(length(alleles) == 0){
              alleles = '.'
            }else{
              
              temp_alleles = data.frame(position = str_extract(alleles, '\\d+'),
                                        allele = gsub('\\d+', '', alleles))
              
              temp_alleles[["position"]] = as.integer(temp_alleles[["position"]]) - 
                (as.integer(end_intronic_pos) - as.integer(start_intronic_pos) + 1)
              
              alleles = paste0(temp_alleles[["position"]], temp_alleles[["allele"]])
              
            }
          }
          
          dna_alleles = NULL
          aa_alleles = NULL
          
          
          if(clone == "."){ # if allele equals to wild type
            
            codons = "."
            
          }else{
            
            if("." %in% alleles){
              
              codons = "."
              
            }else{
              
              codons = data.frame(alleles = gsub('[0-9]','',alleles), cds_position = sapply(alleles, function(allele){
                
                # position in the mhap
                mhap_position = as.integer(str_extract(allele, '\\d+'))
                
                # variant or nucleotide found
                #mhap_variant = gsub('[0-9]+', '', allele)
                
                # calculate position in the CDS
                
                cds_position = markers_of_interest[markers_of_interest$amplicon == amplicon,'start_cds'] + mhap_position - 1
                
              }))
              
              
              # Identify nucleotide in the reference strain
              codons$ref_variant = sapply(1:nrow(codons), function(x){
                
                if(grepl('D=',codons[x, 'alleles'])){
                  
                  gsub('D=', '',codons[x, 'alleles'])
                  
                }else if(grepl('I=',codons[x, 'alleles'])){
                  
                  ''
                  
                }else{
                  
                  as.character(subseq(ref_seqs[which(names(ref_seqs) ==
                                                       markers_of_interest[markers_of_interest$amplicon == amplicon, 'gene_id'])],
                                      start = codons[x, 'cds_position'],
                                      end = codons[x,'cds_position']))
                  
                }
                
              })
              
              # calculate the aminoacid position
              codons$aa_position = ceiling(codons$cds_position/3)
              codons$first_nucleotide = 3*codons$aa_position - 2
              codons$last_nucleotide = 3*codons$aa_position
              
            }
            
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
                                      as.character(Biostrings::translate(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_id'])],
                                                                                start = unique(codons[codons$aa_position == codon,][['first_nucleotide']]),
                                                                                end = unique(codons[codons$aa_position == codon,][['last_nucleotide']])))),
                                      as.character(Biostrings::translate(reverseComplement(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_id'])],
                                                                                                  start = unique(codons[codons$aa_position == codon,][['first_nucleotide']]),
                                                                                                  end = unique(codons[codons$aa_position == codon,][['last_nucleotide']]))))))
              
              
              codon_positons = unique(codons[codons$aa_position == codon,][['first_nucleotide']]):unique(codons[codons$aa_position == codon,][['last_nucleotide']])
              
              positions_from_ref = codon_positons[!(unique(codons[codons$aa_position == codon,][['first_nucleotide']]):unique(codons[codons$aa_position == codon,][['last_nucleotide']]) %in% 
                                                      codons[codons$aa_position == codon,][['cds_position']])]
              
              positions_from_ref_df = NULL
              
              for(positions_from_ref1 in  positions_from_ref){
                positions_from_ref_df = rbind(positions_from_ref_df,
                                              data.frame(alleles = as.character(subseq(ref_seqs[which(names(ref_seqs)==markers_of_interest[markers_of_interest$amplicon == amplicon,'gene_id'])],
                                                                                       start = positions_from_ref1,
                                                                                       end = positions_from_ref1)),
                                                         cds_position = positions_from_ref1))
              }
              
              
              sample_codon = rbind(codons[codons$aa_position == codon,c('alleles', 'cds_position')],
                                   positions_from_ref_df
              )
              
              sample_codon = sample_codon[order(sample_codon$cds_position),]
              
              sample_codon$alleles = gsub('(D|I)=', '', sample_codon$alleles)
              
              mhap_aa_variant = ifelse(markers_of_interest[markers_of_interest$amplicon == amplicon,'strand'] == "+",
                                       as.character(Biostrings::translate(DNAString(paste0(sample_codon$alleles, collapse = "")))),
                                       as.character(Biostrings::translate(reverseComplement(DNAString(paste0(sample_codon$alleles, collapse = ""))))))
              
              
              if(markers_of_interest[markers_of_interest$amplicon == amplicon,'strand'] == "+"){
                dna_alleles = c(dna_alleles, paste(paste0('c.',codons[codons$aa_position == codon,][['cds_position']],
                                                          codons[codons$aa_position == codon,][['ref_variant']],
                                                          '>',
                                                          codons[codons$aa_position == codon,][['alleles']]), collapse = " "))
                
                if(ref_aa_variant != mhap_aa_variant){
                  aa_alleles = c(aa_alleles, paste0(ref_aa_variant, codon, mhap_aa_variant))
                }
                
              }else{
                dna_alleles = c(dna_alleles, paste(paste0('c.',
                                                          markers_of_interest[markers_of_interest$amplicon == amplicon,'ref_length'] - 
                                                            codons[codons$aa_position == codon,][['cds_position']] + 1,
                                                          reverseComplement(DNAStringSet(codons[codons$aa_position == codon,][['ref_variant']])), '>',
                                                          reverseComplement(DNAStringSet(codons[codons$aa_position == codon,][['alleles']]))), collapse = ' '))
                
                if(ref_aa_variant != mhap_aa_variant){
                  aa_alleles = c(aa_alleles, paste0(ref_aa_variant,
                                                    ceiling((markers_of_interest[markers_of_interest$amplicon == amplicon,'ref_length'] - 3*codon + 1)/3),
                                                    mhap_aa_variant))
                }
                
                
              }
              
              
            }
            
            
          }
          
          dna_clones[[clone]] = paste(dna_alleles, collapse = " ")
          aa_clones[[clone]] = paste(aa_alleles, collapse = " ")
          
        }
        
        moi_loci_dna_table[sample, amplicon] = paste(unique(dna_clones), collapse = " / ")
        moi_loci_aa_table[sample, amplicon] = paste(unique(aa_clones), collapse = " / ")
        
      }
    }
  }
  
  
  moi_loci_aa_table[moi_loci_aa_table == '' |
                      moi_loci_aa_table == ' / p.(=)' |
                      moi_loci_aa_table == 'p.(=) / '] = 'p.(=)'
  
  if(plot_haplo_freq){
    
    
    # Empty table to fill cigar outputs
    aacigar_table = matrix(NA,
                           nrow = nrow(moi_loci_aa_table),
                           ncol = length(gene_ids),
                           dimnames = list(rownames(moi_loci_aa_table),
                                           gene_ids))
    
    
    for(gene in 1:length(gene_ids)){
      
      gene_aa = matrix(moi_loci_aa_table[, colnames(moi_loci_aa_table) %in% 
                                           markers_of_interest[markers_of_interest$gene_id == gene_ids[gene],][['amplicon']]],
                       ncol = sum(colnames(moi_loci_aa_table) %in% 
                                    markers_of_interest[markers_of_interest$gene_id == gene_ids[gene],][['amplicon']]),
                       dimnames = list(rownames(moi_loci_aa_table),
                                       colnames(moi_loci_aa_table)[colnames(moi_loci_aa_table) %in% 
                                                                     markers_of_interest[markers_of_interest$gene_id == gene_ids[gene],][['amplicon']]]))
      
      
      # filter amplicons for the gene of interest
      gene_of_interest_info = markers_of_interest[markers_of_interest[['gene_id']] == gene_ids[gene],]
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
    
    if(length(gene_ids) > 1){
      mon_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) == 0),]
      poly_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) != 0),]  
    }else{
      mon_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) == 0),]
      mon_aacigar_table = matrix(mon_aacigar_table, ncol = 1,
                                 dimnames = list(
                                   names(mon_aacigar_table),
                                   gene_ids
                                 ))
      poly_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) != 0),]
      
      if(!isEmpty(poly_aacigar_table)){
        poly_aacigar_table = matrix(poly_aacigar_table, ncol = 1,
                                    dimnames = list(
                                      names(mon_aacigar_table),
                                      gene_ids
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
    
    extended_aacigar_table %<>% pivot_longer(cols = all_of(gene_ids),
                                             names_to = 'gene_ids',
                                             values_to = 'haplotype')
    
    extended_aacigar_table$gene_names = NA
    
    for(gene in 1:length(gene_ids)){
      extended_aacigar_table[extended_aacigar_table$gene_ids == gene_ids[gene],][['gene_names']] = gene_names[gene]
      
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
  
  if(is.null(gene_names)){
    gene_names = gene_ids
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
    
    for(amplicon in drug_markers[drug_markers$gene_id == gene, 'amplicon']){ # for each amplicon in the gene
      
      if(amplicon %in% colnames(loci_aa_table)){ # Apply only for amplicons present in loci table (some amplicons might be removed in previous steps)
        
        # Get all polymorphic positions in the population
        found_positions = loci_aa_table[,amplicon]
        
        found_positions = as.character(unlist(strsplit(as.character(unlist(strsplit(found_positions, " / "))), ' ')))
        
        
        found_positions = unique(as.character(gsub('^[A-Z]',
                                                   '',
                                                   gsub('([A-Z]$|\\*$)',
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
          
          
          if(length(positions) > 0){
            
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
                                          gsub('([A-Z]|\\*)$',
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
            

            
            # Update for each amplicon in the corresponding gene
            if(is.na(aacigar_table[sample, gene])){
              
              aacigar_haplotype = paste(aacigar_haplotype[order(aacigar_haplotype$position),'aacigar_haplotype'], collapse = " ")
              
              aacigar_table[sample, gene] = aacigar_haplotype
              phenotype_table[sample, gene] = sample_phenotype
              
            }else{
              
              previous_aa_cigar = unlist(strsplit(aacigar_table[sample, gene], ' '))
              
              previous_aa_cigar_table = data.frame(position = as.integer(str_extract(previous_aa_cigar, '\\d+')),
                                                   aacigar_haplotype = previous_aa_cigar
                                                   )
              
              aacigar_haplotype = rbind(previous_aa_cigar_table,
                                              aacigar_haplotype)
              
              aacigar_haplotype = paste(aacigar_haplotype[order(aacigar_haplotype$position),'aacigar_haplotype'], collapse = " ")
              
              aacigar_table[sample, gene] = aacigar_haplotype
              
              #aacigar_table[sample, gene] = paste(aacigar_table[sample, gene], aacigar_haplotype, sep = ' ')
              phenotype_table[sample, gene] = paste(phenotype_table[sample, gene], sample_phenotype, sep = '; ')
              
            }
            
          }else{
            
            sample_clones = loci_aa_table[sample,
                                          amplicon]
            
            if(!is.na(sample_clones) & sample_clones == 'p.(=)'){
              
              # Update for each amplicon in the corresponding gene
              if(is.na(aacigar_table[sample, gene])){
                
                phenotype_table[sample, gene] = 'Sensitive phenotype'
                
              }else{
                
                phenotype_table[sample, gene] = paste(phenotype_table[sample, gene], 'Sensitive phenotype', sep = '; ')
                
              }
              
            }else if(is.na(sample_clones)){
              
              # Update for each amplicon in the corresponding gene
              if(is.na(aacigar_table[sample, gene])){
                
                phenotype_table[sample, gene] = paste0(amplicon, ' amplicon did not amplify')
                
              }else{
                
                phenotype_table[sample, gene] = paste(phenotype_table[sample, gene], paste0(amplicon, ' amplicon did not amplify'), sep = '; ')
                
              }
              
            }
            
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
    summarise(allele_count = n(), .by = c(gene_names, gene_ids, var1, var2, haplotype))
  
  haplotype_counts2 = extended_aacigar_table %>%
    summarise(count = 1, .by = c(Sample_id, gene_names, gene_ids, var1, var2, haplotype))
  
  haplotype_counts2 = haplotype_counts2 %>%
    summarise(sample_count = n(), .by = c(gene_names, gene_ids, var1, var2, haplotype))
  
  haplotype_counts = left_join(haplotype_counts,
                               haplotype_counts2,
                               by = c('gene_names', 'gene_ids', 'var1', 'var2', 'haplotype'))
  
  genotype_phenotype_table_w_meta = left_join(genotype_phenotype_table,
                                       metadata,
                                       by = join_by('Sample' == 'Sample_id')
                                       )
  
  print('Calculating haplotype frequencies')
  
  haplotype_counts$genotyped_samples = NA
  haplotype_counts$total_nalleles = NA
  
  haplotype_counts$prev = NA
  haplotype_counts$prev_lower = NA
  haplotype_counts$prev_upper = NA
  haplotype_counts$freq = NA
  haplotype_counts$freq_lower = NA
  haplotype_counts$freq_upper = NA
  
  
  for(gene in levels(as.factor(haplotype_counts$gene_ids))){
    for(Pop in levels(as.factor(haplotype_counts[haplotype_counts$gene_ids == gene,][['var1']]))){
      for(date in levels(as.factor(haplotype_counts[haplotype_counts$gene_ids == gene&
                                                    haplotype_counts$var1 == Pop,][['var2']]))){
        
        total_number_of_alleles = sum(haplotype_counts[haplotype_counts$gene_ids == gene&
                                                         haplotype_counts$var1 == Pop&
                                                         haplotype_counts$var2 == date,][['allele_count']])
        allele_count = haplotype_counts[haplotype_counts$gene_ids == gene&
                                          haplotype_counts$var1 == Pop&
                                          haplotype_counts$var2 == date,][['allele_count']]
        
        total_number_of_genotyped_samples = nrow(genotype_phenotype_table_w_meta[genotype_phenotype_table_w_meta[['var1']] == Pop &
                                                                                  genotype_phenotype_table_w_meta[['var2']] == date &
                                                                                  genotype_phenotype_table_w_meta[['Gene']] == gene &
                                                                                  !grepl('Gene .+ did not amplify', genotype_phenotype_table_w_meta[['Phenotype']])
                                                                                ,])
        number_of_samples_with_allele = haplotype_counts2[haplotype_counts2$gene_ids == gene&
                                                            haplotype_counts2$var1 == Pop&
                                                            haplotype_counts2$var2 == date,][['sample_count']]
        
        temp_freq = binconf(allele_count,
                            total_number_of_alleles,
                            method = 'exact'
        )
        
        temp_prev = binconf(number_of_samples_with_allele,
                            total_number_of_genotyped_samples,
                            method = 'exact'
        )
        
        haplotype_counts[haplotype_counts$gene_ids == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['genotyped_samples']] = total_number_of_genotyped_samples
        
        haplotype_counts[haplotype_counts$gene_ids == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['total_nalleles']] = total_number_of_alleles
        
        haplotype_counts[haplotype_counts$gene_ids == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['freq']] = temp_freq[,1]
        
        haplotype_counts[haplotype_counts$gene_ids == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['freq_lower']] = temp_freq[,2]
        
        haplotype_counts[haplotype_counts$gene_ids == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['freq_upper']] = temp_freq[,3]
        
        haplotype_counts[haplotype_counts$gene_ids == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['prev']] = temp_prev[,1]
        
        haplotype_counts[haplotype_counts$gene_ids == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['prev_lower']] = temp_prev[,2]
        
        haplotype_counts[haplotype_counts$gene_ids == gene&
                           haplotype_counts$var1 == Pop&
                           haplotype_counts$var2 == date,][['prev_upper']] = temp_prev[,3]
        
      }
    }
  }
  
  haplotype_counts$compact_haplotype = NA
  
  unique_haplotypes = unique(haplotype_counts$haplotype)
  
  unique_haplotypes = unique_haplotypes[!is.na(unique_haplotypes)]
  
  for(haplo in unique_haplotypes){
    
    compact_haplo = unlist(strsplit(haplo, ' '))
    
    ref = gsub('\\d+\\w+$', '', compact_haplo)
    alt = gsub('^\\w+\\d+', '', compact_haplo)
    
    compact_haplo = compact_haplo[ref != alt]
    
    if(length(compact_haplo) == 0){
      compact_haplo = 'p.(=)'
    }
    
    compact_haplo = paste(compact_haplo, collapse = ' ')
    
    haplotype_counts[haplotype_counts[['haplotype']] == haplo &
                       !is.na(haplotype_counts[['haplotype']]), ][['compact_haplotype']] = compact_haplo
    
  }
  
  haplotype_counts %<>% mutate(gene_haplo = paste(gene_names, haplotype, sep = ": "),
                               gene_haplo_compact = paste(gene_names, compact_haplotype, sep = ": "))
  
  haplotype_counts$phenotype = NA
  
  for(haplotype in haplotype_counts$haplotype){
    
    if(!is.na(haplotype)){
      haplotype_counts[haplotype_counts$haplotype == haplotype & !is.na(haplotype_counts$haplotype), ][['phenotype']] =
        genotype_phenotype_match[genotype_phenotype_match$Genotype == haplotype & !is.na(genotype_phenotype_match$Genotype), ][['Phenotype']]
    }
    
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
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| |,)+ resistance')), '\\d+')))) == 1){
        reds[2]
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| |,)+ resistance')), '\\d+')))) == 2){
        reds[3]
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| |,)+ resistance')), '\\d+')))) == 3){
        reds[5]
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| |,)+ resistance')), '\\d+')))) >= 4){
        reds[6]
      }else if(grepl('Sensitive', Phenotype)){
        blue
      }else if(grepl('polymorphism', Phenotype) &
               !grepl('resistance', Phenotype)){
        'gold3'
      }else if(grepl('variant unreported', Phenotype) & !grepl('resistance', Phenotype)){
        orange
      }else if(grepl("^Amplicon.+amplify(; .+)?$", Phenotype)){
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
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| |,)+ resistance')), '\\d+')))) == 1){
        4
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| |,)+ resistance')), '\\d+')))) == 2){
        3
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| |,)+ resistance')), '\\d+')))) == 3){
        2
      }else if(grepl('resistance', Phenotype) & 
               sum(as.integer(unlist(str_extract_all(unlist(str_extract_all(Phenotype, '\\d+ ([A-z]| |,)+ resistance')), '\\d+')))) >= 4){
        1
      }else if(grepl('Sensitive', Phenotype)){
        8
      }else if(grepl('polymorphism', Phenotype) &
               !grepl('resistance', Phenotype)){
        7
      }else if(grepl('variant unreported', Phenotype) & !grepl('resistance', Phenotype)){
        6
      }else if(grepl("^Amplicon.+amplify(; .+)?$", Phenotype)){
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
  
  genotype_phenotype_match_sorted$compact_haplotype = NA
  
  unique_genotpyes = unique(genotype_phenotype_match_sorted$Genotype)
  
  unique_genotpyes = unique_genotpyes[!is.na(unique_genotpyes)]
  
  for(haplo in unique_genotpyes){
    
    compact_haplo = unlist(strsplit(haplo, ' '))
    
    ref = gsub('\\d+\\w+$', '', compact_haplo)
    alt = gsub('^\\w+\\d+', '', compact_haplo)
    
    compact_haplo = compact_haplo[ref != alt]
    
    if(length(compact_haplo) == 0){
      compact_haplo = 'p.(=)'
    }
    
    compact_haplo = paste(compact_haplo, collapse = ' ')
    
    genotype_phenotype_match_sorted[genotype_phenotype_match_sorted[['Genotype']] == haplo &
                                      !is.na(genotype_phenotype_match_sorted[['Genotype']]), ][['compact_haplotype']] = compact_haplo
    
  }
  
  
  genotype_phenotype_match_sorted %<>% mutate(
    gene_haplo = paste(Gene_name, Genotype, sep = ': '),
    gene_haplo_compact = paste(Gene_name, compact_haplotype, sep = ': ')
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
    ggplot(aes(y = prev, x = var2, fill  = gene_haplo, alpha = gene_haplo)) +
    geom_col(color = 'gray75')+
    facet_grid(var1 ~ gene_names)+
    scale_fill_manual(values = genotype_phenotype_match_sorted$color_pal, 
                      labels = genotype_phenotype_match_sorted$gene_haplo_compact)+
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
    ggplot(aes(y = prev, x = var2, group  = gene_haplo, color = gene_haplo)) +
    geom_point()+
    geom_errorbar(aes(ymin = prev_lower, ymax = prev_upper), alpha = .5, width = .2)+
    geom_line()+
    facet_grid(var1 ~ gene_names)+
    scale_color_manual(values = genotype_phenotype_match_sorted$color_pal, 
                       labels = genotype_phenotype_match_sorted$gene_haplo_compact)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = 'bottom') +
    labs(y = 'Frequency in population',
         x = 'Date of Collection',
         color = 'Gene: Haplotype')+
    guides(alpha = "none",
           color=guide_legend(ncol=3))
  
  #names(haplotype_counts) = c(names(haplotype_counts)[1], variables[2:3], names(haplotype_counts)[-1:-3])
  
  
  print('generating table for minimal report')
  
  
  
  
  loci_aa_table_long = as.data.frame(loci_aa_table) %>%
    mutate(Sample_id = rownames(loci_aa_table)) %>%
    pivot_longer(cols = all_of(colnames(loci_aa_table)),
                 names_to = 'Marker',
                 values_to = 'Polymorphism'
                )
  
  loci_aa_table_long = left_join(loci_aa_table_long,
                                 drug_markers[,c('amplicon', 'gene_id')],
                                 by = join_by('Marker' == 'amplicon')
                                 )
  
  loci_aa_table_long[['gene_name']] = NA
  
  for(gene in 1:length(gene_ids)){
    loci_aa_table_long[loci_aa_table_long[['gene_id']] == gene_ids[gene], ][['gene_name']] = gene_names[gene]
  }
  
  
  
  minimal_aa_table = NULL

  
  for(gene in unique(loci_aa_table_long[['gene_id']])){
    
    polymorphisms = unique(unlist(strsplit(loci_aa_table_long[loci_aa_table_long[['gene_id']] == gene,][['Polymorphism']], '( / | )')))
    
    polymorphisms = polymorphisms[!is.na(polymorphisms)]
    
    for(polymorphism in polymorphisms){
      
      temp_summary_aa_table = data.frame(Gene_id = gene,
                                         Gene_name = gene_names[which(gene_ids == gene)],
                                         Polymorphism = polymorphism
                                         )
      
      Annotation = drugR_reference_alleles[drugR_reference_alleles$Gene_Id == gene &
                                             drugR_reference_alleles$Mutation == polymorphism, ][['Annotation']]
      
      if(length(Annotation) == 1){
        temp_summary_aa_table[['Annotation']] = Annotation
      }else if(polymorphism != 'p.(=)'){
        temp_summary_aa_table[['Annotation']] = 'Polymorphism respect to the reference strain'
      }else{
        temp_summary_aa_table[['Annotation']] = 'Sensitive phenotype'
        
      }
      
      minimal_aa_table = rbind(minimal_aa_table, temp_summary_aa_table)
      
    }
  }
  
  
  
  minimal_aa_table[['Samples_with_polymorphism']] = NA
  minimal_aa_table[['Genotyped_Samples']] = NA

  
  for(gene in minimal_aa_table[['Gene_id']]){
    polymorphisms = unique(minimal_aa_table[minimal_aa_table[['Gene_id']] == gene,][['Polymorphism']])
    
    minimal_aa_table[minimal_aa_table[['Gene_id']] == gene, ][['Genotyped_Samples']] =
      length(unique(loci_aa_table_long[loci_aa_table_long$gene_id == gene &
                           !is.na(loci_aa_table_long$Polymorphism)
                           , ][['Sample_id']]))
    
    
    for(polymorphism in polymorphisms){
      
      if(polymorphism == 'p.(=)'){
        polymorphism = 'p.\\(=\\)'
      }
      
      minimal_aa_table[minimal_aa_table[['Gene_id']] == gene &
                         grepl(polymorphism, minimal_aa_table[['Polymorphism']])
                         , ][['Samples_with_polymorphism']] = 
        length(unique(loci_aa_table_long[loci_aa_table_long$gene_id == gene &
             !is.na(loci_aa_table_long$Polymorphism) & 
               grepl(polymorphism, loci_aa_table_long$Polymorphism)
             , ][['Sample_id']]))
      
    }
  }
  
  
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
                                   haplotype_freq_barplot = haplotype_freq_barplot,
                                   drug_markers_table = drug_markers,
                                   minimal_aa_table = minimal_aa_table)
    
  }else{
    
    drug_resistant_hap_list = list(aa_mutations = haps_respect_to_ref$loci_aa_table,
                                   dna_mutations = haps_respect_to_ref$loci_dna_table,
                                   genotype_phenotype_table = genotype_phenotype_table,
                                   drug_phenotype_table = drug_phenotype_table,
                                   drug_phenotyope_lineplot = drug_phenotyope_lineplot,
                                   drug_phenotype_barplot = drug_phenotype_barplot,
                                   haplotypes_freq_lineplot = haplotypes_freq_lineplot,
                                   haplotype_freq_barplot = haplotype_freq_barplot,
                                   drug_markers_table = drug_markers,
                                   minimal_aa_table = minimal_aa_table)
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
      ampseq_object = filter_samples(ampseq_object, v = !(is.na(ampseq_object@metadata[[strata]]) | grepl('NA',ampseq_object@metadata[[strata]])), update_cigars = F)
    }else if(length(ampseq_object@metadata[is.na(ampseq_object@metadata[[strata]]) | grepl('NA',ampseq_object@metadata[[strata]]),][[strata]])>0){
      ampseq_object@metadata[is.na(ampseq_object@metadata[[strata]]) | grepl('NA',ampseq_object@metadata[[strata]]),][[strata]] = 'missing data'
    }
    
    if(!is.null(filters)){
      ampseq_object = filter_samples(ampseq_object, v = grepl(filters,ampseq_object@metadata[[strata]]), update_cigars = F)
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
  
  if(length(polyclonals) > 0 & length(polyclonals) < nrow(polygenomic)){
    polygenomic$Clonality = NA
    polygenomic[polyclonals,][['Clonality']] = 'Polyclonal'
    polygenomic[-1*polyclonals,][['Clonality']] = 'Monoclonal'
  }else if(length(polyclonals) == 0){
    polygenomic[['Clonality']] = 'Monoclonal'
  }else if(length(polyclonals) == nrow(polygenomic)){
    polygenomic[['Clonality']] = 'Polyclonal'
  }
  
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
  
  if(!is.null(loci_performance)){
    loci_performance[['prop_poly_detected']] = prop_poly_detected
    
  }else{
    
    loci_performance = data.frame(amplicon = ampseq_object@markers$amplicon,
                                  prop_poly_detected = prop_poly_detected
                                  )
    
  }
  
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

draw_haplotypes = function(ampseq_object = NULL, comparison_list = NULL, max_coi = NULL, combine_lists = F){
  
  if(!combine_lists){
    
    coded_haplotypes = NULL
    
    for(comparison in 1:length(comparison_list)){
      
      temp_ampseq = filter_samples(ampseq_object, comparison_list[[comparison]], update_cigars = FALSE)
      
      gt = temp_ampseq@gt
      gt = gsub(':\\d+', '', gt)
      
      for(mhap in 1:ncol(gt)){
        
        alleles = gt[,mhap]
        alleles = alleles[!is.na(alleles)]
        alleles = unlist(strsplit(alleles, '_')) 
        alleles = unique(alleles)
        
        #if( "4G60C1" %in% alleles){stop()}
        
        if(length(alleles) > 0){
          alleles_labels = 1:length(alleles)
          
          for(allele in alleles_labels){
            
            if(alleles[allele] != '.'){
              
              for(sample in 1:nrow(gt)){
                
                pattern_to_replace = paste0('(^|_)',alleles[allele], '(_|$)')
                
                if(grepl(pattern_to_replace, gt[sample, mhap])){
                  
                  pattern_to_replace = str_extract(gt[sample, mhap], pattern_to_replace)
                  replacement_pattern = gsub(alleles[allele], allele, pattern_to_replace)
                  
                  if(!grepl('^_', pattern_to_replace)){
                    pattern_to_replace = paste0('^', pattern_to_replace)
                  }
                  
                  if(!grepl('_$', pattern_to_replace)){
                    pattern_to_replace = paste0(pattern_to_replace, '$')
                  }
                  
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
          
          if(is.null(max_coi)){
            sample_clonality = clonality[poly_sample]
          }else{
            sample_clonality = max_coi
          }
          
          temp_gt_polyclonal = matrix(NA, 
                                      nrow = sample_clonality,
                                      ncol = ncol(gt),
                                      dimnames = list(
                                        paste0(comparison_list[[comparison]][poly_sample], '_C',1:sample_clonality),
                                        colnames(gt)
                                      ))
          
          alleles = strsplit(gt[poly_sample,], '_')
          
          for(clone in 1:sample_clonality){
            
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
      
      if(!is.null(names(comparison_list))){
        coded_haplotypes_temp$Comparison = names(comparison_list)[comparison]
      }else{
        coded_haplotypes_temp$Comparison = paste0('Comparison ', comparison)
      }

      coded_haplotypes = rbind(coded_haplotypes, coded_haplotypes_temp)
      
    }
    
    haplotypes_plot = coded_haplotypes %>%
      filter(!is.na(Allele))%>%
      ggplot(aes(x = Marker, y = Haplotype, fill = Allele)) +
      geom_tile()+
      scale_fill_brewer(palette = 'Set1')+
      facet_grid(Comparison~chromosome, scales = 'free', space = 'free', switch = 'y')+
      theme(axis.text.x = element_blank())
    
  }else{
    
    if(is.null(names(comparison_list))){
      names(comparison_list) = paste0('Comparison ', length(comparison_list))
    }
    
    comparison_df = data.frame(Strata = c(sapply(names(comparison_list), function(Strata){
      rep(Strata, length(comparison_list[[Strata]]))})),
               Sample =unlist(comparison_list))
    
    temp_ampseq = filter_samples(ampseq_object, comparison_df[['Sample']], update_cigars = FALSE)
    
    gt = temp_ampseq@gt
    gt = gsub(':\\d+', '', gt)
    
    for(mhap in 1:ncol(gt)){
      
      alleles = gt[,mhap]
      alleles = alleles[!is.na(alleles)]
      alleles = unlist(strsplit(alleles, '_')) 
      alleles = unique(alleles)
      
      #if( "4G60C1" %in% alleles){stop()}
      
      if(length(alleles) > 0){
        alleles_labels = 1:length(alleles)
        
        for(allele in alleles_labels){
          
          if(alleles[allele] != '.'){
            
            for(sample in 1:nrow(gt)){
              
              pattern_to_replace = paste0('(^|_)',alleles[allele], '(_|$)')
              
              if(grepl(pattern_to_replace, gt[sample, mhap])){
                
                pattern_to_replace = str_extract(gt[sample, mhap], pattern_to_replace)
                replacement_pattern = gsub(alleles[allele], allele, pattern_to_replace)
                
                if(!grepl('^_', pattern_to_replace)){
                  pattern_to_replace = paste0('^', pattern_to_replace)
                }
                
                if(!grepl('_$', pattern_to_replace)){
                  pattern_to_replace = paste0(pattern_to_replace, '$')
                }
                
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
    
    for(sample in comparison_df[['Sample']]){
      
      if(monoclonals[which(comparison_df[['Sample']] == sample)]){
        temp_gt_polyclonal = matrix(gt[sample,],
                                    ncol = ncol(gt),
                                    nrow = 1,
                                    dimnames = list(
                                      sample,
                                      colnames(gt)
                                    ))
        
      }else if(polyclonals[which(comparison_df[['Sample']] == sample)]){
        
        poly_sample = which(comparison_df[['Sample']] == sample)
        
        if(is.null(max_coi)){
          sample_clonality = clonality[poly_sample]
        }else{
          sample_clonality = max_coi
        }
        
        temp_gt_polyclonal = matrix(NA, 
                                    nrow = sample_clonality,
                                    ncol = ncol(gt),
                                    dimnames = list(
                                      paste0(comparison_df[['Sample']][poly_sample], '_C',1:sample_clonality),
                                      colnames(gt)
                                    ))
        
        alleles = strsplit(gt[poly_sample,], '_')
        
        
        
        for(clone in 1:sample_clonality){
          
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
    
    coded_haplotypes = data.frame(Marker = colnames(gt_final),
                                  t(gt_final))
    
    colnames(coded_haplotypes) = c('Marker', rownames(gt_final))
    
    coded_haplotypes = left_join(coded_haplotypes, 
                                 ampseq_object@markers[,c('amplicon', 'chromosome', 'pos')],
                                 by = join_by('Marker' == 'amplicon'))
    
    coded_haplotypes %<>% pivot_longer(cols = all_of(rownames(gt_final)), names_to = 'Haplotype', values_to = 'Allele')
    
    
    
    coded_haplotypes$Haplotype = 
      factor(coded_haplotypes$Haplotype,
             levels = rownames(gt_final))
    
    coded_haplotypes$Marker = 
      factor(coded_haplotypes$Marker,
             levels = unique(coded_haplotypes$Marker))
    
    coded_haplotypes$chromosome = gsub('(^[A-z]+\\d+_|_v\\d+$)', '', coded_haplotypes$chromosome)
    
    
    coded_haplotypes$Allele = as.character(coded_haplotypes$Allele)
    
 
    coded_haplotypes$Comparison = NA
    
    for(Sample in comparison_df$Sample){
      
      coded_haplotypes$Comparison[
        grepl(Sample, coded_haplotypes$Haplotype)
      ] = comparison_df$Strata[comparison_df$Sample == Sample]
      
    }  
    
    haplotypes_plot = coded_haplotypes %>%
      filter(!is.na(Allele))%>%
      ggplot(aes(x = Marker, y = Haplotype, fill = Allele)) +
      geom_tile()+
      scale_fill_brewer(palette = 'Set1')+
      facet_grid(Comparison ~ chromosome, scales = 'free', space = 'free', switch = 'y')+
      theme(axis.text.x = element_blank(),
            strip.placement = 'outside')
    
  }
  
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


pairwise_hmmIBD = function(obj = NULL, parallel = TRUE, w = 1, n = 1){
  library(parallel)
  library(doMC)
  library(svMisc)
  
  if(class(obj) == 'ampseq'){
    loci_object = ampseq2loci(ampseq_object)
  }else if(class(obj) == 'loci'){
    loci_object = obj
  }else{
    stop("This function requires an object of the class loci or ampseq")
  }
  
  loci_table = loci_object@loci_table
  freq_table = loci_object@freq_table
  markers = loci_object@markers
  
  polymorphic_sites = which(rowSums(freq_table == 1) == 0)
  
  loci_table = loci_table[,polymorphic_sites]
  freq_table = freq_table[polymorphic_sites,]
  markers = markers[polymorphic_sites,]
  
  pairs = as.data.frame(t(combn(rownames(loci_table), 2)))
  
  s = round(seq(1,nrow(pairs)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  pairs = pairs[low:high,]
  
  fx_get_relatedness = function(Yi, Yj, freq_table, markers){
    
    Ys = cbind(Yi, Yj)
    
    rownames(Ys) = colnames(loci_table)
    
    fs = freq_table[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    
    no_na_markers = markers[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    
    Ys =Ys[!is.na(Ys[,"Yi"]) & !is.na(Ys[,"Yj"]),]
    fs = fs[Ys[,1] <= 9 & Ys[,2] <= 9,]
    no_na_markers = no_na_markers[Ys[,1] <= 9 & Ys[,2] <= 9,]
    Ys = Ys[Ys[,1] <= 9 & Ys[,2] <= 9,]
    
    no_na_markers[["distance"]] = Inf
    no_na_markers = no_na_markers[fs[,1] != 1,]
    Ys = Ys[fs[,1] != 1,]
    fs = fs[fs[,1] != 1,]
    
    for(chromosome in levels(as.factor(no_na_markers[["chromosome"]]))){
      for(amplicon in 1:(nrow(no_na_markers[no_na_markers[["chromosome"]] == chromosome,])-1)){
        no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon, "distance"] = 
          no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - no_na_markers[no_na_markers[["chromosome"]] == chromosome,][amplicon, "pos"]
      }
    }
    
    if(sum(is.na(no_na_markers$distance)) > 0){
      no_na_markers[is.na(no_na_markers$distance),][['distance']] = Inf  
    }
    
    if(sum(no_na_markers$distance < 0) > 0){
      no_na_markers$distance[no_na_markers$distance < 0 ] = Inf
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
      
      Yi_id = pairs[pair, 1]
      Yj_id = pairs[pair, 2]
      
      Yi = split_clones(loci_table[Yi_id,], Yi_id)
      
      Yi = matrix(as.integer(Yi), 
                  nrow = nrow(Yi),
                  ncol = ncol(Yi),
                  dimnames = list(
                    rownames(Yi),
                    colnames(Yi)
                  ))
      
      Yj = split_clones(loci_table[Yj_id,], Yj_id)
      
      
      Yj = matrix(as.integer(Yj), 
                  nrow = nrow(Yj),
                  ncol = ncol(Yj),
                  dimnames = list(
                    rownames(Yj),
                    colnames(Yj)
                  ))
      
      estimate = NULL
      
      for(haplotype_i in 1:nrow(Yi)){
        for(haplotype_j in 1:nrow(Yj)){
          
          if(sum(is.na(Yi[haplotype_i,] - Yj[haplotype_j,]))/length(Yi[haplotype_i,]) < 1){
            
            estimate = rbind(estimate, fx_get_relatedness(Yi[haplotype_i,], Yj[haplotype_j,], freq_table, markers))
            
          }else{
            
            estimate = c(NA,NA)
            names(estimate) = c('khat', 'rhat')
            
          }
          
        }
      }
      
      if(sum(is.na(estimate[,'rhat'])) > 0){
        
        estimate = estimate[1,]
        
      }else{
        
        estimate = estimate[which.min(estimate[,'rhat']),]
        
      }
      
      
      data.frame(Yi = Yi_id, Yj = Yj_id, t(estimate))
    }
    
  }else{
    
    pairwise_df = NULL
    
    for(pair in 1:nrow(pairs)){
      
      Yi_id = pairs[pair, 1]
      Yj_id = pairs[pair, 2]
      
      Yi = split_clones(loci_table[Yi_id,], Yi_id)
      
      Yi = matrix(as.integer(Yi), 
                  nrow = nrow(Yi),
                  ncol = ncol(Yi),
                  dimnames = list(
                    rownames(Yi),
                    colnames(Yi)
                  ))
      
      Yj = split_clones(loci_table[Yj_id,], Yj_id)
      
      Yj = matrix(as.integer(Yj), 
                  nrow = nrow(Yj),
                  ncol = ncol(Yj),
                  dimnames = list(
                    rownames(Yj),
                    colnames(Yj)
                  ))
      
      estimate = NULL
      
      for(haplotype_i in 1:nrow(Yi)){
        for(haplotype_j in 1:nrow(Yj)){
          
          if(sum(is.na(Yi[haplotype_i,] - Yj[haplotype_j,]))/length(Yi[haplotype_i,]) < 1){
            
            estimate = rbind(estimate, fx_get_relatedness(Yi[haplotype_i,], Yj[haplotype_j,], freq_table, markers))
            
          }else{
            
            estimate = c(NA,NA)
            names(estimate) = c('khat', 'rhat')
            
          }
          
        }
      }
      
      if(sum(is.na(estimate[,'rhat'])) > 0){
        
        estimate = estimate[1,]
        
      }else{
        
        estimate = estimate[which.min(estimate[,'rhat']),]
        
      }
      
      progress(round(100*pair/nrow(pairs)))
      
      pairwise_df = rbind(pairwise_df, data.frame(Yi = Yi_id, Yj = Yj_id, t(estimate)))
      
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

### Relatedness_heatmap----



sort_Long_matrix = function(pairwise_relatedness_l, Var1, Var2){
  
  pairwise_relatedness_l = pairwise_relatedness_l[,c('Yi', 'Yj', 'khat', 'rhat', Var1, Var2)]
  
  pairwise_relatedness_l_sorted = NULL
  
  comparisons = cbind(rbind(
    sort(unique(c(unique(pairwise_relatedness_l[[Var1]]), unique(pairwise_relatedness_l[[Var2]])))),
    sort(unique(c(unique(pairwise_relatedness_l[[Var1]]), unique(pairwise_relatedness_l[[Var2]]))))),
    combn(sort(unique(c(unique(pairwise_relatedness_l[[Var1]]), unique(pairwise_relatedness_l[[Var2]])))), 2))
  
  for(comparison in 1:ncol(comparisons)){
    var1 = comparisons[,comparison][1]
    var2 = comparisons[,comparison][2]
    if(var1 == var2){
      
      pairwise_relatedness_l_sorted = rbind(
        pairwise_relatedness_l_sorted,
        pairwise_relatedness_l[pairwise_relatedness_l[[Var1]] == var1 &
                                 pairwise_relatedness_l[[Var2]] == var2, ])
      
    }else{
      
      temp1 = pairwise_relatedness_l[pairwise_relatedness_l[[Var1]] == var1 &
                                       pairwise_relatedness_l[[Var2]] == var2, ]
      
      temp2 = pairwise_relatedness_l[pairwise_relatedness_l[[Var1]] == var2 &
                                       pairwise_relatedness_l[[Var2]] == var1, ]
      
      Yi = temp2[['Yj']]
      Population_Var1 = temp2[[Var2]]
      
      Yj = temp2[['Yi']]
      Population_Var2 = temp2[[Var1]]
      
      temp2[['Yi']] = Yi
      temp2[[Var1]] = Population_Var1
      
      temp2[['Yj']] = Yj
      temp2[[Var2]] = Population_Var2
      
      pairwise_relatedness_l_sorted = rbind(pairwise_relatedness_l_sorted,
                                            temp1,
                                            temp2)
      
    }
    
  }
  
  return(pairwise_relatedness_l_sorted)
  
}

### plot_network----

plot_network = function(pairwise_relatedness,  
                           threshold,
                           metadata,
                           sample_id,
                           group_by,
                           levels,
                           colors,
                        vertex.size = 4,
                        method = 'mst'){
  
  if(sum(is.na(metadata[[group_by]])) > 0){
    metadata[is.na(metadata[[group_by]]),][[group_by]] = 'missing data'
  }
  
  pairwise_relatedness_l = pairwise_relatedness
  
  pairwise_relatedness_l %<>% filter(Yi %in% metadata$Sample_id,
                                Yj %in% metadata$Sample_id)
  
  variable = 'rhat'
  cols = c("Yi", "Yj")
  
  library(igraph)
  if(!require(GGally)){
    install.packages("GGally")
    library(GGally)
  }else{library(GGally)}
  
  if(!require(network)){
    install.packages("network")
    library(network)
  }else{library(network)}
  
  if(!require(sna)){
    install.packages("sna")
    library(sna)
  }else{library(sna)}
  
  
  if(!require(ggnet)){
    devtools::install_github("briatte/ggnet")
    library(ggnet)
  }else{library(ggnet)}
  
  
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
  
  
  if(method %in% c('msn', 'mst')){
    
    network_object <- mst(network_object)
    
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
    
  }else{
    
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
    
  }
  
  clusters = NULL
  for(node in 1:length(network_object)){
    
    # Set of samples linked with the node
    set_of_samples = c(names(network_object[[node]]), names(network_object[[node]][[1]]))
    
    proportion_of_samples_in_the_cluster = sapply(clusters, function(cluster){
      sum(set_of_samples %in% cluster)/length(set_of_samples)
    })
    
    cluster = which(proportion_of_samples_in_the_cluster > 0 & proportion_of_samples_in_the_cluster < 1 )
    
    if(length(cluster) == 1){
      
      clusters[[paste0('Cluster_', cluster)]] =
        c(clusters[[paste0('Cluster_', cluster)]],
          set_of_samples[!(set_of_samples %in% clusters[[paste0('Cluster_', cluster)]])]
        )
      
    }else if(length(cluster) > 1){
      # Combine clusters
      merged_cluster = unlist(clusters[cluster])
      
      names(merged_cluster) = NULL
      
      merged_cluster = c(merged_cluster,
                         set_of_samples[!(set_of_samples %in% merged_cluster)]
      )
      
      clusters[[cluster[1]]] = merged_cluster
      
      clusters[cluster[-1]] = NULL
      
      names(clusters) = paste0('Cluster_', 1:length(clusters))
      
      
    }else if(length(cluster) == 0 & isEmpty(proportion_of_samples_in_the_cluster)){
      cluster = length(clusters) + 1
      clusters[[paste0('Cluster_', cluster)]]  = set_of_samples
    }else if(length(cluster) == 0 & sum(proportion_of_samples_in_the_cluster == 1) == 0){
      cluster = length(clusters) + 1
      clusters[[paste0('Cluster_', cluster)]]  = set_of_samples
      
    }
    
  }
  
  
  cluster_df = NULL
  
  for(cluster in names(clusters)){
    
    cluster_name = cluster
    
    if(length(clusters[[cluster]]) == 1){
      cluster_name = gsub('Cluster', 'Singleton', cluster_name)
    }else if(length(clusters[[cluster]]) == 2){
      cluster_name = gsub('Cluster', 'Doubleton', cluster_name)
    }else if(length(clusters[[cluster]]) == 3){
      cluster_name = gsub('Cluster', 'Tripleton', cluster_name)
    }
    
    
    cluster_df = rbind(cluster_df, data.frame(Cluster = cluster_name, Sample_id = clusters[[cluster]]))
  }
  
  return(list(network_object = network_object,
              plot_network = plot_network,
              clusters = cluster_df))
  
}


### plot_ggnetwork----

plot_ggnetwork = function(pairwise_relatedness,  
                          threshold,
                          metadata,
                          sample_id,
                          color_by = NULL,
                          shape_by = NULL,
                          alpha_by = NULL,
                          alpha = NULL,
                          vertex.size = 4,
                          mode = 'fruchtermanreingold'){
  
  
  library(igraph)
  if(!require(GGally)){
    install.packages("GGally")
    library(GGally)
  }else{library(GGally)}
  
  if(!require(network)){
    install.packages("network")
    library(network)
  }else{library(network)}
  
  if(!require(sna)){
    install.packages("sna")
    library(sna)
  }else{library(sna)}
  
  if(!require(ggnet)){
    devtools::install_github("briatte/ggnet")
    library(ggnet)
  }else{library(ggnet)}
  
  
  if(!is.null(color_by)){
    if(sum(is.na(metadata[[color_by]])) > 0){
      metadata[is.na(metadata[[color_by]]),][[color_by]] = 'missing data'
    }
  }
  
  if(!is.null(shape_by)){
    if(sum(is.na(metadata[[shape_by]])) > 0){
      metadata[is.na(metadata[[shape_by]]),][[shape_by]] = 'missing data'
    }
  }
  
  
  if(!is.null(alpha_by)){
    if(sum(is.na(metadata[[alpha_by]])) > 0){
      metadata[is.na(metadata[[alpha_by]]),][[alpha_by]] = 'missing data'
    }
  }
  
  
  pairwise_relatedness_l = pairwise_relatedness
  
  pairwise_relatedness_l %<>% filter(Yi %in% metadata$Sample_id,
                                     Yj %in% metadata$Sample_id)
  
  pairwise_relatedness_matrix = matrix(data = NA,
                                       ncol = nrow(metadata),
                                       nrow = nrow(metadata),
                                       dimnames = list(metadata$Sample_id,
                                                       metadata$Sample_id))
  
  for(pair in 1:nrow(pairwise_relatedness_l)){
    
    pairwise_relatedness_matrix[pairwise_relatedness_l[pair,][['Yi']],
                                pairwise_relatedness_l[pair,][['Yj']]] = 
      pairwise_relatedness_l[pair,][['rhat']]
    
    pairwise_relatedness_matrix[pairwise_relatedness_l[pair,][['Yj']],
                                pairwise_relatedness_l[pair,][['Yi']]] = 
      pairwise_relatedness_l[pair,][['rhat']]
  }
  
  pairwise_relatedness_matrix[pairwise_relatedness_matrix < threshold] = 0
  
  
  relatedness_network = network(pairwise_relatedness_matrix, directed = FALSE)
  
  network.vertex.names(relatedness_network) = metadata$Sample_id
  
  
  if(!is.null(color_by)){
    relatedness_network %v% color_by = metadata[[color_by]]
  }else{
    color_by = 'gray75'
  }
  
  if(!is.null(shape_by)){
    relatedness_network %v% shape_by = metadata[[shape_by]]  
  }else{
    shape_by = 16
  }
  
  if(!is.null(alpha_by)){
    relatedness_network %v% alpha_by = metadata[[alpha_by]]
  }else if(!is.null(alpha)){
    alpha_by = alpha
  }else{
    alpha_by = 1
  }
  
  
  plot_network = ggnet2(relatedness_network,
                        size = vertex.size,
                        color = color_by,
                        shape = shape_by,
                        palette = 'Paired',
                        alpha = alpha_by,
                        mode = mode)

  
  return(list(network_object = relatedness_network,
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

GRM_evectors = function(gt = NULL, dist_table = NULL, k = nrow(gt), metadata, Pop = 'Population', q = 2){
  
  if(!is.null(gt)){
    
    X = t(gsub('_.+$','',gt))
    
    X = matrix(as.numeric(X), ncol = ncol(X),
               nrow = nrow(X), 
               dimnames = list(
                 rownames(X),
                 colnames(X)
               ))
    
    X[is.na(X)] = 0
    
    
  }else if(!is.null(dist_table)){
    
    
    dist_matrix = matrix(data = NA,
                         ncol = nrow(metadata),
                         nrow = nrow(metadata),
                         dimnames = list(metadata$Sample_id,
                                         metadata$Sample_id))
    
    # #### Add metadata to the PCA
    # Pop_col = merge(data.frame(Sample_id = gsub('_C[1,2]$','',colnames(dist_matrix)),
    #                            order = 1:ncol(dist_matrix)), metadata[,c('Sample_id', Pop)], 
    #                 by = 'Sample_id', all.y = T)
    # 
    # Pop_col = Pop_col[order(Pop_col$order),]
    
    
    dist_table %<>% filter(Yi %in% metadata$Sample_id & Yj %in% metadata$Sample_id)
    
    
    for(pair in 1:nrow(dist_table)){
      
      dist_matrix[dist_table[pair,][['Yi']],
                  dist_table[pair,][['Yj']]] = 
        dist_table[pair,][['euDist']]
      
      dist_matrix[dist_table[pair,][['Yj']],
                  dist_table[pair,][['Yi']]] = 
        dist_table[pair,][['euDist']]
    }
    
    dist_matrix[is.na(dist_matrix)] = 1
    
    X = dist_matrix
    
  }
  
  
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
IBD_evectors = function(ampseq_object, relatedness_table, k = NULL, Pop = 'Population', q = 2,
                        method = 'princomp', # fastGRM
                        cor = T
                        ){
  
  metadata = ampseq_object@metadata
  
  
  
  pairwise_relatedness_matrix = matrix(data = NA,
                                       ncol = nrow(ampseq_object@metadata),
                                       nrow = nrow(ampseq_object@metadata),
                                       dimnames = list(ampseq_object@metadata$Sample_id,
                                                       ampseq_object@metadata$Sample_id))
  
  #### Add metadata to the PCA
  Pop_col = merge(data.frame(Sample_id = gsub('_C[1,2]$','',colnames(pairwise_relatedness_matrix)),
                             order = 1:ncol(pairwise_relatedness_matrix)), metadata[,c('Sample_id', Pop)], 
                  by = 'Sample_id', all.y = T)
  
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
  
  if(method == 'fastGRM'){
    
    
    evector = fastSVDCpp(pairwise_relatedness_matrix, k, q)

    for(i in 1:k){
      evector[,i] = sign(evector[1,i])*evector[,i]
    }

    evalues = NULL

    for(i in 1:k){
      evalues = c(evalues, unlist((pairwise_dist_matrix %*% evector[,i])/evector[,i])[1])
    }

    contrib = 100*(evalues)^2/sum((evalues)^2)

    evector = data.frame(Pop_col, evector)
    names(evector) = c(colnames(Pop_col), paste0(rep('PC', k), 1:k))

    ibd_pca = list(eigenvector=evector, eigenvalues = evalues, contrib = contrib)

  }else if(method == 'princomp'){
    
    ## Using princomp from R
    
    ibd_pca = princomp(pairwise_dist_matrix, cor = cor)
    ibd_evector = ibd_pca$scores
    ibd_evalues = ibd_pca$sdev
    
    ibd_contrib = 100*(ibd_evalues)^2/sum((ibd_evalues)^2)
    
    ibd_evector = data.frame(Pop_col, ibd_evector)
    names(ibd_evector) = c(colnames(Pop_col), paste0(rep('PC', k), 1:k))
    
    ibd_pca = list(eigenvector=ibd_evector, eigenvalues = ibd_evalues, contrib = ibd_contrib)
    
  }else{
    print("No Method detected")
  }

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

## get_gene_description ----

get_gene_description = function(obj = NULL, gff = NULL){
  
  if(class(obj) == 'ampseq'){
    data = obj@markers
  }else{
    data = obj
  }
  
  ref_gff = ape::read.gff(gff)
  coding_regions = ref_gff[grepl('gene', ref_gff$type)&
                             !grepl('^Transfer',ref_gff$seqid),
                           c('seqid', 'start', 'end', 'attributes')]
  
  coding_regions$gene_id = gsub('ID=','',str_extract(coding_regions$attributes, 'ID=(PVP01|PF3D7)_([0-9]+|MIT[0-9]+|API[0-9]+)'))
  
  coding_regions$gene_description = gsub('description=','',str_extract(coding_regions$attributes, '(description=.+$|description=.+;)'))
  
  coding_regions$gene_name = gsub('(Name=|;)','',str_extract(coding_regions$attributes, 'Name=\\w+;'))
  
  coding_regions$gene_name = ifelse(is.na(coding_regions$gene_name), coding_regions$gene_id, coding_regions$gene_name)
  
  coding_regions = coding_regions[order(coding_regions$start),]
  coding_regions = coding_regions[order(coding_regions$seqid),]
  rownames(coding_regions) = 1:nrow(coding_regions)
  
  
  data$gene_id = NA
  data$gene_name = NA
  data$gene_description = NA
  
  
  for(gene in 1:nrow(coding_regions)){
    
    if(nrow(data[data$chromosome == coding_regions[gene, ][['seqid']]&
                 data$start >= coding_regions[gene, ][['start']]&
                 data$start <= coding_regions[gene, ][['end']],]) != 0){
      
      data[data$chromosome == coding_regions[gene, ][['seqid']]&
             data$start >= coding_regions[gene, ][['start']]&
             data$start <= coding_regions[gene, ][['end']],][['gene_id']] = coding_regions[gene, ][['gene_id']]
      
      data[data$chromosome == coding_regions[gene, ][['seqid']]&
             data$start >= coding_regions[gene, ][['start']]&
             data$start <= coding_regions[gene, ][['end']],][['gene_name']] = coding_regions[gene, ][['gene_name']]
      
      data[data$chromosome == coding_regions[gene, ][['seqid']]&
             data$start >= coding_regions[gene, ][['start']]&
             data$start <= coding_regions[gene, ][['end']],][['gene_description']] = coding_regions[gene, ][['gene_description']]
    }
    
  }
  
  # indels that start out of the gene region but end 
  
  for(pos in rownames(data[is.na(data$gene_id), ])){
    
    if(nrow(coding_regions[coding_regions[['seqid']] == data[pos,][['chromosome']] &
                           coding_regions[['start']] <= data[pos,][['end']] &
                           coding_regions[['end']] >= data[pos,][['end']],]) != 0){
      data[pos,][['gene_id']] = coding_regions[coding_regions[['seqid']] == data[pos,][['chromosome']] &
                                                 coding_regions[['start']] <= data[pos,][['end']] &
                                                 coding_regions[['end']] >= data[pos,][['end']],][['gene_id']]
      
      data[pos,][['gene_name']] = coding_regions[coding_regions[['seqid']] == data[pos,][['chromosome']] &
                                                          coding_regions[['start']] <= data[pos,][['end']] &
                                                          coding_regions[['end']] >= data[pos,][['end']],][['gene_name']]
      
      data[pos,][['gene_description']] = coding_regions[coding_regions[['seqid']] == data[pos,][['chromosome']] &
                                                          coding_regions[['start']] <= data[pos,][['end']] &
                                                          coding_regions[['end']] >= data[pos,][['end']],][['gene_description']]
    }
    
  }
  
  return(data[, c('gene_id', 'gene_name', 'gene_description')])
}

remove_replicates = function(ampseq_object, v){
  
  gt = ampseq_object@gt
  metadata = ampseq_object@metadata
  
  # Remove duplicates

  duplicated_samples = metadata[duplicated(metadata[[v]]),][[v]]
  
  duplicated_samples = unique(duplicated_samples)
  
  gt_replicates = gt[duplicated(metadata[[v]]),]
  gt = gt[!duplicated(metadata[[v]]),]
  
  metadata_replicates = metadata[duplicated(metadata[[v]]),]
  metadata = metadata[!duplicated(metadata[[v]]),]
  
  for(sample in duplicated_samples){
    
    temp_data = rbind(gt[metadata[metadata[[v]] == sample, ][['Sample_id']],],
                      gt_replicates[metadata_replicates[metadata_replicates[[v]] == sample, ][["Sample_id"]],])
    
    temp_metadata = rbind(metadata[metadata[[v]] == sample, ],
                          metadata_replicates[metadata_replicates[[v]]== sample, ])
    
    temp_test = rowSums(!is.na(temp_data), na.rm = T)
    
    gt[metadata[metadata[[v]] == sample, ][['Sample_id']],] = temp_data[which.max(temp_test),]
    metadata[metadata[[v]] == sample, ] = temp_metadata[which.max(temp_test),]
    
  }  
  
  ampseq_object@gt = gt
  ampseq_object@metadata = metadata
  
  
  cigars_gt = get_cigar_alleles(ampseq_object, 
                                from = 'gt',
                                as = 'vector')
  
  cigars_asvtab = get_cigar_alleles(ampseq_object, 
                                    from = 'asv_table',
                                    as = 'vector',
                                    col = 'CIGAR_masked')
  
  if(sum(!(cigars_gt %in% cigars_asvtab)) > 0 & 
     sum(!(cigars_asvtab %in% cigars_gt)) > 0){
    
    cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
               paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
               '\n\n',
               'Cigar strings in asv_table that are not present in gt matrix are:\n',
               paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
               '\n'
    ))
    
    stop('There are cigar strings that are different in the gt and the asv_table')
    
  }else if(sum(!(cigars_gt %in% cigars_asvtab)) > 0){
    
    cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
               paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
               '\n'))
    
    stop('There are cigar strings in the gt that are not present in the asv_table')
    
  }else if(sum(!(cigars_asvtab %in% cigars_gt)) > 0){
    
    cat(paste0('The following cigar strings in asv_table are not longer present in gt matrix:\n',
               paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
               '\n'))
    
    cigar_strings_to_remove = cigars_asvtab[!(cigars_asvtab %in% cigars_gt)]
    
    for(cigar_string_to_remove in unique(cigar_strings_to_remove)){
      
      Amplicon = gsub(';.+$', '', cigar_string_to_remove)
      CIGAR = gsub('^.+;', '', cigar_string_to_remove)
      
      ampseq_object@asv_table = 
        ampseq_object@asv_table[!(ampseq_object@asv_table[['Amplicon']] == Amplicon &
                                ampseq_object@asv_table[['CIGAR_masked']] == CIGAR),]
      
      ampseq_object@asv_seqs = ampseq_object@asv_seqs[
        names(ampseq_object@asv_seqs) %in% ampseq_object@asv_table$hapid]
      
      ampseq_object@asv_seqs_masked = ampseq_object@asv_seqs_masked[
        names(ampseq_object@asv_seqs_masked) %in% ampseq_object@asv_table$hapid]
      
      ampseq_object@asv_table$hapid = paste0('ASV', 1:nrow(ampseq_object@asv_table))
      names(ampseq_object@asv_seqs) = ampseq_object@asv_table$hapid
      names(ampseq_object@asv_seqs_masked) = ampseq_object@asv_table$hapid
      
    }
    
    
  }else{
    print('cigar strings are consistent between gt and asv_table')
  }
  
  consistency_between_gt_and_asvtab(ampseq_object)
  
  return(ampseq_object)
  
}


## rGenome2ampseq----
rGenome2ampseq = function(rGenome_object,
                          markers,
                          ref_seqs = '~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq249_refseqs.fasta'){
  
  
  markers[['amplicon']] = gsub('(-|:|/)', '_', markers[['amplicon']])
  
  ampseq_gt_table = matrix(NA,
                           ncol = length(markers[['amplicon']]),
                           nrow = ncol(rGenome_object@gt),
                           dimnames = list(
                             colnames(rGenome_object@gt),
                             markers[['amplicon']]
                           )
                           )
  
  
  amplicons_ref_seqs = Biostrings::readDNAStringSet(ref_seqs)
  
  names(amplicons_ref_seqs) = gsub('(-|:|/)', '_', names(amplicons_ref_seqs))
  
  for(amplicon in markers$amplicon){
    
    amplicon_start = markers[markers[['amplicon']] == amplicon,][['start']]
    amplicon_end = markers[markers[['amplicon']] == amplicon,][['end']]
    amplicon_chromosome = markers[markers[['amplicon']] == amplicon,][['chromosome']]
    amplicon_positions = amplicon_start:amplicon_end
    
    rGenome_positions = which(rGenome_object@loci_table[['CHROM']] == amplicon_chromosome &
            rGenome_object@loci_table[['POS']] %in% amplicon_positions)
    
    if(length(rGenome_positions) > 0){
      
      if(length(rGenome_positions) == 1){
        
        rGenome_gt_amplicon = 
          matrix(rGenome_object@gt[rGenome_object@loci_table[['CHROM']] == amplicon_chromosome &
                                     rGenome_object@loci_table[['POS']] %in% amplicon_positions,],
                 nrow = length(rGenome_positions),
                 dimnames = list(
                   rownames(rGenome_object@loci_table[rGenome_object@loci_table[['CHROM']] == amplicon_chromosome &
                                                rGenome_object@loci_table[['POS']] %in% amplicon_positions,]),
                   names(rGenome_object@gt[rGenome_object@loci_table[['CHROM']] == amplicon_chromosome &
                                                rGenome_object@loci_table[['POS']] %in% amplicon_positions,])
                 )
          )
        
      }else if(length(rGenome_positions) >1){
        rGenome_gt_amplicon = 
          matrix(rGenome_object@gt[rGenome_object@loci_table[['CHROM']] == amplicon_chromosome &
                                     rGenome_object@loci_table[['POS']] %in% amplicon_positions,],
                 nrow = length(rGenome_positions),
                 dimnames = list(
                   rownames(rGenome_object@gt[rGenome_object@loci_table[['CHROM']] == amplicon_chromosome &
                                                rGenome_object@loci_table[['POS']] %in% amplicon_positions,]),
                   colnames(rGenome_object@gt[rGenome_object@loci_table[['CHROM']] == amplicon_chromosome &
                                                rGenome_object@loci_table[['POS']] %in% amplicon_positions,])
                 )
          )
        
      }
      
      
      rGenome_loci_table_amplicon = rGenome_object@loci_table[rGenome_object@loci_table[['CHROM']] == amplicon_chromosome &
                                                                rGenome_object@loci_table[['POS']] %in% amplicon_positions, c('CHROM', 'POS', 'REF', 'ALT')]
      
      
      
      rGenome_gt_amplicon_alleles = gsub(':\\d+', '', rGenome_gt_amplicon)
      rGenome_gt_amplicon_rd = gsub('\\d+:', '', rGenome_gt_amplicon)
      
      rGenome_gt_amplicon[,Sample_id]
      rGenome_gt_amplicon_alleles[,Sample_id]
      rGenome_gt_amplicon_rd[,Sample_id]
      # 
      # amplicon_position = 2442390
      # amplicon_position = 459959
      
      # rGenome_gt_amplicon_alleles[grepl(amplicon_position, rownames(rGenome_gt_amplicon_alleles)),]
      
      for(amplicon_position in rGenome_loci_table_amplicon$POS){ # Replace alternative alleles
        
        alt_alleles = strsplit(rGenome_loci_table_amplicon[rGenome_loci_table_amplicon$POS == amplicon_position, ][['ALT']], ',')[[1]]
        ref_allele = rGenome_loci_table_amplicon[rGenome_loci_table_amplicon$POS == amplicon_position, ][['REF']]
        
        #alt_allele = 1
        
        for(alt_allele in 1:length(alt_alleles)){
          
          if(nchar(ref_allele) == 1 & nchar(alt_alleles[alt_allele]) == 1 & alt_alleles[alt_allele] != '*'){# it is an SNV
            
            # For monoclonals
            rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                        grepl(paste0('^', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
            ] =
              gsub(paste0('^', alt_allele, '$'), 
                   paste0(amplicon_position - amplicon_start + 1, alt_alleles[alt_allele]), 
                   rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                               grepl(paste0('^', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
                   ])
            
            # For first clone
            rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                        grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
            ] =
              gsub(paste0('^', alt_allele, '/'), 
                   paste0(amplicon_position - amplicon_start + 1, alt_alleles[alt_allele], '/'), 
                   rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                               grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
                   ])
            
            # For second clone
            rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                        grepl(paste0('/', alt_allele ,'$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
            ] =
              gsub(paste0('/', alt_allele, '$'), 
                   paste0('/', amplicon_position - amplicon_start + 1, alt_alleles[alt_allele]), 
                   rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                               grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
                   ])
            
            
          }else if(nchar(ref_allele) > nchar(alt_alleles[alt_allele]) | alt_alleles[alt_allele] == '*'){# if it is a deletion
            
            if(alt_alleles[alt_allele] == '*'){ # The current position is deleted
              
              # stop()
              # deletion_region = paste0('D=', ref_allele)
              # 
              # # For monoclonals
              # 
              # rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
              #                             grepl(paste0('^', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              # ] =
              #   gsub(paste0('^', alt_allele, '$'), 
              #        paste0(amplicon_position - amplicon_start + 1, deletion_region), 
              #        rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
              #                                    grepl(paste0('^', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])
              # 
              # # For first clone
              # rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
              #                             grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              # ] =
              #   gsub(paste0('^', alt_allele, '/'), 
              #        paste0(amplicon_position - amplicon_start + 1, deletion_region, '/'), 
              #        rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
              #                                    grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])
              # 
              #               
              # # For second clone
              # rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
              #                             grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              # ] =
              #   gsub(paste0('/', alt_allele, '$'), 
              #        paste0('/', amplicon_position - amplicon_start + 1, deletion_region), 
              #        rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
              #                                    grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])
              
              
              rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                          grepl(paste0('^', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ] =
                gsub(paste0('^', alt_allele, '$'),
                     '0',
                     rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                                 grepl(paste0('^', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])

              # For first clone
              
              read1 = gsub('/\\d+', '', rGenome_gt_amplicon_rd[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                     grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ])
              
              read2 = gsub('\\d+/', '', rGenome_gt_amplicon_rd[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                                             grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ])
              
              
              rGenome_gt_amplicon_rd[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                          grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ] = as.integer(read1) + as.integer(read2)
              
              rm(read1)
              rm(read2)
              
              rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                          grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ] =
                gsub(paste0('^', alt_allele, '/'),
                     '0/',
                     rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                                 grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])

              
              
              # For second clone
              
              read1 = gsub('/\\d+', '', rGenome_gt_amplicon_rd[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                          grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ])
              
              read2 = gsub('\\d+/', '', rGenome_gt_amplicon_rd[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                                               grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ])
              
              
              rGenome_gt_amplicon_rd[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                          grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ] = as.integer(read1) + as.integer(read2)
              
              rm(read1)
              rm(read2)
              
              rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                          grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ] =
                gsub(paste0('/', alt_allele, '$'),
                     '/0',
                     rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                                 grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])
              
              # Edit homozygous sites
              
              rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                          grepl('0/0',rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ] = '0'
              
            }else{ # The next position is deleted
              
              deletion_region = paste0('D=', gsub(alt_alleles[alt_allele], '',ref_allele))
              
              # For monoclonals
              rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                          grepl(paste0('^', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ] =
                gsub(paste0('^', alt_allele, '$'), 
                     paste0(amplicon_position - amplicon_start + 2, deletion_region), 
                     rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                                 grepl(paste0('^', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])
     
              # For first clone
              rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                          grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ] =
                gsub(paste0('^', alt_allele, '/'), 
                     paste0(amplicon_position - amplicon_start + 2, deletion_region, '/'), 
                     rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                                 grepl(paste0('^', alt_allele, '/'), rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])
              
                       
              # For second clone
              rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                          grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
              ] =
                gsub(paste0('/', alt_allele, '$'), 
                     paste0('/', amplicon_position - amplicon_start + 2, deletion_region), 
                     rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                                 grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])
              
            }
            
            
          }else if(nchar(ref_allele) == nchar(alt_alleles[alt_allele]) & alt_alleles[alt_allele] != '*'){
            
            ref_allele2 = strsplit(ref_allele, '')[[1]]
            
            alt_allele2 = strsplit(alt_alleles[alt_allele], '')[[1]]
            
            polymorphic_sites = which(ref_allele2 != alt_allele2)
            
            alt_allele_replacement = NULL
            
            # site = 1
            
            for(site in polymorphic_sites){
              
              alt_allele_replacement = paste0(alt_allele_replacement,
                                              paste0(amplicon_position - amplicon_start + site, alt_allele2[site])
                                              )
              
            }
            
            
            # For monoclonals
            rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                        grepl(paste0('^', alt_allele, '$'),
                                              rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
            ] =
              gsub(paste0('^', alt_allele, '$'), 
                   alt_allele_replacement, 
                   rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                               grepl(paste0('^', alt_allele, '$'),
                                                     rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
                   ])
            
            # For first clone
            rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                        grepl(paste0('^', alt_allele, '/'),
                                              rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
            ] =
              gsub(paste0('^', alt_allele, '/'), 
                   paste0(alt_allele_replacement, '/'), 
                   rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                               grepl(paste0('^', alt_allele, '/'),
                                                     rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
                   ])
            
            # For second clone
            rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                        grepl(paste0('/', alt_allele ,'$'),
                                              rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
            ] =
              gsub(paste0('/', alt_allele, '$'), 
                   paste0('/', alt_allele_replacement), 
                   rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                               grepl(paste0('/', alt_allele, '$'),
                                                     rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
                   ])
            
            
            
          }else if(nchar(ref_allele) < nchar(alt_alleles[alt_allele])){# if it is an insertion
            
            #inserted_region = paste0('I=', gsub(paste0('^', ref_allele), '',alt_alleles[alt_allele]))
            inserted_region = paste0('I=', alt_alleles[alt_allele])
            
            # For monoclonals
            rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                        grepl(paste0('^', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
            ] =
              gsub(paste0('^', alt_allele, '$'), 
                   paste0(amplicon_position - amplicon_start + 1,# + nchar(ref_allele), 
                          inserted_region), 
                   rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                               grepl(paste0('^', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])
            
            # For first clone
            rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                        grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
            ] =
              gsub(paste0('^', alt_allele, '/'), 
                   paste0(amplicon_position - amplicon_start + 1,# + nchar(ref_allele), 
                          inserted_region, '/'), 
                   rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                               grepl(paste0('^', alt_allele, '/'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])
            
            # For second clone
            rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                        grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])
            ] =
              gsub(paste0('/', alt_allele, '$'), 
                   paste0('/', amplicon_position - amplicon_start + 1, # + nchar(ref_allele), 
                          inserted_region), 
                   rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),
                                               grepl(paste0('/', alt_allele, '$'),rGenome_gt_amplicon_alleles[which(rGenome_loci_table_amplicon$POS == amplicon_position),])])
            
          }
          
        }
        
      }
      
      amplicon_cigarstrings = NULL

      for(Sample_id in colnames(rGenome_gt_amplicon_alleles)){

        
        sample_positions = rGenome_gt_amplicon_alleles[,Sample_id]
        
        if(sum(is.na(sample_positions)) == 0){ # if samples do not contain missing information
          
          if(sum(grepl('/', sample_positions)) == 0){ # monoclonal samples
            
            if(sum(sample_positions != 0) > 0){# samples with mutations respect to reference strain
              
              sample_polymorphic_positions_alleles = sample_positions[sample_positions != 0]
              
              sample_polymorphic_positions = unlist(str_extract_all(sample_polymorphic_positions_alleles, '\\d+'))
              
              if(sum(duplicated(sample_polymorphic_positions)) > 0){
                
                duplicated_positions = sample_polymorphic_positions[duplicated(sample_polymorphic_positions)]
                
                for(duplicated_position in duplicated_positions){
                  
                  duplicated_positions_alleles = sample_polymorphic_positions_alleles[grepl(paste0('^', duplicated_position, '[DIAGTC]'), sample_polymorphic_positions_alleles)]
                  
                  if(sum(grepl('D', duplicated_positions_alleles)) == length(duplicated_positions_alleles)){
                    
                    correct_allele = paste0(duplicated_positions_alleles[1], gsub('^\\d+D=','',duplicated_positions_alleles[-1]))
                    
                    sample_polymorphic_positions_alleles = sample_polymorphic_positions_alleles[!(duplicated(sample_polymorphic_positions) & 
                                                                                                  grepl(paste0('^', duplicated_position, '[DIAGTC]'), sample_polymorphic_positions_alleles))]
                    
                    if(length(correct_allele) != length(sample_polymorphic_positions_alleles[grepl(paste0('^', duplicated_position, '[DIAGTC]'), sample_polymorphic_positions_alleles)])){
                      stop('1')
                    }
                    
                    sample_polymorphic_positions_alleles[grepl(paste0('^', duplicated_position, '[DIAGTC]'), sample_polymorphic_positions_alleles)] = correct_allele
                    
                    sample_polymorphic_positions = unlist(str_extract_all(sample_polymorphic_positions_alleles, '\\d+'))
                    
                  }
                  
                }
                
              }
              
              if(sum(grepl("109I=TG116D=A116G141C", paste(paste0(paste(sample_polymorphic_positions_alleles, collapse = ''),
                                                     ":",
                                                     min(as.integer(rGenome_gt_amplicon_rd[,Sample_id]))), collapse = ''))) > 0){
                stop('Problematic string 1')
              }
              
              amplicon_cigarstrings = 
                c(amplicon_cigarstrings,
                  paste0(paste(sample_polymorphic_positions_alleles, collapse = ''),
                         ":",
                         min(as.integer(rGenome_gt_amplicon_rd[,Sample_id])))
                )
              
            }else{ # Samples equals to reference strain
              
              amplicon_cigarstrings = 
                c(amplicon_cigarstrings,
                  paste0('.',
                         ":",
                         min(as.integer(rGenome_gt_amplicon_rd[,Sample_id])))
                )
            }
            
          }else{ # samples with more than two alleles (polyclonal samples)
            
            clones = list(gsub('/.+', '', sample_positions),
                          gsub('.+/', '', sample_positions))
            
            clones_rds = list(gsub('/.+', '', rGenome_gt_amplicon_rd[,Sample_id]),
                              gsub('.+/', '', rGenome_gt_amplicon_rd[,Sample_id]))
            clone = 2
            clone_cigarstring = NULL
            
            for(clone in 1:length(clones)){
              
              if(sum(clones[[clone]] != 0) > 0){# clones with mutations respect to reference strain
                
                
                sample_polymorphic_positions_alleles = clones[[clone]][clones[[clone]] != 0] # Alleles present in all polymorphic positions in the clone of the sample
                
                sample_polymorphic_positions = unlist(str_extract_all(sample_polymorphic_positions_alleles, '\\d+')) # polymorphic positions in the clone of the sample
                
                if(sum(duplicated(sample_polymorphic_positions)) > 0){ # if there are duplicated positions
                  
                  duplicated_positions = sample_polymorphic_positions[duplicated(sample_polymorphic_positions)]
                  
                  for(duplicated_position in duplicated_positions){
                    
                    duplicated_positions_alleles = sample_polymorphic_positions_alleles[grepl(paste0('^', duplicated_position, '[DIAGTC]'), sample_polymorphic_positions_alleles)] # extract alleles in the selected duplicated position
                    
                    if(sum(grepl('D', duplicated_positions_alleles)) == length(duplicated_positions_alleles)){ # Duplication is because two or more consecutive deletions
                      
                      correct_allele = paste0(duplicated_positions_alleles[1], gsub('^\\d+D=', '', duplicated_positions_alleles[-1])) # Concatenate duplicated deletions
                      
                      sample_polymorphic_positions_alleles = sample_polymorphic_positions_alleles[!(duplicated(sample_polymorphic_positions) & # Delete the duplicated position
                                                                                                      grepl(paste0('^', duplicated_position, '[DIAGTC]'), 
                                                                                                            sample_polymorphic_positions_alleles))]
                      
                      if(length(correct_allele) != length(sample_polymorphic_positions_alleles[grepl(paste0('^', duplicated_position, '[DIAGTC]'), sample_polymorphic_positions_alleles)])){
                        stop('2')
                      }
                      
                      sample_polymorphic_positions_alleles[grepl(paste0('^', duplicated_position, '[DIAGTC]'), sample_polymorphic_positions_alleles)] = correct_allele # Replace the correct allele in the position
                      
                      sample_polymorphic_positions = unlist(str_extract_all(sample_polymorphic_positions_alleles, '\\d+'))
                      
                    }else if(grepl('D', duplicated_positions_alleles[1]) & grepl('I', duplicated_positions_alleles[2]) & length(duplicated_positions_alleles) == 2){# Duplication is because a deletion and then insertion
                      
                      allele_in_first_replicate = duplicated_positions_alleles[1]
                      
                      allele_in_second_replicate = duplicated_positions_alleles[2]
                      
                      if(nchar(gsub('\\d+D=','',allele_in_first_replicate)) == 1){ # If only one position is deleted
                        
                        remainning_nucleotides = gsub(gsub('\\d+D=','',allele_in_first_replicate), '', gsub('\\d+I=','',allele_in_second_replicate))
                        
                        if(nchar(remainning_nucleotides) == 1){ # and if only one position is inserted
                          
                          correct_allele = paste0(duplicated_position, remainning_nucleotides)
                          
                          sample_polymorphic_positions_alleles = sample_polymorphic_positions_alleles[!(duplicated(sample_polymorphic_positions) & # Delete the duplicated position
                                                                                                          grepl(paste0('^', duplicated_position, '[DIAGTC]'), sample_polymorphic_positions_alleles))]
                          
                          sample_polymorphic_positions_alleles[grepl(paste0('^', duplicated_position, 'D'), sample_polymorphic_positions_alleles)] = correct_allele # Replace the correct allele in the position
                          
                          sample_polymorphic_positions = unlist(str_extract_all(sample_polymorphic_positions_alleles, '\\d+'))
                          
                        }else{ # if 2 or more positions are inserted
                          
                          ## CHECK HERE
                          
                          anchor_position = as.integer(duplicated_position) - 1
                          
                          anchor_nucleotide = (Biostrings::subseq(amplicons_ref_seqs[[amplicon]],
                                 start = anchor_position,
                                 end = anchor_position))
                          
                          if(!(as.character(anchor_position) %in% sample_polymorphic_positions)){
                            
                            correct_allele = paste0(anchor_position, 'I=', anchor_nucleotide, remainning_nucleotides)
                            
                            sample_polymorphic_positions_alleles = sample_polymorphic_positions_alleles[!(duplicated(sample_polymorphic_positions) & # Delete the duplicated position
                                                                                                            grepl(paste0('^', duplicated_position, '[DIAGTC]'), sample_polymorphic_positions_alleles))]
                            
                            sample_polymorphic_positions_alleles[grepl(paste0('^', duplicated_position, 'D'), sample_polymorphic_positions_alleles)] = correct_allele # Replace the correct allele in the position
                            
                            sample_polymorphic_positions = unlist(str_extract_all(sample_polymorphic_positions_alleles, '\\d+'))
                            
                          }else{
                            
                            stop('Position already exist')
                            
                          }
                          
                          
                        }
                        
                      }
                      
                      
                    }else if(sum(grepl('D', duplicated_positions_alleles[1])) > 0 & sum(!grepl('I|D', duplicated_positions_alleles[2])) > 0 & length(duplicated_positions_alleles) == 2){# Deletion and then substitution
                      
                      stop('4 or more duplicated sites with deletions and insertions')
                      
                    }else if(sum(grepl('D', duplicated_positions_alleles)) > 0 & sum(grepl('I', duplicated_positions_alleles)) > 0 & length(duplicated_positions_alleles) > 2){ 
                      
                      stop('3 or more duplicated sites with deletions and insertions')
                      
                    }
                    
                  }
                  
                }
                
                
                clone_cigarstring = 
                  c(clone_cigarstring,
                    paste0(paste(sample_polymorphic_positions_alleles, collapse = ''),
                           ":",
                           min(as.integer(clones_rds[[clone]]), na.rm = T))
                  )
                
              }else{ # clones equals to reference strain
                
                clone_cigarstring = 
                  c(clone_cigarstring,
                    paste0('.',
                           ":",
                           min(as.integer(clones_rds[[clone]]), na.rm = T))
                  )
              }
              
            }
            
            if(sum(grepl("109I=TG116D=A116G141C", paste0(clone_cigarstring, collapse = "_"))) > 0){
              stop('Problematic string 2')
            }
            
            amplicon_cigarstrings = c(amplicon_cigarstrings, paste0(clone_cigarstring, collapse = "_"))
            
          }
          
        }else{ # samples with at least one missing position
          amplicon_cigarstrings = c(amplicon_cigarstrings, NA)
        }
        
        
      } # Finish getting cigarstring for the amplicon in all samples
      
      if(sum(grepl("109I=TG116D=A116G141C", amplicon_cigarstrings)) > 0){
        stop('Problematic string 3')
      }
      
      ampseq_gt_table[, amplicon] = amplicon_cigarstrings
      
    }
    
  }
  
  asv_table = NULL
  
  for(amplicon in colnames(ampseq_gt_table)){
    
    cigar_strings = unique(unlist(strsplit(gsub(':\\d+','' , ampseq_gt_table[,amplicon]), '_')))
    
    cigar_strings = cigar_strings[!is.na(cigar_strings)]
    
    if(length(cigar_strings) > 0){
      
      asv_table = rbind(asv_table,
                        data.frame(hapid = paste0('ASV', (length(asv_table) + 1):(length(asv_table) + length(cigar_strings))),
                                   Amplicon = amplicon,
                                   CIGAR = cigar_strings,
                                   CIGAR_masked = cigar_strings
                        ))
      
    }
    
  }

  

  
  asv_table$haplength = NA
  asv_table$total_reads = NA
  asv_table$total_samples = NA
  asv_table$refid_PvP01 = asv_table$Amplicon
  asv_table$snv_dist_from_PvP01 = 0
  asv_table$indel_dist_from_PvP01 = 0
  asv_table$bimera = FALSE
  asv_table$strain = NA
  
  
  asv_table$hapid = paste0('ASV', 1:nrow(asv_table))
  
  asv_seqs = NULL
  
  #asv = 'ASV4780'
  
  for(asv in asv_table$hapid){
    
    amplicon = asv_table[asv_table$hapid == asv,][['Amplicon']]
    cigar_string = asv_table[asv_table$hapid == asv,][['CIGAR']]
    cigar_string_masked = asv_table[asv_table$hapid == asv,][['CIGAR_masked']]
    
    asv_table[asv_table$hapid == asv,][['total_samples']] = 
      sum(grepl(paste0('(^|_)',cigar_string_masked,':'), ampseq_gt_table[,amplicon]), na.rm = T)
    
    asv_table[asv_table$hapid == asv,][['total_reads']] =
      sum(as.integer(gsub('_.+$', '', gsub(paste0('(^|^.+_)',cigar_string_masked,':'),
                                           '',
                                           ampseq_gt_table[grepl(paste0('(^|_)',cigar_string_masked,':'),
                                                                 ampseq_gt_table[,amplicon]),amplicon]))))
    
    asv_seq = as.character(amplicons_ref_seqs[[amplicon]])
    
    if(cigar_string != '.'){
      
      positions = as.integer(str_extract_all(cigar_string, '\\d+')[[1]])
      variants = strsplit(cigar_string, '\\d+')[[1]]
      variants = variants[variants != '']
      
      variants = variants[!duplicated(positions)]
      positions = positions[!duplicated(positions)]
      
      asv_seq = strsplit(asv_seq, '')[[1]]
      
      for(position in positions){
        variant = variants[which(positions == position)]
        
        if(!grepl('(I|D)=', variant)){
          
          asv_seq[position] = variant
          
        }else if(grepl('I=', variant)){
          
          asv_seq[position] = paste0(asv_seq[position], gsub('I=', '', variant))
          
        }else if(grepl('D=', variant)){
          
          variant = gsub('D=', '', variant)
          
          asv_seq[position:(position + nchar(variant))] = NA
          
        }
        
      }
      
      asv_seq = asv_seq[!is.na(asv_seq)]
      
      asv_seq = paste(asv_seq, collapse = '')
      
      asv_seqs[[asv]] = asv_seq
      
      
      
    }else{
      
      asv_seqs[[asv]] = asv_seq
      
    }
    
    asv_table[asv_table$hapid == asv,][['haplength']] = nchar(asv_seqs[[asv]])
    
  }
  

  modified_gt = ampseq_gt_table
  
  duplicated_asvs = names(asv_seqs)[(duplicated(unlist(asv_seqs)))]
  
  for(duplicated_asv in duplicated_asvs){

    duplicated_asv_ids = names(asv_seqs)[asv_seqs %in% asv_seqs[duplicated_asv]]
    
    if(sum(asv_table[asv_table$hapid %in% duplicated_asv_ids,][['total_reads']]) == 0){
      
      for(hap_id_to_remove in duplicated_asv_ids){
        asv_seqs[[hap_id_to_remove]] = NULL
        asv_table = asv_table[!(asv_table$hapid %in% hap_id_to_remove),]
      }
      
    }else{
      
      duplicated_asv_ids_with_zero_reads = which(asv_table[asv_table$hapid %in% duplicated_asv_ids,][['total_reads']] == 0)
      duplicated_asv_ids_with_nozero_reads = which(asv_table[asv_table$hapid %in% duplicated_asv_ids,][['total_reads']] != 0)
      
      if(length(duplicated_asv_ids_with_nozero_reads) == 1){
        
        for(hap_id_to_remove in duplicated_asv_ids[duplicated_asv_ids_with_zero_reads]){
          asv_seqs[[hap_id_to_remove]] = NULL
          asv_table = asv_table[!(asv_table$hapid %in% hap_id_to_remove),]
        }
        
      }else{
        
        hap_id_to_keep = duplicated_asv_ids[
          which.max(asv_table[(asv_table$hapid %in% duplicated_asv_ids),][['total_reads']])]
        
        hap_ids_to_remove = duplicated_asv_ids[duplicated_asv_ids != hap_id_to_keep]
        
        hap_id_to_remove = hap_ids_to_remove[1]
        
        for(hap_id_to_remove in hap_ids_to_remove){
          
          reads_in_id_to_remove = asv_table[(asv_table$hapid == hap_id_to_remove),][['total_reads']]
          
          
          if(reads_in_id_to_remove > 0){
            
            amplicon_of_id_to_remove = asv_table[(asv_table$hapid == hap_id_to_remove),][['Amplicon']]
            samples_in_id_to_remove = asv_table[(asv_table$hapid == hap_id_to_remove),][['total_samples']]
            cigar_masked_in_id_to_remove = asv_table[(asv_table$hapid == hap_id_to_remove),][['CIGAR_masked']]
            

            
            reads_in_id_to_keep = asv_table[(asv_table$hapid == hap_id_to_keep),][['total_reads']]
            samples_in_id_to_keep = asv_table[(asv_table$hapid == hap_id_to_keep),][['total_samples']]
            cigar_masked_in_id_to_keep = asv_table[(asv_table$hapid == hap_id_to_keep),][['CIGAR_masked']]
            
            modified_gt[grepl(paste0('^', cigar_masked_in_id_to_remove, ':'), modified_gt[,amplicon_of_id_to_remove]),
                        amplicon_of_id_to_remove
            ] = gsub(paste0('^', cigar_masked_in_id_to_remove, ':'),
                     paste0(cigar_masked_in_id_to_keep, ':'),
                     modified_gt[grepl(paste0('^', cigar_masked_in_id_to_remove, ':'),
                                       modified_gt[,amplicon_of_id_to_remove]),
                                 amplicon_of_id_to_remove])
            
            modified_gt[grepl(paste0('_', cigar_masked_in_id_to_remove, ':'), modified_gt[,amplicon_of_id_to_remove]),
                        amplicon_of_id_to_remove
            ] = gsub(paste0('_', cigar_masked_in_id_to_remove, ':'),
                     paste0('_',cigar_masked_in_id_to_keep, ':'),
                     modified_gt[grepl(paste0('_', cigar_masked_in_id_to_remove, ':'),
                                       modified_gt[,amplicon_of_id_to_remove]),
                                 amplicon_of_id_to_remove])
            
            
            asv_table[(asv_table$hapid == hap_id_to_keep),][['total_reads']] =
              reads_in_id_to_keep + reads_in_id_to_remove
            asv_table[(asv_table$hapid == hap_id_to_keep),][['total_samples']] =
              samples_in_id_to_keep + samples_in_id_to_remove
            
          }
          
          asv_seqs[[hap_id_to_remove]] = NULL
          asv_table = asv_table[!(asv_table$hapid %in% hap_id_to_remove),]
          
          
        }
        
      }
      
    }
    
    
  }
  

  

  
  asv_seqs = Biostrings::DNAStringSet(unlist(asv_seqs))
  
  asv_table = asv_table[, c(
    "hapid",
    "haplength",
    "total_reads",
    "total_samples",
    "strain",
    "refid_PvP01",
    "snv_dist_from_PvP01",
    "indel_dist_from_PvP01",
    "bimera",
    "Amplicon",
    "CIGAR",
    "CIGAR_masked" 
  )]
  

  ampseq_object = create_ampseq(gt = modified_gt,
                metadata = rGenome_object@metadata,
                markers = markers,
                asv_table = asv_table,
                asv_seqs = asv_seqs,
                asv_seqs_masked = asv_seqs
                )
  
  
  return(ampseq_object)
}



## run_vcftools ----
run_vcftools = function(vcf = NULL,
                        
                        bash_file = NULL, 
                        
                        out = NULL,
                        
                        keep_regexp = NULL, # regular expression pattern that identify samples
                        remove_regexp = NULL,
                        
                        keep = NULL,
                        remove = NULL,
                        
                        chr = NULL,
                        not_chr = NULL,
                        
                        bed = NULL,
                        exclude_bed = NULL,
                        
                        positions = NULL,
                        exclude_positions = NULL,
                        
                        keep_only_indels = FALSE,
                        remove_indels = FALSE,
                        
                        remove_filtered_all = FALSE,
                        
                        maf = NULL,
                        max_maf = NULL,
                        
                        non_ref_af = NULL,
                        max_non_ref_af = NULL,
                        non_ref_ac = NULL,
                        max_non_ref_ac = NULL,
                        
                        non_ref_af_any = NULL,
                        max_non_ref_af_any = NULL,
                        non_ref_ac_any = NULL,
                        max_non_ref_ac_any = NULL,
                        
                        mac = NULL,
                        max_mac = NULL,
                        
                        min_alleles = NULL,
                        max_alleles = NULL,
                        
                        # Output options
                        freq = FALSE,
                        counts = FALSE,
                        
                        depth = FALSE,
                        site_depth = FALSE,
                        site_mean_depth = FALSE,
                        geno_depth = FALSE,
                        
                        hap_r2 = FALSE,
                        geno_r2 = FALSE,
                        geno_chisq = FALSE,
                        hap_r2_positions = NULL,
                        geno_r2_positions = NULL,
                        ld_window = NULL,
                        ld_window_bp = NULL,
                        ld_window_min = NULL,
                        ld_window_bp_min = NULL,
                        min_r2 = NULL,
                        interchrom_hap_r2 = FALSE,
                        interchrom_geno_r2 = FALSE,
                        
                        TsTv = NULL,
                        TsTv_by_count = FALSE,
                        TsTv_by_qual = FALSE,
                        
                        site_pi = FALSE,
                        window_pi = NULL,
                        window_pi_step = NULL,
                        
                        weir_fst_pop = NULL,
                        fst_window_size = NULL,
                        fst_window_step = NULL,
                        
                        het = FALSE,
                        TajimaD = NULL,
                        relatedness = FALSE,
                        relatedness2 = FALSE,
                        
                        recode = FALSE,
                        recode_bcf = FALSE,
                        recode_INFO_all = FALSE
){
  
  print('Starting VCFTools')
  
  vcf_run_file = c('#!/bin/bash',
                   'source /broad/software/scripts/useuse',
                   'use .vcftools-0.1.14',
                   'use Tabix')
  
  vcf_arguments = 'vcftools'
  
  if(grepl('.vcf$', vcf)){
    vcf_arguments = paste0(vcf_arguments, ' --vcf ', vcf)
  }else if(grepl('.vcf.gz$', vcf)){
    vcf_arguments = paste0(vcf_arguments, ' --gzvcf ', vcf)
  }else if(grepl('.bcf$', vcf)){
    vcf_arguments = paste0(vcf_arguments, ' --bcf ', vcf)}
  
  # Filters
  
  if(!is.null(keep_regexp)){
    system(paste0("zgrep '^#[A-Z]' ",  vcf, " > ", "temp_samples.indv"))
    temp_samples = as.character(read.csv("temp_samples.indv", header = FALSE, sep = '\t'))
    samples = temp_samples[grepl(keep_regexp,temp_samples)]
    write.table(samples, 'samples.indv', sep = '\t', quote = FALSE, row.names = FALSE, col.names = F)
    system(paste0('rm ', "temp_samples.indv"))
    vcf_arguments = paste0(vcf_arguments, ' --keep samples.indv')
  }
  
  if(!is.null(remove_regexp)){
    system(paste0("zgrep '^#[A-Z]' ",  vcf, " > ", "temp_samples.indv"))
    temp_samples = as.character(read.csv("temp_samples.indv", header = FALSE, sep = '\t'))
    samples = temp_samples[grepl(keep_regexp,temp_samples)]
    write.table(samples, 'rsamples.indv', sep = '\t', quote = FALSE, row.names = FALSE, col.names = F)
    system(paste0('rm ', "temp_samples.indv"))
    vcf_arguments = paste0(vcf_arguments, ' --remove rsamples.indv')
  }
  
  if(!is.null(keep)){
    vcf_arguments = paste0(vcf_arguments, ' --keep ', keep)
  }
  
  if(!is.null(remove)){
    vcf_arguments = paste0(vcf_arguments, ' --remove ', remove)
  }
  
  if(!is.null(chr)){
    for(chromosome in chr){
      vcf_arguments = paste0(vcf_arguments, ' --chr ', chromosome)
    }
  }
  
  if(!is.null(not_chr)){
    for(chromosome in not_chr){
      vcf_arguments = paste0(vcf_arguments, ' --not-chr ', chromosome)
    }
  }
  
  if(!is.null(bed)){vcf_arguments = paste0(vcf_arguments, ' --bed ', bed)
  print('Adding bed argument')
  }
  if(!is.null(exclude_bed)){vcf_arguments = paste0(vcf_arguments, ' --exclude-bed ', exclude_bed)}
  
  if(!is.null(positions)){vcf_arguments = paste0(vcf_arguments, ' --positions ', positions)}
  if(!is.null(exclude_positions)){vcf_arguments = paste0(vcf_arguments, ' --exclude-positions ', exclude_positions)}
  
  if(keep_only_indels){vcf_arguments = paste0(vcf_arguments, ' --keep-only-indels')}
  if(remove_indels){vcf_arguments = paste0(vcf_arguments, ' --remove-indels')}
  
  if(remove_filtered_all){vcf_arguments = paste0(vcf_arguments, ' --remove-filtered-all')}
  
  if(!is.null(maf)){vcf_arguments = paste0(vcf_arguments, ' --maf ', maf)}
  if(!is.null(max_maf)){vcf_arguments = paste0(vcf_arguments, ' --max-maf ', max_maf)}
  
  if(!is.null(non_ref_af)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-af ', non_ref_af)}
  if(!is.null(max_non_ref_af)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-af ', max_non_ref_af)}
  if(!is.null(non_ref_ac)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-ac ', non_ref_ac)}
  if(!is.null(max_non_ref_ac)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-ac ', max_non_ref_ac)}
  
  if(!is.null(non_ref_af_any)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-af-any ', non_ref_af_any)}
  if(!is.null(max_non_ref_af_any)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-af-any ', max_non_ref_af_any)}
  if(!is.null(non_ref_ac_any)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-ac-any ', non_ref_ac_any)}
  if(!is.null(max_non_ref_ac_any)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-ac-any ', max_non_ref_ac_any)}
  if(!is.null(mac)){vcf_arguments = paste0(vcf_arguments, ' --mac ', mac)}
  if(!is.null(max_mac)){vcf_arguments = paste0(vcf_arguments, ' --max-mac ', max_mac)}
  
  if(!is.null(min_alleles)){vcf_arguments = paste0(vcf_arguments, ' --min-alleles ', min_alleles)}
  if(!is.null(max_alleles)){vcf_arguments = paste0(vcf_arguments, ' --max-alleles ', max_alleles)}
  
  # Output options
  
  if(freq){vcf_arguments = paste0(vcf_arguments, ' --freq')}
  if(counts){vcf_arguments = paste0(vcf_arguments, ' --counts')}
  
  if(depth){vcf_arguments = paste0(vcf_arguments, ' --depth')}
  if(site_depth){vcf_arguments = paste0(vcf_arguments, ' --site-depth')}
  if(site_mean_depth){vcf_arguments = paste0(vcf_arguments, ' --site-mean-depth')}
  if(geno_depth){vcf_arguments = paste0(vcf_arguments, ' --geno_depth')}
  
  if(hap_r2){vcf_arguments = paste0(vcf_arguments, ' --hap-r2')}
  if(geno_r2){vcf_arguments = paste0(vcf_arguments, ' --geno-r2')}
  if(geno_chisq){vcf_arguments = paste0(vcf_arguments, ' --geno-chisq')}
  if(!is.null(ld_window)){vcf_arguments = paste0(vcf_arguments, ' --ld-window ', ld_window)}
  if(!is.null(ld_window_bp)){vcf_arguments = paste0(vcf_arguments, ' --ld-window-bp ', ld_window_bp)}
  if(!is.null(ld_window_min)){vcf_arguments = paste0(vcf_arguments, ' --ld-window-min ', ld_window_min)}
  if(!is.null(ld_window_bp_min)){vcf_arguments = paste0(vcf_arguments, ' --ld-window-bp-min ', ld_window_bp_min)}
  if(!is.null(min_r2)){vcf_arguments = paste0(vcf_arguments, ' --min-r2 ', min_r2)}
  if(interchrom_hap_r2){vcf_arguments = paste0(vcf_arguments, ' --interchrom-hap-r2')}
  if(interchrom_geno_r2){vcf_arguments = paste0(vcf_arguments, ' --interchrom-geno-r2')}
  
  if(!is.null(TsTv)){vcf_arguments = paste0(vcf_arguments, ' --TsTv ', TsTv)}
  if(TsTv_by_count){vcf_arguments = paste0(vcf_arguments, ' --TsTv-by-count')}
  if(TsTv_by_qual){vcf_arguments = paste0(vcf_arguments, ' --TsTv-by-qual')}
  
  if(site_pi){vcf_arguments = paste0(vcf_arguments, ' --site-pi')}
  if(!is.null(window_pi)){vcf_arguments = paste0(vcf_arguments, ' --window-pi ', window_pi)}
  if(!is.null(window_pi_step)){vcf_arguments = paste0(vcf_arguments, ' --window-pi-step ', window_pi_step)}
  
  if(!is.null(weir_fst_pop)){
    for(pop in weir_fst_pop){
      vcf_arguments = paste0(vcf_arguments, ' --weir-fst-pop ', pop)
    }
  }
  
  if(!is.null(fst_window_size)){vcf_arguments = paste0(vcf_arguments, ' --fst-window-size ', fst_window_size)}
  if(!is.null(fst_window_step)){vcf_arguments = paste0(vcf_arguments, ' --fst-window-step ', fst_window_step)}
  
  if(het){vcf_arguments = paste0(vcf_arguments, ' --het')}
  if(!is.null(TajimaD)){vcf_arguments = paste0(vcf_arguments, ' --TajimaD ', TajimaD)}
  if(relatedness){vcf_arguments = paste0(vcf_arguments, ' --relatedness')}
  if(relatedness2){vcf_arguments = paste0(vcf_arguments, ' --relatedness2')}
  
  if(recode){vcf_arguments = paste0(vcf_arguments, ' --recode')}
  if(recode_bcf){vcf_arguments = paste0(vcf_arguments, ' --recode-bcf')}
  if(recode_INFO_all){vcf_arguments = paste0(vcf_arguments, ' --recode-INFO-all')}
  
  if(!is.null(out)){vcf_arguments = paste0(vcf_arguments, ' --out ', out)}
  
  vcf_run_file = c(vcf_run_file, vcf_arguments)
  
  print(vcf_arguments)
  
  write.table(vcf_run_file, bash_file, row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  system(paste0('chmod 777 ', bash_file))
  system(paste0('./', bash_file))
  system(paste0('rm ', bash_file))
  
}

## load_vcf----
load_vcf = function(vcf = NULL,
                    gzvcf = NULL,
                    na.rm = TRUE,
                    start = NULL,
                    end = NULL
){
  
  if(!is.null(vcf)){
    
    temp_tsv_name = gsub('vcf', 'tsv', vcf)
    system(paste0("grep -v '^##' ", vcf, " > ", temp_tsv_name))
  }
  
  if(!is.null(gzvcf)){
    system(paste0("zgrep -v '^##' ", vcf, " > ", temp_tsv_name))
  }
  
  col_names = strsplit(system(paste0("grep '^#' ", temp_tsv_name), intern = TRUE), '\t')[[1]]
  col_names[1] = 'CHROM'
  
  temp_tsv_name2 = gsub('.tsv$','_2.tsv',temp_tsv_name)
  
  if(is.null(start) | is.null(end)){
    system(paste0("grep -v '^#' ", temp_tsv_name, " > ", temp_tsv_name2))
  }else if(!is.null(start) & !is.null(end)){
    system(paste0("grep -v '^#' ", 
                  temp_tsv_name, 
                  " | sed -n ", 
                  as.character(as.integer(start)), 
                  ",", 
                  as.character(as.integer(end)), 
                  "p > ", 
                  temp_tsv_name2))
  }
  
  vcf = read.table(temp_tsv_name2, header = FALSE)
  names(vcf) = col_names
  
  if(na.rm){
    vcf %<>% filter(ALT != '.')
  }
  
  system(paste0('rm ', temp_tsv_name))
  system(paste0('rm ', temp_tsv_name2))
  
  return(vcf)
}


## rGenome S4class and vcf2rGenome----

## rGenome S4 class

setClass('rGenome', slots = c(
  gt = "ANY",
  loci_table = "ANY",
  metadata = "ANY"
))

## rGenome constructor
rGenome = function(gt = NULL,
                   loci_table = NULL,
                   metadata = NULL){
  obj = new('rGenome')
  obj@gt = gt
  obj@loci_table = loci_table
  obj@metadata = metadata
  
  return(obj)
}

## vcf2rGenome----

vcf2rGenome = function(vcf, n = 500, threshold = 5) {
  
  # Generate metadata
  metadata = data.frame(Sample_id = names(vcf)[-1:-9])
  rownames(metadata) = metadata[['Sample_id']]
  
  # generate loci_table
  loci_table = vcf[,1:9]
  rownames(loci_table) = paste(loci_table$CHROM, loci_table$POS, sep = '_')
  
  # generate a haplotype table (gt)
  
  gt = NULL
  for(w in 1:n){
    start = Sys.time()
    gt = rbind(gt, get_GTAD_matrix(vcf, w = w, n = n, threshold = threshold))
    end = Sys.time()
    print(w)
    print(end-start)
  }
  
  obj = rGenome(gt = gt, loci_table = loci_table, metadata = metadata)
  
  return(obj)
}

## rGenome2vcf----
rGenome2vcf = function(rGenome_object){
  
  loci_table = rGenome_object@loci_table
  
  vcf_names1 = c("CHROM", "POS", "ID",  "REF", "ALT",  "QUAL", "FILTER", "INFO", "FORMAT")
  
  loci_table$INFO = paste(loci_table$INFO, apply(loci_table[, !(colnames(loci_table) %in% vcf_names1)], 1, function(x){paste(paste(names(x), x, sep = '='), collapse = ';')}), sep = ';')
  
  loci_table$FORMAT = 'GT:AD'
  
  loci_table = loci_table[,vcf_names1]
  
  colnames(loci_table) = c("#CHROM", "POS", "ID",  "REF", "ALT",  "QUAL", "FILTER", "INFO", "FORMAT")
  
  gt = rGenome_object@gt
  
  ad = matrix(gsub('\\d+:', '', gt), nrow = nrow(gt), ncol = ncol(gt),
              dimnames = list(rownames(gt), colnames(gt)))
  
  ad = gsub('/', ',', ad)
  
  ad[is.na(ad)] = 0
  
  gt = matrix(gsub(':\\d+', '', gt), nrow = nrow(gt), ncol = ncol(gt),
              dimnames = list(rownames(gt), colnames(gt)))
  
  
  gt1 = gsub('/\\d+', '', gt)
  gt2 = gsub('\\d+/', '', gt)
  
  gt3 = matrix(paste(gt1, gt2 , sep = '/'), nrow = nrow(gt), ncol = ncol(gt),
               dimnames = list(rownames(gt), colnames(gt)))
  
  gt3 = gsub('NA', '.', gt3)
  
  final_gt = matrix(paste(gt3, ad , sep = ':'), nrow = nrow(gt), ncol = ncol(gt),
                    dimnames = list(rownames(gt), colnames(gt)))
  
  final_vcf = cbind(loci_table, final_gt)
  
  return(final_vcf)
  
}


get_GTAD_matrix = function(vcf, w = 1, n = 100, threshold = 5){
  start = Sys.time()
  s = round(seq(1,nrow(vcf)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  w_gt_table = vcf[low:high,-1:-9]
  
  gt_temp = gsub(':.+', '', as.matrix(w_gt_table))
  ad_temp = gsub(':.+', '', gsub('^\\d+(/|\\|)\\d+:', '', as.matrix(w_gt_table)))
  
  if(!is.null(dim(w_gt_table))){
    
    gt = t(sapply(1:nrow(gt_temp), function(variant) {
      sapply(1:ncol(gt_temp), function(sample){
        gt = strsplit(as.character(gt_temp[variant, sample]), '/|\\|')[[1]]
        ad = strsplit(as.character(ad_temp[variant, sample]), ',')[[1]]
        
        if(gt[1] != '.'){
          ad = ad[as.numeric(gt) + 1]
          gt_df = data.frame(gt = gt, ad = as.numeric(ad))
          gt_df = gt_df[order(gt_df$ad, decreasing = T),]
          gt_df = gt_df[gt_df$ad >= threshold,]
          gt = paste(gt_df$gt, gt_df$ad, sep = ':')
          gt = unique(gt)
          
          if(length(gt) > 0){
            gt = paste(gt, collapse = '/')
          }else{
            gt = NA
          }
          
        }else{
          gt = NA
        }
        
      })
    }))
    
  }else{
    
    gt = matrix(sapply(1:nrow(gt_temp), function(variant) {
      sapply(1:ncol(gt_temp), function(sample){
        gt = strsplit(as.character(gt_temp[variant, sample]), '/|\\|')[[1]]
        ad = strsplit(as.character(ad_temp[variant, sample]), ',')[[1]]
        
        if(gt[1] != '.'){
          ad = ad[as.numeric(gt) + 1]
          gt_df = data.frame(gt = gt, ad = as.numeric(ad))
          gt_df = gt_df[order(gt_df$ad, decreasing = T),]
          gt_df = gt_df[gt_df$ad >= threshold,]
          gt = paste(gt_df$gt, gt_df$ad, sep = ':')
          gt = unique(gt)
          
          if(length(gt) > 0){
            gt = paste(gt, collapse = '/')
          }else{
            gt = NA
          }
          
        }else{
          gt = NA
        }
        
      })
    }), ncol = 1)
    
  }
  
  
  
  end = Sys.time()
  end - start
  
  rownames(gt) = vcf[low:high,1:2] %>% mutate(Locus = paste(CHROM, POS, sep = '_')) %>% dplyr::select(Locus) %>% unlist
  colnames(gt) = colnames(vcf)[-1:-9]
  
  return(gt)
  
}

## get_ampseq_limit_of_detection----

get_ampseq_limit_of_detection = function(ampseq_object,
                                         read_depth_thresholds,
                                         AmpRate_Thresholds,
                                         sample_concentration,
                                         for_markers = FALSE,
                                         for_procedure = TRUE,
                                         strata,
                                         upper_limit_marker_model = 1000,
                                         nObs_marker_model = 1000,
                                         digits = 6,
                                         x_lable = 'DNA parasite concentration by qPCR',
                                         y_lable = 'Fraction of amplified loci',
                                         color_lable = 'Amp. Rate\nThreshold',
                                         model = 'glm' # 'lm'
){
  
  registerDoMC(detectCores())
  
  coverage_by_marker_by_sample = get_ReadDepth_coverage(ampseq_object, variable = strata, plot = FALSE)[[1]]
  
  names(coverage_by_marker_by_sample) = c('Sample_id', 'Strata', 'Markers', 'Read_depth')
  
  coverage_by_marker_by_sample = left_join(coverage_by_marker_by_sample, ampseq_object@metadata[,c('Sample_id', sample_concentration)], by = 'Sample_id')
  
  coverage_by_marker_by_sample %<>% dplyr::rename('Concentration' = all_of(sample_concentration))
  
  coverage_by_marker_by_sample %<>%
    mutate(Read_depth = case_when(
      is.na(Read_depth) ~ 0,
      !is.na(Read_depth) ~ Read_depth
    ))
  
  
  # Performance across markers
  
  if(for_markers){
    
    markers_performance = NULL
    
    for(threshold in read_depth_thresholds){
      
      markers_performance = rbind(markers_performance,
                                  coverage_by_marker_by_sample %>%
                                    mutate(RD_Threshold = threshold,
                                           AmpSuccess = as.integer(Read_depth >= threshold)))
      
      # for(marker in unique(coverage_by_marker_by_sample$Markers)){
      #   
      #   markers_performance = rbind(markers_performance,
      #                               coverage_by_marker_by_sample %>%
      #                                 filter(Markers == marker) %>%
      #                                 mutate(RD_Threshold = threshold,
      #                                        median_RD = median(Read_depth),
      #                                        mean_RD = mean(Read_depth),
      #                                        Total_RD = sum(Read_depth),
      #                                        AmpRate = sum(Read_depth >= threshold)/n(),
      #                                        AmpSuccess = as.integer(Read_depth >= threshold)))
      #   
      # }
      
    }
    
    # Defining logistic model accross markers
    
    # Pop = 'Guyana'
    # RD_Thres = 10
    # marker = 'PVP01_0109700_1'
    
    marker_predicted_model = 
      foreach(Pop = unique(markers_performance$Strata),
              .combine = rbind) %do% {
                foreach(RD_Thres = read_depth_thresholds,
                        .combine = rbind)%do%{
                          foreach(marker = ampseq_object@markers$amplicon,
                                  .combine = rbind)%dopar% {
                                    
                                    # Fit the model and get coefficients
                                    
                                    
                                    glm_model = glm(AmpSuccess ~ Concentration, family = quasibinomial, data = markers_performance %>% filter(RD_Threshold == RD_Thres, Strata == Pop, Markers == marker))
                                    
                                    if(is.null(upper_limit_marker_model)){
                                      
                                      upper_limit_marker_model = max(markers_performance$Concentration, na.rm = T)
                                      
                                    }
                                    
                                    if(is.null(nObs_marker_model)){
                                      
                                      nObs_marker_model = 100
                                      
                                    }
                                    
                                    newdata = data.frame(Concentration = seq(0, upper_limit_marker_model, length.out = nObs_marker_model))
                                    
                                    preddat = predict(glm_model, newdata = newdata,  se.fit=TRUE)
                                    
                                    coefs = data.frame(t(coefficients(glm_model)))
                                    
                                    data.frame(
                                      Marker = marker,
                                      Concentration = seq(0, upper_limit_marker_model, length.out = nObs_marker_model),
                                      AmpSuccess = exp(preddat$fit) / (1 + exp(preddat$fit)),
                                      lower = exp(preddat$fit - 1.96 * preddat$se.fit) / (1 + exp(preddat$fit - 1.96 * preddat$se.fit)),
                                      upper= exp(preddat$fit + 1.96 * preddat$se.fit) / (1 + exp(preddat$fit + 1.96 * preddat$se.fit)),
                                      RD_Threshold = RD_Thres,
                                      Strata = Pop,
                                      Coef_Intercept = coefs$X.Intercept.,
                                      Coef_Concentration = coefs$Concentration
                                      
                                    )
                                    
                                  }
                        }
              }
    
    # Calculating LOD and confidence intervals
    marker_LOD_Values = foreach(
      Pop = unique(markers_performance$Strata), 
      .combine = rbind
    ) %do% {
      foreach(
        RD_Thres = read_depth_thresholds,
        .combine = rbind
      )%do%{
        foreach(
          AmpRate = AmpRate_Thresholds,
          .combine = rbind
        )%do%{
          foreach(
            marker = ampseq_object@markers$amplicon,
            .combine = rbind
          )%dopar%{
            
            marker_summary = markers_performance %>% 
              filter(RD_Threshold == RD_Thres, Strata == Pop, Markers == marker) %>%
              summarise(AmpRate = sum(AmpSuccess)/n(),
                        Median = median(Read_depth),
                        Mean = mean(Read_depth),
                        Total_RD = sum(Read_depth),
                        q25 = quantile(Read_depth, .25),
                        q75 = quantile(Read_depth, .75),
                        sd = sd(Read_depth),
                        IQR = IQR(Read_depth))
            
            temp_marker_predicted_model = marker_predicted_model %>% 
              filter(Strata == Pop,
                     RD_Threshold == RD_Thres,
                     Marker == marker)
            
            Intercept = temp_marker_predicted_model %>% 
              filter(Strata == Pop,
                     RD_Threshold == RD_Thres,
                     Marker == marker) %>%
              select(Coef_Intercept) %>% unlist()
            
            Intercept = mean(Intercept, na.rm = TRUE)
            
            Concentration = temp_marker_predicted_model %>%
              filter(Strata == Pop,
                     RD_Threshold == RD_Thres,
                     Marker == marker) %>%
              select(Coef_Concentration) %>% unlist
            
            Concentration = mean(Concentration, na.rm = T)
            
            # Get LOD
            
            LOD = ((log(AmpRate / (1 - AmpRate)) - Intercept) / Concentration)
            
            
            
            if(LOD > 0.001 & LOD < 10000 & Concentration > 0){
              
              glm_model = glm(AmpSuccess ~ Concentration, family = quasibinomial, data = markers_performance %>% filter(RD_Threshold == RD_Thres, Strata == Pop, Markers == marker))
              
              # Get the lower limit of the LOD Value
              
              LowerLim_LOD = ceiling(LOD)
              
              predSup = predict(glm_model, newdata = data.frame('Concentration' = (LowerLim_LOD)), se.fit = TRUE)
              LowerLim_AmpRate = exp(predSup$fit + 1.96 * predSup$se.fit) / (1 + exp(predSup$fit + 1.96 * predSup$se.fit))
              
              if(min(temp_marker_predicted_model$upper, na.rm = T) < AmpRate){
                
                for(precision in 10^(-(0:digits))){
                  
                  while(LowerLim_AmpRate >= AmpRate & LowerLim_LOD > 0){
                    LowerLim_LOD = LowerLim_LOD - precision
                    predSup = predict(glm_model, newdata = data.frame('Concentration' = (LowerLim_LOD)), se.fit = TRUE)
                    LowerLim_AmpRate = exp(predSup$fit + 1.96 * predSup$se.fit) / (1 + exp(predSup$fit + 1.96 * predSup$se.fit))
                    
                  }
                  
                  LowerLim_LOD = LowerLim_LOD + precision
                  
                  predSup = predict(glm_model, newdata = data.frame('Concentration' = (LowerLim_LOD)), se.fit = TRUE)
                  LowerLim_AmpRate = exp(predSup$fit + 1.96 * predSup$se.fit) / (1 + exp(predSup$fit + 1.96 * predSup$se.fit))
                  
                }
                
              }else{
                
                LowerLim_LOD = NA
                
              }
              
              # Get the upper limit of the LOD Value
              
              UpperLim_LOD = ceiling(LOD)
              
              predInf = predict(glm_model, newdata = data.frame('Concentration' = (UpperLim_LOD)), se.fit = TRUE)
              UpperLim_AmpRate = exp(predInf$fit - 1.96 * predInf$se.fit) / (1 + exp(predInf$fit - 1.96 * predInf$se.fit))
              
              
              if(max(temp_marker_predicted_model$lower, na.rm = T) > AmpRate){
                
                for(precision in 10^(-(0:digits))){
                  
                  while(UpperLim_AmpRate <= AmpRate){
                    UpperLim_LOD = UpperLim_LOD + precision
                    predInf = predict(glm_model, newdata = data.frame('Concentration' = (UpperLim_LOD)), se.fit = TRUE)
                    UpperLim_AmpRate = exp(predInf$fit - 1.96 * predInf$se.fit) / (1 + exp(predInf$fit - 1.96 * predInf$se.fit))
                    
                  }
                  
                  UpperLim_LOD = UpperLim_LOD - precision
                  
                  predInf = predict(glm_model, newdata = data.frame('Concentration' = (UpperLim_LOD)), se.fit = TRUE)
                  UpperLim_AmpRate = exp(predInf$fit - 1.96 * predInf$se.fit) / (1 + exp(predInf$fit - 1.96 * predInf$se.fit))
                  
                }
                
              }else{
                
                UpperLim_LOD = NA
                
              }
              
              
            }else{
              
              LOD = NA
              LowerLim_LOD = NA
              UpperLim_LOD = NA
              
            }
            
            # Store LOD values
            
            data.frame(
              Marker = marker,
              Strata = Pop,
              RD_Threshold = RD_Thres,
              AmpRate_Threshold = AmpRate,
              LOD = LOD,
              LowerLim_LOD = LowerLim_LOD,
              UpperLim_LOD = UpperLim_LOD,
              AmpRate = marker_summary$AmpRate,
              Median = marker_summary$Median,
              Mean = marker_summary$Mean,
              Total_RD = marker_summary$Total_RD,
              q25 = marker_summary$q25,
              q75 = marker_summary$q75,
              sd = marker_summary$sd,
              IQR = marker_summary$IQR
            )
            
          }
        }}}
    
    # Plotting LOD across markers
    
    plot_LOD_by_marker = marker_LOD_Values %>%
      ggplot(aes(x = factor(RD_Threshold), y = log2(LOD))) +
      geom_jitter(aes(color = factor(RD_Threshold), size = Median), width = .1, alpha = .25) +
      geom_line(aes(group = Marker), linetype = 2, alpha = 0.025)+
      geom_hline(yintercept = log2(c(1, 10, 100, 500, 1000))) +
      scale_y_continuous(breaks = log2(c(1, 10, 100, 500, 1000)), labels = c(1, 10, 100, 500, 1000))+
      theme_minimal()+
      theme(legend.position = 'none')+
      facet_grid(AmpRate_Threshold ~ Strata)
    
  }
  
  
  # Performance across samples
  
  if(for_procedure){
    
    sample_performance = NULL
    
    for(threshold in read_depth_thresholds){
      
      sample_performance = rbind(sample_performance,
                                 coverage_by_marker_by_sample %>%
                                   dplyr::summarise(RD_Threshold = threshold,
                                                    Strata = unique(.[.[['Sample_id']] == Sample_id,][['Strata']]),
                                                    Concentration = unique(Concentration),
                                                    median_RD = median(Read_depth),
                                                    AmpRate = sum(Read_depth >= threshold) / nrow(ampseq_object@markers),
                                                    .by = Sample_id))
    }
    
    sample_performance %<>% 
      filter(!is.na(Concentration))%>%
      mutate(
        log_Concentration = log2(Concentration),
        log_OR_AmplRate = log2((AmpRate + 1 / nrow(ampseq_object@markers)) / (1 - (AmpRate + 1 / nrow(ampseq_object@markers)))),
        
        log_OR_AmplRate = case_when(
          is.na(log_OR_AmplRate) ~ nrow(ampseq_object@markers),
          is.infinite(log_OR_AmplRate) ~ nrow(ampseq_object@markers),
          .default = log_OR_AmplRate
        )
      )
    
    predicted_model = NULL
    LOD_Values = NULL
    
    if(model == 'lm'){
      
      for(Pop in unique(sample_performance$Strata)){
        for(RD_Thres in read_depth_thresholds){
          
          # Fit the model and get coeficients
          
          lm_model = lm(log_OR_AmplRate~log_Concentration, data = sample_performance %>% 
                          filter(RD_Threshold == RD_Thres, 
                                 Strata == Pop#, log_OR_AmplRate >= nrow(ampseq_object@markers) + 1
                          )
          )
          
          newdata = data.frame(log_Concentration = seq(1, max(sample_performance$log_Concentration), 0.1))
          
          preddat = predict(lm_model, newdata = newdata,  se.fit=TRUE)
          
          coefs = data.frame(t(coefficients(lm_model)))
          
          predicted_model = rbind(predicted_model, 
                                  data.frame(log_Concentration = seq(1, max(sample_performance$log_Concentration), 0.1),
                                             AmpRate = exp(preddat$fit) / (1 + exp(preddat$fit)),
                                             lower = exp(preddat$fit - 1.96 * preddat$se.fit) / (1 + exp(preddat$fit - 1.96 * preddat$se.fit)),
                                             upper= exp(preddat$fit + 1.96 * preddat$se.fit) / (1 + exp(preddat$fit + 1.96 * preddat$se.fit)),
                                             RD_Threshold = RD_Thres,
                                             Strata = Pop
                                  ))
          
          # Get LOD values and cofident intervals 
          
          for(AmpRate in AmpRate_Thresholds){
            
            # Get LOD
            
            LOD = 2^((log(AmpRate / (1 - AmpRate)) - coefs$X.Intercept.)/coefs$log_Concentration)
            
            
            if(LOD > 0.001 & LOD < 10000 & coefs$log_Concentration > 0.1){
              
              # Get the lower limit of the LOD Value
              
              LowerLim_LOD = ceiling(LOD)
              
              predSup = predict(lm_model, newdata = data.frame(log_Concentration = log2(LowerLim_LOD)), se.fit = TRUE)
              LowerLim_AmpRate = exp(predSup$fit + 1.96 * predSup$se.fit) / (1 + exp(predSup$fit + 1.96 * predSup$se.fit))
              
              for(precision in 10^(-(0:digits))){
                
                while(LowerLim_AmpRate >= AmpRate & LowerLim_LOD > 0){
                  LowerLim_LOD = LowerLim_LOD - precision
                  predSup = predict(lm_model, newdata = data.frame(log_Concentration = log2(LowerLim_LOD)), se.fit = TRUE)
                  LowerLim_AmpRate = exp(predSup$fit + 1.96 * predSup$se.fit) / (1 + exp(predSup$fit + 1.96 * predSup$se.fit))
                  
                }
                
                LowerLim_LOD = LowerLim_LOD + precision
                
                predSup = predict(lm_model, newdata = data.frame(log_Concentration = log2(LowerLim_LOD)), se.fit = TRUE)
                LowerLim_AmpRate = exp(predSup$fit + 1.96 * predSup$se.fit) / (1 + exp(predSup$fit + 1.96 * predSup$se.fit))
                
              }
              
              # Get the upper limit of the LOD Value
              
              UpperLim_LOD = ceiling(LOD)
              
              predInf = predict(lm_model, newdata = data.frame(log_Concentration = log2(UpperLim_LOD)), se.fit = TRUE)
              UpperLim_AmpRate = exp(predInf$fit - 1.96 * predInf$se.fit) / (1 + exp(predInf$fit - 1.96 * predInf$se.fit))
              
              
              for(precision in 10^(-(0:digits))){
                
                while(UpperLim_AmpRate <= AmpRate){
                  UpperLim_LOD = UpperLim_LOD + precision
                  predInf = predict(lm_model, newdata = data.frame(log_Concentration = log2(UpperLim_LOD)), se.fit = TRUE)
                  UpperLim_AmpRate = exp(predInf$fit - 1.96 * predInf$se.fit) / (1 + exp(predInf$fit - 1.96 * predInf$se.fit))
                  
                }
                
                UpperLim_LOD = UpperLim_LOD - precision
                
                predInf = predict(lm_model, newdata = data.frame(log_Concentration = log2(UpperLim_LOD)), se.fit = TRUE)
                UpperLim_AmpRate = exp(predInf$fit - 1.96 * predInf$se.fit) / (1 + exp(predInf$fit - 1.96 * predInf$se.fit))
                
              }
              
              # Store LOD values
              
              LOD_Values = rbind(LOD_Values,
                                 data.frame(
                                   Strata =  Pop,
                                   RD_Threshold = RD_Thres,
                                   AmpRate_Threshold = AmpRate,
                                   LOD = LOD,
                                   LowerLim_LOD = LowerLim_LOD,
                                   UpperLim_LOD = UpperLim_LOD
                                 )
              )
              
            }else{
              
              LOD_Values = rbind(LOD_Values,
                                 data.frame(
                                   Strata =  Pop,
                                   RD_Threshold = RD_Thres,
                                   AmpRate_Threshold = AmpRate,
                                   LOD = NA,
                                   LowerLim_LOD = NA,
                                   UpperLim_LOD = NA
                                 )
              )
              
            }
            
          }
          
          
        }
      }
      
    }else if(model == 'glm'){
      
      for(Pop in unique(sample_performance$Strata)){
        for(RD_Thres in read_depth_thresholds){
          
          for(AmpRate_Threshold in AmpRate_Thresholds){
            
            #stop()
            # Fit the model and get coeficients
            
            glm_model = glm(AmpRate ~ Concentration, 
                            family = 'quasibinomial', 
                            data = sample_performance %>% 
                              filter(RD_Threshold == RD_Thres, Strata == Pop) %>% 
                              mutate(AmpRate = as.integer(AmpRate >= AmpRate_Threshold))
            )
            
            
            newdata = data.frame(Concentration = seq(1, max(sample_performance$Concentration), 0.1))
            
            preddat = predict(glm_model, newdata = newdata,  se.fit=TRUE)
            
            coefs = data.frame(t(coefficients(glm_model)))
            
            predicted_model = 
              rbind(predicted_model,
                    data.frame(
                      Concentration = seq(1, max(sample_performance$Concentration), 0.1),
                      AmpRate = exp(preddat$fit) / (1 + exp(preddat$fit)),
                      lower = exp(preddat$fit - 1.96 * preddat$se.fit) / (1 + exp(preddat$fit - 1.96 * preddat$se.fit)),
                      upper= exp(preddat$fit + 1.96 * preddat$se.fit) / (1 + exp(preddat$fit + 1.96 * preddat$se.fit)),
                      RD_Threshold = RD_Thres,
                      Strata = Pop,
                      Coef_Intercept = coefs$X.Intercept.,
                      Coef_Concentration = coefs$Concentration))
            
            # Get LOD values and confident intervals 
            
            # Get LOD
            
            LOD = ((log(AmpRate_Threshold / (1 - AmpRate_Threshold)) - coefs$X.Intercept.) / coefs$Concentration)
            
            
            if(LOD > 0.001 & LOD < 10000 & coefs$Concentration > 0){
              
              # Get the lower limit of the LOD Value
              
              LowerLim_LOD = ceiling(LOD)
              
              predSup = predict(glm_model, newdata = data.frame(Concentration = LowerLim_LOD), se.fit = TRUE)
              LowerLim_AmpRate = exp(predSup$fit + 1.96 * predSup$se.fit) / (1 + exp(predSup$fit + 1.96 * predSup$se.fit))
              
              for(precision in 10^(-(0:digits))){
                
                while(LowerLim_AmpRate >= AmpRate_Threshold & LowerLim_LOD > 0){
                  LowerLim_LOD = LowerLim_LOD - precision
                  predSup = predict(glm_model, newdata = data.frame(Concentration = LowerLim_LOD), se.fit = TRUE)
                  LowerLim_AmpRate = exp(predSup$fit + 1.96 * predSup$se.fit) / (1 + exp(predSup$fit + 1.96 * predSup$se.fit))
                  
                }
                
                LowerLim_LOD = LowerLim_LOD + precision
                
                predSup = predict(glm_model, newdata = data.frame(Concentration = LowerLim_LOD), se.fit = TRUE)
                LowerLim_AmpRate = exp(predSup$fit + 1.96 * predSup$se.fit) / (1 + exp(predSup$fit + 1.96 * predSup$se.fit))
                
              }
              
              # Get the upper limit of the LOD Value
              
              UpperLim_LOD = floor(LOD)
              
              predInf = predict(glm_model, newdata = data.frame(Concentration = UpperLim_LOD), se.fit = TRUE)
              UpperLim_AmpRate = exp(predInf$fit - 1.96 * predInf$se.fit) / (1 + exp(predInf$fit - 1.96 * predInf$se.fit))
              
              
              for(precision in 10^(-(0:digits))){
                
                while(UpperLim_AmpRate <= AmpRate_Threshold){
                  UpperLim_LOD = UpperLim_LOD + precision
                  predInf = predict(glm_model, newdata = data.frame(Concentration = UpperLim_LOD), se.fit = TRUE)
                  UpperLim_AmpRate = exp(predInf$fit - 1.96 * predInf$se.fit) / (1 + exp(predInf$fit - 1.96 * predInf$se.fit))
                  
                }
                
                UpperLim_LOD = UpperLim_LOD - precision
                
                predInf = predict(glm_model, newdata = data.frame(Concentration = UpperLim_LOD), se.fit = TRUE)
                UpperLim_AmpRate = exp(predInf$fit - 1.96 * predInf$se.fit) / (1 + exp(predInf$fit - 1.96 * predInf$se.fit))
                
              }
              
              # Store LOD values
              
              LOD_Values = rbind(LOD_Values,
                                 data.frame(
                                   Strata =  Pop,
                                   RD_Threshold = RD_Thres,
                                   AmpRate_Threshold = AmpRate_Threshold,
                                   LOD = LOD,
                                   LowerLim_LOD = LowerLim_LOD,
                                   UpperLim_LOD = UpperLim_LOD
                                 ))
              
            }else{
              
              LOD_Values = rbind(LOD_Values,
                                 data.frame(
                                   Strata =  Pop,
                                   RD_Threshold = RD_Thres,
                                   AmpRate_Threshold = AmpRate_Threshold,
                                   LOD = NA,
                                   LowerLim_LOD = NA,
                                   UpperLim_LOD = NA
                                 )
              )
              
            }
            
          }
          
          
        }
      }
      
    }
    
    
    
    # Plot LOD
    
    if(model == 'lm'){
      
      plot_LOD = sample_performance %>%
        filter(RD_Threshold %in% read_depth_thresholds) %>%
        ggplot(aes(x = log2(Concentration), y = AmpRate)) +
        geom_point() +
        geom_line(data = predicted_model, aes(x = log_Concentration, y = AmpRate)) +
        geom_ribbon(data = predicted_model, aes(x = log_Concentration, ymin = lower, ymax = upper), 
                    linetype = 2, color = 'firebrick3', fill = 'firebrick3', alpha = .3) +
        geom_vline(data = LOD_Values, aes(xintercept = log2(LOD), color = as.character(AmpRate_Threshold)))+
        facet_grid(Strata ~ RD_Threshold) +
        theme_bw() +
        scale_x_continuous(breaks = log2(c(1, 10, 100, 1000)), labels = c(1, 10, 100, 1000)) +
        labs(x = x_lable,
             y = y_lable,
             color = color_lable)
      
    }else if(model == 'glm'){
      
      plot_LOD = sample_performance %>%
        filter(RD_Threshold %in% read_depth_thresholds) %>%
        ggplot(aes(x = log2(Concentration), y = AmpRate)) +
        geom_point() +
        geom_line(data = predicted_model, aes(x = log2(Concentration), y = AmpRate)) +
        geom_ribbon(data = predicted_model, aes(x = log2(Concentration), ymin = lower, ymax = upper), 
                    linetype = 2, color = 'firebrick3', fill = 'firebrick3', alpha = .3) +
        geom_vline(data = LOD_Values, aes(xintercept = log2(LOD), color = as.character(AmpRate_Threshold)))+
        facet_grid(Strata ~ RD_Threshold) +
        theme_bw() +
        scale_x_continuous(breaks = log2(c(1, 10, 100, 1000)), labels = c(1, 10, 100, 1000)) +
        labs(x = x_lable,
             y = y_lable,
             color = color_lable)
      
    }
    
  }
  
  if(for_markers & !for_procedure){
    
    return(list(
      markers_performance = markers_performance,
      marker_LOD_Values = marker_LOD_Values,
      marker_predicted_model = marker_predicted_model,
      plot_LOD_by_marker = plot_LOD_by_marker))
    
  }else if(!for_markers & for_procedure){
    
    return(list(
      sample_performance = sample_performance,
      predicted_model = predicted_model,
      LOD_Values = LOD_Values,
      plot_LOD = plot_LOD))
    
  }else if(for_markers & for_procedure){
    
    return(list(
      markers_performance = markers_performance,
      marker_LOD_Values = marker_LOD_Values,
      marker_predicted_model = marker_predicted_model,
      plot_LOD_by_marker = plot_LOD_by_marker,
      sample_performance = sample_performance,
      predicted_model = predicted_model,
      LOD_Values = LOD_Values,
      plot_LOD = plot_LOD))
    
  }else{
    
    return(coverage_by_marker_by_sample)
    
  }
  
}


## detect_erroneuns_cigars----
detect_erroneuns_cigars = function(ampseq_object){
  
  erroneuns_cigars = NULL
  
  for(mhap in colnames(ampseq_object@gt)){
    
    mhap_gt = ampseq_object@gt[,mhap]
    
    nclones = sapply(str_extract_all(mhap_gt, '_'), function(x){
      if(length(x) == 0){ 0 }else if(sum(is.na(x)) > 0){ 0 }else{ length(x) }
    }, simplify = T) + 1
    
    max_nalleles = max(nclones)
    
    for(clone in 1:max_nalleles){
      
      duplicated_positions = sapply(str_extract_all(gsub(':\\d+', '',gsub('_.+', '',mhap_gt)), '\\d+'),
                                    function(x){
                                      sum(duplicated(x))
                                    }, simplify = T)
      
      if(sum(duplicated_positions) > 0){
        
        erroneuns_cigars = rbind(erroneuns_cigars, data.frame(Amplicon = mhap,
                                                              Sample_id = names(mhap_gt[duplicated_positions >= 1]),
                                                              Genotype = mhap_gt[duplicated_positions >= 1]))
        
        mhap_gt[nclones == clone] = NA
        
      }
      
    } 
    
  }
  
  if(is.null(erroneuns_cigars)){
    print('There are no erroneous cigars in gt matrix')
  }else{
    
    cat('Erroneous cigars are:\n')
    
    return(erroneuns_cigars)
    
  }
  
}




are_alleles_in_asv_table = function(obj, allele_table){
  
  missalleles = NULL
  
  for(mhap in unique(allele_table$MHap)){
    
    mhap_alleles_1 = unique(allele_table[allele_table$MHap == mhap,][['Allele']])
    
    mhap_alleles_2 = obj@asv_table[obj@asv_table$Amplicon == mhap,][['CIGAR_masked']]
    
    if(sum(!(mhap_alleles_1 %in% mhap_alleles_2)) > 0){
      
      missalleles = rbind(missalleles,
                          data.frame(
                            Mhap = mhap,
                            Allele = mhap_alleles_1[!(mhap_alleles_1 %in% mhap_alleles_2)]
                          ))
      
    }
    
  }
  
  if(is.null(missalleles)){
    print('All alleles are in ASV table')
  }else{
    return(missalleles)
  }
  
}


## remove_filtered_asvs----
setGeneric("remove_filtered_asvs", function(obj = NULL
) standardGeneric("remove_filtered_asvs"))

setMethod("remove_filtered_asvs", signature(obj = "ampseq"),
          
          function(obj = NULL){
            
            gt = obj@gt
            
            mhaps = obj@markers
            
            asv_table = obj@asv_table
            
            
            alt = sapply(colnames(gt), function(mhap){
              alt = unique(unlist(strsplit(gsub(':\\d+', '', gt[,mhap]), '_')))
              
              alt = paste(alt[!is.na(alt)], collapse = ',')
            })
            
            for(mhap in names(alt)){
              
              # if(mhap == 'ApiAP2_1'){
              #   stop('error')
              # }
              
              alleles_in_asv_table = asv_table[asv_table[['Amplicon']] == mhap, ][['CIGAR_masked']]
              alleles_in_gt = unlist(strsplit(alt[[mhap]], ','))
              
              if(sum(!(alleles_in_asv_table %in% alleles_in_gt)) > 0){
                
                missing_alleles = alleles_in_asv_table[!(alleles_in_asv_table %in% alleles_in_gt)]
                
                asv_table[asv_table[['Amplicon']] == mhap &
                            asv_table[['CIGAR_masked']] %in% 
                            missing_alleles &
                            !is.na(asv_table[['CIGAR_masked']]), ][['CIGAR_masked']] = NA
                
              }
              
            }
              
            obj@asv_table = asv_table
            
            return(obj)
            
          }
)



## cigar_strings2fasta ----

cigar_strings2fasta = function(obj, 
                               ref_fasta, 
                               cigar_string_col = 'CIGAR', 
                               amplicon_col = 'Amplicon', 
                               format = 'data.frame' # 'list', 'DNAStringSet'
                               ){
  
  if(class(ref_fasta) == "DNAStringSet"){
    
    ref_fasta_seqs = ref_fasta
    
  }else{
    
    ref_fasta_seqs = Biostrings::readDNAStringSet(ref_fasta)
    
  }
  
  if(sum(grepl('(/|-|:)', names(ref_fasta_seqs))) > 0){
    
    names(ref_fasta_seqs) = gsub('(/|-|:)', '_', names(ref_fasta_seqs))
    print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the cigar table.')
    
  }
  
  if(format %in% c('data.frame', 'DNAStringSet')){
    
    seq_format = 'compact'
    
  }else if(format == 'list'){
    
    seq_format = 'split'
    
  }
  
  cigar_string2fasta = function(cigar_string, ref_fasta_seq, seq_format = 'compact'){
    
    if(is.na(cigar_string)){
      
      seq_of_interest = NA
      
    }else if(cigar_string == '.'){ # Sequence of interest is equal to the reference
      
      if(seq_format == 'split'){
        
        seq_of_interest = unlist(strsplit(ref_fasta_seq, ''))
        
      }else if(seq_format == 'compact'){
        
        seq_of_interest = ref_fasta_seq
        
      }
      
      
    }else{ # Sequence of interest is different to the reference
      
      positions = unlist(str_extract_all(cigar_string, '\\d+')) # collect all polymorphic positions from the cigar string
      
      if(sum(duplicated(positions)) > 0){
        
        print(cigar_string)
        stop('Wrong cigar string pattern, a position is duplicated')
      }
      
      variants = unlist(strsplit(cigar_string, '\\d+'))[-1] # collect all variants from cigar string
      
      seq_of_interest = unlist(strsplit(ref_fasta_seq, ''))
      
      for(position in positions){
        
        if(!grepl('^D', variants[which(positions == position)])){ # If polymorphism is a SNV or an Insertion
          
          if(position == '0'){ # Insertion at position 0
            
            seq_of_interest[1] = paste0(gsub('^I=', '', variants[which(positions == position)]),
                                        seq_of_interest[1])
            
          }else if(position == '1' & '0' %in% positions){ # SNV at position 1 after an Insertion at position 0
            
            substr(seq_of_interest[1], 
                   start = nchar(seq_of_interest[1]),
                   stop = nchar(seq_of_interest[1])) = gsub('^I=', '', variants[which(positions == position)])
            
          }else{ # SNVs or Insertions at other positions
            
            seq_of_interest[as.integer(position)] = gsub('^I=', '', variants[which(positions == position)])
            
          }
          
        }else{ # If polymorphism is a Deletion
          
          if(position == '1' & '0' %in% positions){ # If deletion is at position 1 after an insertion at position 0
            
            #stop("Del")
            seq_of_interest[1] = gsub('.$', '', seq_of_interest[1])
            
            if(nchar(variants[which(positions == position)]) - 3 > 0){
              deletion_starts = as.integer(position) + 1
              deletion_ends = as.integer(position) + nchar(variants[which(positions == position)]) - 3 # 2 for the "D" and "=" and 1 because is already included at the start
              seq_of_interest[deletion_starts:deletion_ends] = ''
            }
            
            
          }else{ # Deletions at other positions
            
            deletion_starts = as.integer(position)
            deletion_ends = as.integer(position) + nchar(variants[which(positions == position)]) - 3 # 2 for the "D" and "=" and 1 because is already included at the start
            seq_of_interest[deletion_starts:deletion_ends] = ''
            
          }
          
        }
        
      }
      
      if(seq_format == 'compact'){
        
        seq_of_interest = paste(seq_of_interest, collapse = '')
        
      }
      
    }
    
    return(seq_of_interest)
    
  }
  
  fasta_seqs = NULL
  
  for(i in 1:nrow(obj)){
    
    cigar_string = obj[[cigar_string_col]][i]
    Amplicon = obj[[amplicon_col]][i]
    ref_fasta_seq = as.character(ref_fasta_seqs[[Amplicon]])
    
    fasta_seq = cigar_string2fasta(cigar_string, ref_fasta_seq, seq_format)
    
    if(length(fasta_seq) == 0){
      stop(paste0('CIGAR string ', cigar_string, ' in Amplicon ', Amplicon, ' do not produced a Fasta sequence'))
    }
    
    if(format == 'data.frame'){
      
      fasta_seqs =
        rbind(fasta_seqs, data.frame(Amplicon = Amplicon, CIGAR = cigar_string, fasta = fasta_seq))
      
    }else if(format == 'list'){
      
      fasta_seqs[[Amplicon]] = fasta_seq
      
    }else if(format == 'DNAStringSet'){
      
      fasta_seqs = DNAStringSet(c(as.character(fasta_seqs), as.character(fasta_seq)))
      
    }
    
  }
  
  return(fasta_seqs)
  
}


# pairwise_euclidean----
pairwise_euclidean = function(obj = NULL, parallel = TRUE, w = 1, n = 100, alpha = 0.05, method = 'exact', pairs = NULL){
  library(parallel)
  library(doMC)
  library(svMisc)
  
  
  if(class(obj) == 'loci'){
    gt  = obj@loci_table
    
  }else if(class(obj) == 'ampseq'){
    loci_object = ampseq2loci(ampseq_object)
    gt  = loci_object@loci_table
  }else if(class(obj) == 'rGenome'){
    gt = t(obj@gt)
    gt = gsub('/', '_', gt)
    gt[grepl('\\.:\\.', gt)] = NA
  }
  
  gt = gsub(':\\d+', '', gt)
  
  if(is.null(pairs)){
    pairs = as.data.frame(t(combn(rownames(gt), 2)))  
  }
  
  s = round(seq(1,nrow(pairs)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  pairs = pairs[low:high,]
  
  get_euclidean_dist = function(x, y){
    
    dist_x_y = x != y
    
    estimate = binconf(sum(dist_x_y, na.rm = T),
              sum(!is.na(dist_x_y)),
              alpha = alpha,
              method = method, 
              include.x = T, 
              include.n = T, 
              return.df = T
      )
    
    return(estimate)
  }
  
  if(parallel){
    registerDoMC(detectCores())
    pairwise_df = foreach(pair = 1:nrow(pairs), .combine = 'rbind') %dopar% {
      
      Yi_id = pairs[pair, 1]
      Yj_id = pairs[pair, 2]
      
      
      Yi = split_clones(gt[Yi_id,], Yi_id)
      
      Yi = matrix(as.integer(Yi), 
                  nrow = nrow(Yi),
                  ncol = ncol(Yi),
                  dimnames = list(
                    rownames(Yi),
                    colnames(Yi)
                  ))
      
      Yj = split_clones(gt[Yj_id,], Yj_id)
      
      Yj = matrix(as.integer(Yj), 
                  nrow = nrow(Yj),
                  ncol = ncol(Yj),
                  dimnames = list(
                    rownames(Yj),
                    colnames(Yj)
                  ))
      
      estimate = NULL
      for(haplotype_i in 1:nrow(Yi)){
        for(haplotype_j in 1:nrow(Yj)){
          
          estimate = rbind(estimate,
                           get_euclidean_dist(Yi[haplotype_i,],
                             Yj[haplotype_j,]
                             ))
          
          
        }
      }
      
      if(sum(is.na(estimate[,'PointEst'])) > 0){
        
        estimate = estimate[1,]
        
      }else{
        
        estimate = estimate[which.min(estimate[,'PointEst']),]
        
      }
      
      names(estimate) = c('nDiff', 'nLoci', 'euDist', 'Lower_euDist', 'Upper_euDist')
      
      data.frame(Yi = Yi_id, Yj = Yj_id, estimate)
      
    }
    
  }else{
    
    pairwise_df = NULL
    
    for(pair in 1:nrow(pairs)){
      
      Yi_id = pairs[pair, 1]
      Yj_id = pairs[pair, 2]
      
      
      Yi = split_clones(genotype = gt[Yi_id,], genotype_name = Yi_id)
      
      Yi = matrix(as.integer(Yi), 
                  nrow = nrow(Yi),
                  ncol = ncol(Yi),
                  dimnames = list(
                    rownames(Yi),
                    colnames(Yi)
                  ))
      
      Yj = split_clones(gt[Yj_id,], Yj_id)
      
      Yj = matrix(as.integer(Yj), 
                  nrow = nrow(Yj),
                  ncol = ncol(Yj),
                  dimnames = list(
                    rownames(Yj),
                    colnames(Yj)
                  ))
      
      estimate = NULL
      for(haplotype_i in 1:nrow(Yi)){
        for(haplotype_j in 1:nrow(Yj)){
          
          estimate = rbind(estimate,
                           get_euclidean_dist(Yi[haplotype_i,],
                                              Yj[haplotype_j,]
                           ))
          
          
        }
      }
      
      if(sum(is.na(estimate[,'PointEst'])) > 0){
        
        estimate = estimate[1,]
        
      }else{
        
        estimate = estimate[which.min(estimate[,'PointEst']),]
        
      }
      
      names(estimate) = c('nDiff', 'nLoci', 'euDist', 'Lower_euDist', 'Upper_euDist')
      
      
      progress(round(100*pair/nrow(pairs)))
      
      pairwise_df = rbind(pairwise_df, data.frame(Yi = Yi_id, Yj = Yj_id, estimate))
      
    }
    
  }
  
  return(pairwise_df)
  
}


## pairwise_NucDiff----

pairwise_NucDiff = function(obj = NULL, 
                            parallel = TRUE, 
                            w = 1, 
                            n = 1, 
                            ref_amp_fasta,
                            pairs = NULL,
                            by_marker = FALSE
){
  library(parallel)
  library(doMC)
  library(svMisc)
  
  fasta_references = Biostrings::readDNAStringSet(ref_amp_fasta)
  
  if(sum(grepl('(/|-|:)', names(fasta_references))) > 0){
    
    names(fasta_references) = gsub('(/|-|:)', '_', names(fasta_references))
    print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the cigar table.')
    
  }
  
  gt = gsub(':\\d+', '', obj@gt)
  markers_table = obj@markers
  
  if(is.null(pairs)){
    pairs = as.data.frame(t(combn(rownames(gt), 2)))
  }
    
  s = round(seq(1,nrow(pairs)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  pairs = pairs[low:high,]
  
  get_NucDiff = function(x, y, markers_table, fasta_references, by_marker = by_marker){
    
    # Convert all cigar string for each marker in each sample into a fasta sequence
    
    x_cigars = data.frame(CIGAR = x, Amplicon = markers_table$amplicon)
    
    x_fastas = cigar_strings2fasta(obj = x_cigars, ref_fasta = fasta_references, format = 'list')
    
    y_cigars = data.frame(CIGAR = y, Amplicon = markers_table$amplicon)
    
    y_fastas = cigar_strings2fasta(obj = y_cigars, ref_fasta = fasta_references, format = 'list')
    
    # Remove missing fastas in both samples 
    
    x_missings = is.na(x_fastas)
    y_missings = is.na(y_fastas)
    
    x_cigars = x_cigars[!x_missings & !y_missings,]
    y_cigars = y_cigars[!x_missings & !y_missings,]
    
    x_cigar_INDELs = unlist(sapply(x_cigars[,1], 
                            function(marker){gsub('(\\.|(\\d+[ATGC]){0,})', '', marker)}))
    
    #x_cigar_INDELs = str_extract_all(x_cigar_INDELs, '\\d+(D|I)=[ATGC]+')
    
    names(x_cigar_INDELs) = x_cigars[,2]
    
    y_cigar_INDELs = unlist(sapply(y_cigars[,1], 
                                   function(marker){gsub('(\\.|(\\d+[ATGC]){0,})', '', marker)}))
    
    #y_cigar_INDELs = str_extract_all(y_cigar_INDELs, '\\d+(D|I)=[ATGC]+')
    
    names(y_cigar_INDELs) = y_cigars[,2]
    
    x_fastas = x_fastas[!x_missings & !y_missings]
    y_fastas = y_fastas[!x_missings & !y_missings]
    
    if(by_marker){
      
      estimate = NULL
      
      for(marker in names(x_fastas)){
        
        x_fastas_marker  = x_fastas[[marker]]
        y_fastas_marker  = y_fastas[[marker]]
        
        x_cigar_INDELs_marker = x_cigar_INDELs[[marker]]
        y_cigar_INDELs_marker = y_cigar_INDELs[[marker]]
        
        if(nchar(x_cigar_INDELs_marker) > 0 | nchar(y_cigar_INDELs_marker) > 0){
          
          indel_positions = sort(unique(
            c(as.integer(unlist(str_extract_all(x_cigar_INDELs_marker, '\\d+'))),
              as.integer(unlist(str_extract_all(y_cigar_INDELs_marker, '\\d+'))))
          ))
          
          x_cigar_INDELs_marker_pos = NULL
          y_cigar_INDELs_marker_pos = NULL
          
          for(position in indel_positions){
            x_cigar_INDELs_marker_pos = c(x_cigar_INDELs_marker_pos, str_extract(x_cigar_INDELs_marker, paste0(position, '(D|I)=[ATGC]+')))
            y_cigar_INDELs_marker_pos = c(y_cigar_INDELs_marker_pos, str_extract(y_cigar_INDELs_marker, paste0(position, '(D|I)=[ATGC]+')))
            
          }
          
          x_cigar_INDELs_marker_pos[is.na(x_cigar_INDELs_marker_pos)] = ''
          y_cigar_INDELs_marker_pos[is.na(y_cigar_INDELs_marker_pos)] = ''
          
          nDiffINDELs = sum(x_cigar_INDELs_marker_pos != y_cigar_INDELs_marker_pos)
          nINDELs = length(indel_positions)
          
        }else{
          
          nDiffINDELs = 0
          nINDELs = 0
          
        }
        
        # Identify indels (MODIFY THE WAY IT COUNT DELETIONS AND ICLUDE THE FIRST POSITION OF AN INSERTION INTO THE SNV COUNT)
        
        x_indels = x_fastas_marker == '' | nchar(x_fastas_marker) > 1
        y_indels = y_fastas_marker == '' | nchar(y_fastas_marker) > 1
        
        # Measure number of nucleotide differences
        
        NucDiff_x_y = sum((x_fastas_marker != y_fastas_marker) &
                            !x_indels & !y_indels)
        
        nNucSites = sum(!x_indels & !y_indels)
        
        IndelLength = sum(x_indels | y_indels)
        
        nTotalSites = length(x_fastas_marker)
        
        estimate = rbind(
          estimate,
          data.frame(
            Marker = marker,
            nVSites = NucDiff_x_y + nDiffINDELs,
            NucDiff = NucDiff_x_y,
            nNucSites = nNucSites,
            nDiffINDELs = nDiffINDELs,
            nINDELs = nINDELs,
            TotalIndelLength = IndelLength,
            nTotalSites = nTotalSites
        ))
        
      }
      
    }else{
      
      # Concatenate all fastas in both samples
      
      x_merged_fastas = unlist(x_fastas)
      y_merged_fastas = unlist(y_fastas)
      
      # Identify indels (MODIFY THE WAY IT COUNT DELETIONS AND ICLUDE THE FIRST POSITION OF AN INSERTION INTO THE SNV COUNT)
      
      x_indels = x_merged_fastas == '' | nchar(x_merged_fastas) > 1
      y_indels = y_merged_fastas == '' | nchar(y_merged_fastas) > 1
      
      # Measure number of nucleotide differences
      
      NucDiff_x_y = sum((x_merged_fastas != y_merged_fastas) &
                          !x_indels & !y_indels)
      
      nNucSites = sum(!x_indels & !y_indels)
      
      # Measure number of indel differences
      
      markers_with_INDELS = x_cigar_INDELs != '' | y_cigar_INDELs != ''
      
      x_cigar_INDELs = x_cigar_INDELs[markers_with_INDELS]
      y_cigar_INDELs = y_cigar_INDELs[markers_with_INDELS]
      
      x_marker_indel_pos = unlist(sapply(names(x_cigar_INDELs), function(marker){
        paste0(marker, ';', unlist(str_extract_all(x_cigar_INDELs[marker], '\\d+')))
      }))
      
      x_marker_indel_pos = x_marker_indel_pos[!grepl(';$', x_marker_indel_pos)]
      
      x_marker_indel = paste(unlist(sapply(names(x_cigar_INDELs), function(marker){
        
        if(x_cigar_INDELs[marker] != ''){
          
          marker_positions = unlist(str_extract_all(x_cigar_INDELs[marker], '\\d+(I|D)=[ATGC]+'))
          
          paste(paste0(marker, ';', marker_positions), collapse = '/')
        }
        
      })), collapse = '/')
      
      
      y_marker_indel_pos = unlist(sapply(names(y_cigar_INDELs), function(marker){
        paste0(marker, ';', unlist(str_extract_all(y_cigar_INDELs[marker], '\\d+')))
      }))
      
      y_marker_indel_pos = y_marker_indel_pos[!grepl(';$', y_marker_indel_pos)]
      
      y_marker_indel = paste(unlist(sapply(names(y_cigar_INDELs), function(marker){
        
        if(y_cigar_INDELs[marker] != ''){
          
          marker_positions = unlist(str_extract_all(y_cigar_INDELs[marker], '\\d+(I|D)=[ATGC]+'))
          
          paste(paste0(marker, ';', marker_positions), collapse = '/')
        }
        
      })), collapse = '/')
      
      marker_indel_positions = unique(c(
        x_marker_indel_pos,
        y_marker_indel_pos
      ))
      
      marker_indel_positions = marker_indel_positions[!grepl(';$', marker_indel_positions)]
      
      x_cigar_INDELs_marker_pos = NULL
      y_cigar_INDELs_marker_pos = NULL
      
      for(marker_indel_position in marker_indel_positions){
        
        x_cigar_INDELs_marker_pos = c(x_cigar_INDELs_marker_pos, str_extract(x_marker_indel, paste0(marker_indel_position, '(D|I)=[ATGC]+')))
        y_cigar_INDELs_marker_pos = c(y_cigar_INDELs_marker_pos, str_extract(y_marker_indel, paste0(marker_indel_position, '(D|I)=[ATGC]+')))
        
      }
      
      x_cigar_INDELs_marker_pos[is.na(x_cigar_INDELs_marker_pos)] = ''
      y_cigar_INDELs_marker_pos[is.na(y_cigar_INDELs_marker_pos)] = ''
      
      nDiffINDELs = sum(x_cigar_INDELs_marker_pos != y_cigar_INDELs_marker_pos)
      nINDELs = length(marker_indel_positions)
      
      IndelLength = sum(x_indels | y_indels) # This Value is incorrect
      
      nTotalSites = length(x_merged_fastas)
      
      estimate = data.frame(
        nVSites = NucDiff_x_y + nDiffINDELs,
        NucDiff = NucDiff_x_y,
        nNucSites = nNucSites,
        nDiffINDELs = nDiffINDELs,
        nINDELs = nINDELs,
        TotalIndelLength = IndelLength,
        nTotalSites = nTotalSites
      )
      
    }
    
    
    
    return(estimate)
  }
  
  if(parallel){
    registerDoMC(detectCores())
    pairwise_df = foreach(pair = 1:nrow(pairs), .combine = 'rbind') %dopar% {
      
      Yi_id = pairs[pair, 1]
      Yj_id = pairs[pair, 2]
      
      Yi = split_clones(genotype = gt[Yi_id,], genotype_name = Yi_id)
      
      Yj = split_clones(genotype = gt[Yj_id,], genotype_name = Yj_id)
      
      estimate = NULL
      for(haplotype_i in 1:nrow(Yi)){
        for(haplotype_j in 1:nrow(Yj)){
          
          estimate = rbind(estimate,
                           get_NucDiff(x = Yi[haplotype_i,],
                                       y = Yj[haplotype_j,],
                                       markers_table = markers_table,
                                       fasta_references = fasta_references,
                                       by_marker = by_marker
                           ))
          
          
        }
      }
      
      if(sum(is.na(estimate[,'nVSites'])) > 0){
        
        estimate = estimate[1,]
        
      }else{
        
        estimate = estimate[which.min(estimate[,'nVSites']),]
        
      }
      
      data.frame(Yi = Yi_id, Yj = Yj_id, estimate)
      
    }
    
  }else{
    
    pairwise_df = NULL
    
    for(pair in 1:nrow(pairs)){
      
      Yi_id = pairs[pair, 1]
      Yj_id = pairs[pair, 2]
      
      Yi = split_clones(genotype = gt[Yi_id,], genotype_name = Yi_id)
      
      Yj = split_clones(gt[Yj_id,], Yj_id)
      
      estimate = NULL
      for(haplotype_i in 1:nrow(Yi)){
        for(haplotype_j in 1:nrow(Yj)){
          
          if(by_marker){
            
            estimate[[paste('hap_i_', haplotype_i, '_vs_hap_j_', haplotype_j)]] = get_NucDiff(x = Yi[haplotype_i,],
                                                                                              y = Yj[haplotype_j,],
                                                                                              markers_table = markers_table,
                                                                                              fasta_references = fasta_references,
                                                                                              by_marker = by_marker
            )
            
          }else{
            
            estimate = rbind(estimate,
                             get_NucDiff(x = Yi[haplotype_i,],
                                         y = Yj[haplotype_j,],
                                         markers_table = markers_table,
                                         fasta_references = fasta_references,
                                         by_marker = by_marker
                             ))
            
          }
        }
      }
      
      if(by_marker){
        
        if(length(estimate) > 1){
          
          nVSites = unlist(sapply(names(estimate), function(comparison){
            
            sum(estimate[[comparison]][,'nVSites'], na.rm = T)
            
          }))
          
          estimate = estimate[[which.min(nVSites)]]
          
        }else{
          
          estimate = estimate[[1]]
          
        }
        
        
        
      }else{
        
        if(sum(is.na(estimate[,'nVSites'])) > 0){
          
          estimate = estimate[1,]
          
        }else{
          
          estimate = estimate[which.min(estimate[,'nVSites']),]
          
        }
        
      }
      
      
      
      
      progress(round(100*pair/nrow(pairs)))
      
      pairwise_df = rbind(pairwise_df, data.frame(Yi = Yi_id, Yj = Yj_id, estimate))
      
    }
    
  }
  
  return(pairwise_df)
  
}


## split_clones----

split_clones = function(genotype, genotype_name = NULL){
  
  if(is.null(genotype_name)){
    genotype_name = 'Yi'
  }
  
  max_nAlleles = ifelse(sum(grepl("_",genotype)) == 0,
                        1,
                        max(unlist(sapply(
                          sapply(
                            genotype[which(grepl("_",genotype))], function(x){strsplit(x, "_")}),
                          function(x) length(x)))))
  
  splited_genotype = matrix(NA,
                            nrow = max_nAlleles,
                            ncol = length(genotype),
                            dimnames = list(
                              paste0(genotype_name, '_C', 1:max_nAlleles),
                              names(genotype)
                            ))
  
  alleles = unlist(strsplit(genotype, '_'))
  cigar_pattern = grepl('(^\\d+([ATGC]|(I|D)=[ATGC]+)|^\\.$)', alleles)
  
  if(sum(cigar_pattern) > 0 & sum(cigar_pattern) == sum(!is.na(alleles))){
    
    for(clone in 1:max_nAlleles){
      splited_genotype[clone, ] = gsub("_.+", "", genotype)
      
      genotype = gsub('(^(\\d+([ATGC]|(I|D)=[ATGC]+))+|^\\.)', '', genotype)
      genotype = gsub('^_', '', genotype)
      
      if(clone > 1){
        
        splited_genotype[clone, splited_genotype[clone, ] == '' & !is.na(splited_genotype[clone, ])] = splited_genotype[1, splited_genotype[clone, ] == '' & !is.na(splited_genotype[clone, ])] 
        
      }
      
    }
    
  }else{
    
    for(clone in 1:max_nAlleles){
      splited_genotype[clone, ] = gsub("_.+", "", genotype)
      
      genotype = gsub("^\\d+", "", genotype)
      genotype = gsub("^_", "", genotype)
      
      if(clone > 1){
        
        splited_genotype[clone, splited_genotype[clone, ] == '' & !is.na(splited_genotype[clone, ])] = splited_genotype[1, splited_genotype[clone, ] == '' & !is.na(splited_genotype[clone, ])] 
        
      }
      
    }
    
  }
  
  
  return(splited_genotype)
  
}


## get_cigar_alleles----
get_cigar_alleles = function(ampseq_object, 
                             from = 'gt', # 'asv_table'
                             as = 'vector', # 'table', 'list'
                             col = 'CIGAR_masked' # 'CIGAR'
){
  
  if(from == 'gt'){
    gt = ampseq_object@gt
    
    alleles = sapply(colnames(gt), function(mhap){
      alleles = unique(unlist(strsplit(gsub(':\\d+', '',gt[,mhap]), '_')))
      
      alleles = alleles[!is.na(alleles)]
    })
    
    nalleles = sapply(alleles, function(allele){length(allele)})
    alleles = alleles[nalleles > 0]
    
    if(as == 'table'){
      
      result = NULL
      
      for(Amplicon in names(alleles)){
        
        result = rbind(result,
                       data.frame(Amplicon = Amplicon,
                                  Allele = alleles[[Amplicon]]))
      }
      
    }else if(as == 'vector'){
      
      result = NULL
      
      for(Amplicon in names(alleles)){
        
        result = c(result,
                   paste0(Amplicon,
                          ';',
                          alleles[[Amplicon]]))
      }
      
    }else if(as == 'list'){
      
      result = alleles
      
    }
    
  }else if(from == 'asv_table'){
    
    asv_table = ampseq_object@asv_table
    
    alleles = asv_table[,c('Amplicon', col)]
    names(alleles) = c('Amplicon', 'Allele')
    
    alleles = alleles[!is.na(alleles[['Allele']]),]
    
    if(as == 'table'){
      
      result = alleles %>% arrange(Amplicon)
      
    }else if(as == 'vector'){
      
      result = sort(paste0(alleles[['Amplicon']], ';', alleles[['Allele']]))
      
    }else if(as == 'list'){
      
      result = NULL
      
      for(Amplicon in sort(unique(alleles[['Amplicon']]))){
        
        result[[Amplicon]] = alleles[alleles[['Amplicon']] == Amplicon, ][['Allele']]
        
      }
      
    }
    
  }
  
  
  return(result)
  
  
  
}


consistency_between_gt_and_asvtab = function(ampseq_object, col = 'CIGAR_masked'){
  
  cigars_gt = get_cigar_alleles(ampseq_object, 
                                    from = 'gt',
                                    as = 'vector')
  
  cigars_asvtab = get_cigar_alleles(ampseq_object, 
                                        from = 'asv_table',
                                        as = 'vector',
                                        col = col)
  
  if(sum(!(cigars_gt %in% cigars_asvtab)) > 0 & 
     sum(!(cigars_asvtab %in% cigars_gt)) > 0){
    
    
    
    cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                 paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                 '\n\n',
                 'Cigar strings in asv_table that are not present in gt matrix are:\n',
                 paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                 '\n'
    ))
    
    stop('There are cigar strings that are different in the gt and the asv_table')
    
  }else if(sum(!(cigars_gt %in% cigars_asvtab)) > 0){
    
    cat(paste0('Cigar strings in gt matrix that are not present in asv_table are:\n',
                 paste(cigars_gt[!(cigars_gt %in% cigars_asvtab)], collapse = '\n'),
                 '\n'))
    
    stop('There are cigar strings in the gt that are not present in the asv_table')
    
  }else if(sum(!(cigars_asvtab %in% cigars_gt)) > 0){
    
    cat(paste0('Cigar strings in asv_table that are not present in gt matrix are:\n',
                 paste(cigars_asvtab[!(cigars_asvtab %in% cigars_gt)], collapse = '\n'),
                 '\n'))
    
    stop('There are cigar strings in the asv_table that are not present in the gt')
    
  }else{
    print('cigar strings are consistent between gt and asv_table')
  }
  
}

consistency_between_seqs_and_cigars = function(obj, 
                                               ref_fasta, 
                                               cigar_string_col = 'CIGAR',
                                               fasta_from = 'asv_seqs'){
  
  test_fasta = cigar_strings2fasta(obj = obj@asv_table, 
                                   ref_fasta = ref_fasta,
                                   cigar_string_col = cigar_string_col
                                   )
  
  
  cigars_not_in_fastas = which(!(test_fasta$fasta %in% as.character(slot(obj, fasta_from))))
  fastas_not_in_cigars = which(!(as.character(slot(obj, fasta_from)) %in% test_fasta$fasta))
  
  if(length(fastas_not_in_cigars) == 0 & length(cigars_not_in_fastas) == 0){
    
    print('There is consistency between the fasta sequences and their cigar strings')
    
  }else{
    
    if(length(fastas_not_in_cigars) > 0){
      
      cat(paste0('There are fasta sequences not found in the cigar strings:\n',
                 paste(names(slot(obj, fasta_from))[fastas_not_in_cigars]), collapse = ', '))
      
    }
    
    if(length(cigars_not_in_fastas) > 0){
      
      cat(paste0('There are cigar strings not found in the fasta sequences:\n',
                 paste(paste(test_fasta$Amplicon[cigars_not_in_fastas], test_fasta$CIGAR[cigars_not_in_fastas], sep = ','), collapse = '; ')
      ))
      
    }
    
  }
  
  
}


## read_rGenome ----
read_rGenome = function(file = NULL, format = 'excel', sep = '\t'){
  
  rGenome_object = rGenome()
  
  if(format == 'excel'){
    
    temp_wb = loadWorkbook(file)
    
    for(sheet in getSheets(temp_wb)){
      if(sheet == 'gt'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = as.matrix(temp_sheet[,-1])
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet == 'metadata'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet == 'loci_table'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = temp_sheet[,-1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, sheet, check = TRUE) = temp_sheet
        
      }
    }
    
  }else if(format == 'tsv'){
    
    for(sheet in list.files(file)){
      if(sheet == 'gt.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = as.matrix(temp_sheet[,-1])
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet == 'metadata.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet == 'loci_table.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = temp_sheet[,-1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }
    }
    
  }else if(format == 'json'){
    # In development
  }
  
  return(rGenome_object)
  
}



write_rGenome = function(rGenome_object, format = c('excel', 'tsv', 'json'), 
                         name = 'wb.xlsx',
                         sep = '\t'){
  
  if(format == 'excel'){
    
    if(file.exists(name)){
      system(paste0('rm ', name))
    }
    
    excel_wb = loadWorkbook(name, create = T)
    
    for(temp_slot in c('gt', 
                       'metadata', 
                       'loci_table')){
      
      if(temp_slot %in% c('gt', 'loci_table')){
        
        temp_sheet = data.frame(Position_id = rownames(slot(rGenome_object, temp_slot)),
                                as.data.frame(slot(rGenome_object, temp_slot)))
        
      }else{
        
        if(!is.null(slot(rGenome_object, temp_slot))){
          temp_sheet = as.data.frame(slot(rGenome_object, temp_slot))
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
    
  }else if(format == 'tsv'){
    
    if(file.exists(name)){
      system(paste0('rm -r ', name))
    }
    
    system(paste0('mkdir ', name))
    
    for(temp_slot in c('gt', 
                       'loci_table',
                       'metadata')){
      
      
      
      if(temp_slot %in% c('gt', 'loci_table')){
        
        temp_sheet = data.frame(Position_id = rownames(slot(rGenome_object, temp_slot)),
                                as.data.frame(slot(rGenome_object, temp_slot)))
        
      }else{
        
        if(!is.null(slot(rGenome_object, temp_slot))){
          temp_sheet = as.data.frame(slot(rGenome_object, temp_slot))
        }else{
          temp_sheet = NULL
        }
        
        
      }
      
      if(!is.null(temp_sheet)){
        
        write.table(temp_sheet, paste0(file.path(name, temp_slot), '.tsv'), 
                    quote = F, 
                    row.names = F,
                    sep = sep)
        
      }
      
      
    }
    
    
  }else if(format == 'json'){
    # In development
    
  }
  
}

## read_rGenome ----
read_rGenome = function(file = NULL, format = 'excel', sep = '\t'){
  
  rGenome_object = rGenome()
  
  if(format == 'excel'){
    
    temp_wb = loadWorkbook(file)
    
    for(sheet in getSheets(temp_wb)){
      if(sheet == 'gt'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = as.matrix(temp_sheet[,-1])
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet == 'metadata'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet == 'loci_table'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = temp_sheet[,-1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, sheet, check = TRUE) = temp_sheet
        
      }
    }
    
  }else if(format == 'tsv'){
    
    for(sheet in list.files(file)){
      if(sheet == 'gt.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = as.matrix(temp_sheet[,-1])
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet == 'metadata.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet == 'loci_table.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = temp_sheet[,-1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }
    }
    
  }else if(format == 'json'){
    # In development
  }
  
  return(rGenome_object)
  
}


find_homoplymer_regions = function(markers, ref_fasta, homopolymer_length){
  
  mhaps = markers
  
  mhaps$homopolymer_regions = NA
  mhaps$homopolymer_anchors = NA
  
  mhaps$STR_regions = NA
  
  ref_sequences = readDNAStringSet(ref_fasta)
  
  if(sum(grepl('(/|-|:)', names(ref_sequences))) > 0){
    
    names(ref_sequences) = gsub('(/|-|:)', '_', names(ref_sequences))
    print('The cigar and ampseq formats do not allows the symbols "-", "/", nor ":" in the name of the amplicons. All these symbols will be replaced by "_" in the names of the fasta sequences.')
    
  }
  
  homopolymer_pattern = '(A{length,}|T{length,}|G{length,}|C{length,})'
  
  homopolymer_pattern = gsub('length', homopolymer_length, homopolymer_pattern)
  
  
  for(mhap in mhaps$amplicon){
    
    # Identification of Homopolymers
    
    homopolymers = unlist(str_extract_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern))
    
    if(length(homopolymers) > 0){
      
      homopolymers_location = str_locate_all(as.character(ref_sequences[[mhap]]), homopolymer_pattern)
      
      mhaps[mhaps$amplicon == mhap, ][['homopolymer_regions']] = paste(paste(homopolymers,
                                                                             paste(homopolymers_location[[1]][,'start'],
                                                                                   homopolymers_location[[1]][,'end'], sep = '-'), sep = ':'),
                                                                       collapse = ',')
      
      homopolymer_starts = homopolymers_location[[1]][,'start']
      
      anchor_position = NULL
      for(homopolymer in 1:length(homopolymers_location[[1]][,'start'])){
        
        homopolymer_start = homopolymers_location[[1]][,'start'][homopolymer]
        
        if(homopolymer_start > 1){
          
          anchor_position = c(anchor_position, 
                              as.character(Biostrings::substr(ref_sequences[[mhap]], 
                                                              start = homopolymer_start - 1, 
                                                              stop = homopolymer_start - 1)))
        }else{
          anchor_position = c(anchor_position, '')
        }
        
      }
      
      if(!is.null(anchor_position)){
        homopolymer_anchors = paste0(homopolymers_location[[1]][,'start'] - 1, anchor_position)
        mhaps[mhaps$amplicon == mhap, ][['homopolymer_anchors']] = paste(homopolymer_anchors, collapse = ',')  
      }
      
    }
    
    
    
    
  }
  
  return(mhaps)
  
}

# END ----
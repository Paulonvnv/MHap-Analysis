
# cigar S4class and create_cigar----

## cigar S4 class

setClass('cigar', slots = c(
  cigar_table = "ANY",
  metadata = "ANY"
))

## cigar constructor
create_cigar = function(cigar_table = NULL,
                   metadata = NULL){
  obj = new('cigar')
  obj@cigar_table = cigar_table
  obj@metadata = metadata
  
  return(obj)
}


# read_cigar_tables----

read_cigar_tables = function(paths = NULL, files = NULL, sample_id_pattern = '.'){
  
  # Uploading data sets and merge CIGARs tables
  # and generate a metadata table for sample names, runs and sample types (controls vs samples of interest)
  cigar_tables = NULL
  metadata = NULL
  
  if(!is.null(paths)&is.null(files)){
    for(run in list.files(paths)){
      
      if(file.exists(file.path(paths, run, "dada2/run_dada2/CIGARVariants.out.tsv"))){
        cigar_run = read.table(file.path(paths, run, "dada2/run_dada2/CIGARVariants.out.tsv"), header = T)
      }else if(file.exists(file.path(paths, run, "dada2/run_dada2/CIGARVariants_Bfilter.out.tsv"))){
        cigar_run = read.table(file.path(paths, run, "dada2/run_dada2/CIGARVariants_Bfilter.out.tsv"), header = T)
      }else{
        print(paste0('Cigar file ', file.path(paths, run, "dada2/run_dada2/CIGARVariants_Bfilter.out.tsv"), ' not found'))
      }
      samples = gsub('_S\\d+$','', cigar_run[,1])
      samples[!grepl(sample_id_pattern,samples)] = paste(samples[!grepl(sample_id_pattern,samples)], run, sep = "_")
      cigar_run = cigar_run[, -1]
      cigar_run = apply(cigar_run, 2, function(x) as.integer(x))
      rownames(cigar_run) = paste(1:length(samples), run, samples, sep = '/')
      cigar_tables[[run]] = cigar_run
      
      rm(list = c("run", "cigar_run", "samples"))
    }
  }else if(is.null(paths)&!is.null(files)){
    for(file in files){
      cigar_run = read.table(file, header = T)
      samples = gsub('_S\\d+$','', cigar_run[,1])
      samples[!grepl(sample_id_pattern,samples)] = paste0(samples[!grepl(sample_id_pattern,samples)], '_file', which(files == file))
      cigar_run = cigar_run[, -1]
      cigar_run = apply(cigar_run, 2, function(x) as.integer(x))
      rownames(cigar_run) = paste(1:length(samples), file, samples, sep = '/')
      cigar_tables[[file]] = cigar_run
      
      rm(list = c("file", "cigar_run", "samples"))
    }
  }
  
  # Combining cigar tables
  
  cigar_table = NULL
  
  for(run in names(cigar_tables)){
    temp_cigar_table = data.frame(Sample_id = rownames(cigar_tables[[run]]), cigar_tables[[run]])
    names(temp_cigar_table) = c("Sample_id", colnames(cigar_tables[[run]]))
    temp_cigar_table %<>% pivot_longer(cols = all_of(names(temp_cigar_table)[-1]), names_to = "alleles", values_to = "counts")
    cigar_table = rbind(cigar_table, temp_cigar_table)
    rm(run)
  }
  
  cigar_table %<>% pivot_wider(names_from = "alleles", values_from = "counts")
  sample_names = gsub("^.+/","",cigar_table[['Sample_id']])
  
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
  
  metadata = data.frame(Sample_id = gsub("^.+/","",cigar_table[['Sample_id']]),
                        run = gsub('/.+$', '',gsub("^\\d+/","",cigar_table[['Sample_id']])),
                        order_in_plate = gsub("/.+$","",cigar_table[['Sample_id']]), typeofSamp = NA)
  
  
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
                              metadata = metadata)
  
  return(cigar_object)
  
}


# ampseq S4class and create_ampseq----

## cigar S4 class

setClass('ampseq', slots = c(
  gt = "ANY",
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

# cigar2ampseq----

cigar2ampseq = function(cigar_object, min_abd = 1, min_ratio = .1, markers = NULL, markers_pattern = NULL, remove_controls = F){
  
  cigar_table = cigar_object@cigar_table
  
  if(!is.null(markers)){
    ampseq_loci_vector = markers[["amplicon"]]
    
  }else if(!is.null(markers_pattern)){
    all_alleles = gsub(markers_pattern, '', colnames(cigar_table))
    all_alleles = paste0(paste(c(all_alleles[all_alleles != '.'], '\\.'), collapse = '|'), '$')
    ampseq_loci_vector = unique(gsub(all_alleles, '',colnames(cigar_table)))
    
  }else if(is.null(markers) & is.null(markers)){
    ampseq_loci_vector = unique(sapply(strsplit(colnames(cigar_table), ","), function(x) x[1]))
  }
  
  ampseq_loci_abd_table = matrix(NA, nrow = nrow(cigar_table), ncol = length(ampseq_loci_vector), dimnames = list(rownames(cigar_table), ampseq_loci_vector))
  
  for(sample in rownames(ampseq_loci_abd_table)){
    for(locus in colnames(ampseq_loci_abd_table)){
      alleles = cigar_table[sample, grepl(paste("^",locus, sep = ""), colnames(cigar_table))]
      
      if(length(alleles) == 1){
        names(alleles) <- paste(locus, ",.", sep = "")
      }
      
      if(length(alleles[which(alleles > min_abd)]) == 1){
        
        if(!is.null(markers_pattern)){
          ampseq_loci_abd_table[sample, locus] = paste(gsub(locus,'',names(alleles[which(alleles > min_abd)]), ","), alleles[which(alleles > min_abd)], sep = ":")
        }else{
          ampseq_loci_abd_table[sample, locus] = paste(gsub(paste0(locus, '.'),'',names(alleles[which(alleles > min_abd)]), ","), alleles[which(alleles > min_abd)], sep = ":")
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
                                  metadata = metadata,
                                  markers = markers,
                                  loci_performance = NULL,
                                  pop_summary = NULL)
    
    return(ampseq_object)
    
  }
  
}

# plot_coverage----

plot_coverage = function(ampseq, variable){
  coverage = sapply(1:ncol(ampseq@gt), function(mhap){
    
    temp_mhap = strsplit(ampseq@gt[,mhap], '_')
    
    sapply(1:length(temp_mhap), function(sample){
      
      sum(as.numeric(gsub('[^*]+:','', temp_mhap[[sample]])))
      
    })
  })
  
  colnames(coverage) = colnames(ampseq@gt)
  rownames(coverage) = rownames(ampseq@gt)
  
  coverage = as.data.frame(coverage)
  coverage$Sample_id = rownames(coverage)
  
  
  coverage$var = ampseq@metadata[[variable]]
  
  coverage %<>% pivot_longer(cols = all_of(colnames(ampseq@gt)),
                             names_to = 'Markers',
                             values_to = 'Read_depth')
  
  sample_order = coverage %>%
    group_by(Sample_id) %>%
    summarise(Read_depth = sum(Read_depth, na.rm = T)) %>%
    arrange(Read_depth) %>%
    select(Sample_id) %>% unlist
  
  
  plot_read_depth_heatmap = coverage %>% 
    ggplot(aes(x = Markers, y = factor(Sample_id, levels = sample_order), fill = log10(Read_depth+1)))+
    geom_tile()+
    scale_fill_gradient(low="white", high="red",
                        breaks = 1:ceiling(log10(max(coverage$Read_depth, na.rm = T))),
                        labels = 10^(1:ceiling(log10(max(coverage$Read_depth, na.rm = T)))))+
    facet_wrap(var~., ncol = 1, scales = 'free_y')+
    labs(y = 'Samples',
         fill = "Read depth")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 6),
          axis.text.y = element_text(size = 6))
  
  
  plot_read_depth_violin = coverage %>% 
    ggplot(aes(x = var, y = log10(Read_depth + 1), color = var))+
    geom_violin()+
    geom_jitter(alpha = .1, width = .2)+
    scale_y_continuous(breaks = log10(10^(0:ceiling(log10(max(coverage$Read_depth, na.rm = T))))+1),
                       labels = 10^(0:ceiling(log10(max(coverage$Read_depth, na.rm = T)))))+
    labs(x = variable, y = 'Read depth')+
    theme_bw()+
    theme(legend.position = 'none')
  
  
  coverage_by_sample = coverage %>% 
    summarise(Total_Read_depth = sum(Read_depth, na.rm = T), .by = c(var, Sample_id))
  
  plot_read_depth_violin_by_sample = coverage_by_sample %>% 
    ggplot(aes(x = var, y = log10(Total_Read_depth + 1), color = var))+
    geom_violin()+
    geom_jitter(alpha = .1, width = .2)+
    scale_y_continuous(breaks = log10(10^(0:ceiling(log10(max(coverage_by_sample$Total_Read_depth, na.rm = T))))+1),
                       labels = 10^(0:ceiling(log10(max(coverage_by_sample$Total_Read_depth, na.rm = T)))))+
    labs(x = variable, y = 'Sample Read depth')+
    theme_bw()+
    theme(legend.position = 'none')
  
  return(list(plot_read_depth_heatmap = plot_read_depth_heatmap,
              plot_read_depth_violin = plot_read_depth_violin,
              plot_read_depth_violin_by_sample = plot_read_depth_violin_by_sample))
  
}


# filter_samples----


filter_samples = function(ampseq, v){
  
  obj = ampseq
  
  obj@gt = ampseq@gt[v,]
  obj@metadata = ampseq@metadata[v,]
  
  return(obj)
}

# locus_amplification_rate----

locus_amplification_rate = function(ampseq_object, threshold = .65, chr_lengths = c(640851,
                                                                                      947102,
                                                                                      1067971,
                                                                                      1200490,
                                                                                      1343557,
                                                                                      1418242,
                                                                                      1445207,
                                                                                      1472805,
                                                                                      1541735,
                                                                                      1687656,
                                                                                      2038340,
                                                                                      2271494,
                                                                                      2925236,
                                                                                      3291936)){
  
  ampseq_loci_abd_table = ampseq_object@gt
  
  loci_performance = data.frame(loci = colnames(ampseq_loci_abd_table),
                                loci_ampl_rate = apply(ampseq_loci_abd_table, 2, function(x) 1 - sum(is.na(x))/length(x)))
  
  
  all_loci_performance_plot = loci_performance %>% ggplot(aes(loci_ampl_rate))+
    geom_histogram(binwidth = 0.05, alpha = .7, boundary = 1)+
    labs(title = paste(sum(loci_performance[["loci_ampl_rate"]] > threshold), "Loci have an amplification rate >",threshold),
         x = "Amplification rate",
         y = "# of Loci")+
    geom_vline(xintercept = threshold) +
    theme_bw()
  
  amplification_rate_per_locus = ggplot()+
    geom_segment(data = data.frame(x = 0, y = 1:14, xend =chr_lengths, yend = 1:14),
                 aes(x=x, y=y, xend=xend, yend=yend), alpha = .5)+
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
    theme(legend.position = c(.9,.4))+
    scale_y_continuous(breaks = 1:14)+
    scale_color_continuous(type = "viridis")
  
  ampseq_loci_abd_table_discarded_loci = ampseq_loci_abd_table[, !(colnames(ampseq_loci_abd_table) %in% loci_performance[loci_performance[["loci_ampl_rate"]] > threshold,][["loci"]])]
  ampseq_loci_abd_table = ampseq_loci_abd_table[, colnames(ampseq_loci_abd_table) %in% loci_performance[loci_performance[["loci_ampl_rate"]] > threshold,][["loci"]]]
  
  markers = ampseq_object@markers
  
  discarded_markers = markers[loci_performance[["loci_ampl_rate"]] <= threshold ,]
  markers = markers[loci_performance[["loci_ampl_rate"]] > threshold ,]
  
  markers[["distance"]] = Inf
  
  for(chromosome in levels(as.factor(markers[["chromosome"]]))){
    for(amplicon in 1:(nrow(markers[markers[["chromosome"]] == chromosome,])-1)){
      markers[markers[["chromosome"]] == chromosome,][amplicon, "distance"] = markers[markers[["chromosome"]] == chromosome,][amplicon + 1, "pos"] - markers[markers[["chromosome"]] == chromosome,][amplicon, "pos"]
    }
  }
  
  
  loci_performance_complete = loci_performance
  loci_performance = loci_performance[loci_performance[["loci_ampl_rate"]] > threshold ,]
  
  ampseq_object@gt = ampseq_loci_abd_table
  ampseq_object@markers = markers
  ampseq_object@loci_performance = loci_performance
  ampseq_object@discarded_loci = list(gt = ampseq_loci_abd_table_discarded_loci,
                                           loci_performance = loci_performance_complete,
                                           markers = discarded_markers)
  ampseq_object@plots[["all_loci_amplification_rate"]] = all_loci_performance_plot
  ampseq_object@plots[["amplification_rate_per_locus"]] = amplification_rate_per_locus
  return(ampseq_object)
  
}

# sample_amplification_rate----

sample_amplification_rate = function(ampseq_object, threshold = .8){
  
  metadata = ampseq_object@metadata
  ampseq_loci_abd_table = ampseq_object@gt
  loci_performance = ampseq_object@loci_performance
  
  metadata[["sample_ampl_rate"]] = apply(ampseq_loci_abd_table, 1, function(x) 1 - sum(is.na(x))/length(x))
  
  
  all_samples_performance_plot = metadata %>% ggplot(aes(sample_ampl_rate))+
    geom_histogram(binwidth = 0.05, alpha = .7, boundary = 1)+
    labs(title = paste(sum(metadata[["sample_ampl_rate"]] > threshold), "Samples have an amplification rate >",threshold),
         x = "Amplification rate",
         y = "# of samples")+
    geom_vline(xintercept = threshold) +
    theme_bw()
  
  ampseq_loci_abd_table_discarded_samples = ampseq_loci_abd_table[!(rownames(ampseq_loci_abd_table) %in% metadata[metadata[["sample_ampl_rate"]] > threshold ,][["Sample_id"]]),]
  ampseq_loci_abd_table = ampseq_loci_abd_table[rownames(ampseq_loci_abd_table) %in% metadata[metadata[["sample_ampl_rate"]] > threshold ,][["Sample_id"]],]
  
  metadata_complete = metadata
  metadata = metadata[metadata[["sample_ampl_rate"]] > threshold ,]
  
  
  loci_performance[["loci_ampl_rate2"]] = apply(ampseq_loci_abd_table, 2, function(x) 1- sum(is.na(x))/length(x))
  
  
  loci_performance_plot = loci_performance %>% ggplot(aes(loci_ampl_rate2))+
    geom_histogram(binwidth = 0.05, alpha = .7, boundary = 1)+
    labs(title = "Loci amplification rate",
         x = "Amplification rate",
         y = "# of Loci")+
    geom_vline(xintercept = threshold) +
    theme_bw()
  
  
  ampseq_object@gt = ampseq_loci_abd_table
  ampseq_object@metadata = metadata
  ampseq_object@loci_performance = loci_performance
  
  ampseq_object@discarded_samples = list(gt = ampseq_loci_abd_table_discarded_samples,
                                              metadata = metadata_complete)
  
  ampseq_object@plots[["loci_amplification_rate"]] = loci_performance_plot
  ampseq_object@plots[["samples_amplification_rate"]] = all_samples_performance_plot
  
  return(ampseq_object)
  
  
}

# haplotypes_respect_to_reference----

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
                                              fasta_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7_Genome.fasta"){
  library(ape)
  library(Biostrings)
  
  # Call 3D7 genome and gff---
  # Call 3D7 reference genome and its corresponding anotation in the gff file
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
  moi_loci_abd_table = ampseq_object@gt[,grep(paste(gene_names, collapse = "|"), colnames(ampseq_object@gt))]
  
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
  
  haplotypes_respect_to_reference = list(markers_of_interest = markers_of_interest,
                                         loci_dna_table = moi_loci_dna_table,
                                         loci_aa_table = moi_loci_aa_table,
                                         ref_seqs = ref_seqs)
  
  return(haplotypes_respect_to_reference)
  
}

# drug_resistant_haplotypes----

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
                                     gff_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7.gff",
                                     fasta_file = "reference/3D7/PlasmoDB-59_Pfalciparum3D7_Genome.fasta",
                                     variables = c('samples', 'Population', 'quarter_of_collection'),
                                     filters = c('Population;Buenaventura,Quibdo,Guapi',
                                                 'quarter_of_collection;2020-Q4,2021-Q1,2021-Q2,2021-Q3,2021-Q4'),
                                     na.var.rm = FALSE){
  
  # Call reference alleles
  drugR_reference_alleles = read.csv(reference_alleles)
  
  drugR_reference_alleles['reference'] = gsub('[0-9]+[A-Z]$', '', drugR_reference_alleles$Mutation)
  drugR_reference_alleles['mutant'] = gsub('^[A-Z][0-9]+', '', drugR_reference_alleles$Mutation)
  drugR_reference_alleles['position'] = gsub('^[A-Z]', '', gsub('[A-Z]$', '', drugR_reference_alleles$Mutation))
  
  # Define haplotypes respect to a reference genome
  
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
                         ncol = length(unique(drugR_reference_alleles$Gene_Id)),
                         dimnames = list(rownames(loci_aa_table),
                                         unique(drugR_reference_alleles$Gene_Id)))
  
  phenotype_table = matrix(NA,
                           nrow = nrow(loci_aa_table),
                           ncol = length(unique(drugR_reference_alleles$Gene_Id)),
                           dimnames = list(rownames(loci_aa_table),
                                           unique(drugR_reference_alleles$Gene_Id)))
  
  for(gene in unique(drugR_reference_alleles$Gene_Id)){ # For each gene
    
    for(amplicon in drug_markers[drug_markers$gene_ids == gene, 'amplicon']){ # for each amplicon in the gene
      
      if(amplicon %in% colnames(loci_aa_table)){
        
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
        tested_amplicon = drug_markers$amplicon %in% amplicon
        
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
            
            
          }else{
            
            sample_positions = NULL
            sample_alleles = NULL
            
            for(position in positions){
              
              sample_positions = c(sample_positions, position)
              
              sample_alleles = c(sample_alleles, '?')
              
            }
            
          }
          
          reference_alleles = NULL
          
          for(position in sample_positions){
            
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
          
          aacigar_haplotype = data.frame(position = as.integer(sample_positions),
                                         aacigar_haplotype = paste(reference_alleles, sample_positions, sample_alleles, sep = ""))
          
          aacigar_haplotype[['aacigar_haplotype']] = ifelse(aacigar_haplotype[['position']] %in% drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene,'position'],
                                                            aacigar_haplotype[['aacigar_haplotype']],
                                                            tolower(aacigar_haplotype[['aacigar_haplotype']]))
          
          sample_phenotype = NULL
          
          for(position in sort(aacigar_haplotype[['position']])){
            
            if(nrow(drugR_reference_alleles[drugR_reference_alleles[['Gene_Id']] == gene &
                                            drugR_reference_alleles[['position']] == position,]) != 0){
              
              sample_allele = gsub('[A-Z][0-9]+', '', aacigar_haplotype[aacigar_haplotype[['position']] == position, 'aacigar_haplotype'])
              
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
                                                        drugR_reference_alleles[['mutant']] == sample_allele,'Anotation']
                  
                  sample_phenotype = c(sample_phenotype, phenotype)
                  
                }else{
                  
                  if(sample_allele != '?'){
                    
                    sample_phenotype = c(sample_phenotype, paste0(sample_allele, ' variant unreported for position ', position))
                    
                  }else{
                    
                    sample_phenotype = c(sample_phenotype, paste0(amplicon, ' amplicon did not amplify'))
                    
                  }
                  
                }
                
              }
              
            }else{
              
              reference_allele = gsub('[0-9]+[A-Z]', '', aacigar_haplotype[aacigar_haplotype[['position']] == position,'aacigar_haplotype'], ignore.case = T)
              
              alternative_allele = gsub('[A-Z][0-9]+', '', aacigar_haplotype[aacigar_haplotype[['position']] == position,'aacigar_haplotype'], ignore.case = T)
              
              if(reference_allele != alternative_allele){
                
                sample_phenotype = c(sample_phenotype, paste0(aacigar_haplotype[aacigar_haplotype[['position']] == position,'aacigar_haplotype'], ' unreported polymorphims respect to 3D7'))  
                
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
  
  mon_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) == 0),]
  
  poly_aacigar_table = aacigar_table[(apply(aacigar_table, 1, function(i){sum(grepl("\\|",i))}) != 0),]
  
  poly1 = gsub('\\|([a-z])', '', poly_aacigar_table, ignore.case = T)
  
  poly2 = gsub('([a-z]\\|)', '', poly_aacigar_table, ignore.case = T)
  
  extended_aacigar_table = rbind(mon_aacigar_table,
                                 poly1,
                                 poly2)
  
  extended_aacigar_table = data.frame(samples = rownames(extended_aacigar_table), extended_aacigar_table)
  
  metadata = ampseq_object@metadata[,variables]
  
  names(metadata) = c('samples', 'var1', 'var2')
  
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
  
  
  samples_pop_quarter = extended_aacigar_table %>% group_by(var1, var2)%>%
    summarise(count = nlevels(as.factor(samples))) 
  
  haplotype_counts = extended_aacigar_table %>% group_by(gene_names, var1, var2, haplotype)%>%
    summarise(count = n())
  
  if(!is.null(filters)){
    filters = strsplit(filters,';')
    for(temp_filter in 1:length(filters)){
      if(which(variables == filters[[temp_filter]][1]) == 2){
        haplotype_counts %<>% filter(var1 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
        samples_pop_quarter %<>% filter(var1 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
      }else if(which(variables == filters[[temp_filter]][1]) == 3){
        haplotype_counts %<>% filter(var2 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
        samples_pop_quarter %<>% filter(var2 %in% strsplit(filters[[temp_filter]][2],',')[[1]])
      }
    }
  }
  
  
  haplotype_counts$freq = NA
  
  for(gene in levels(as.factor(haplotype_counts$gene_names))){
    for(Pop in levels(as.factor(haplotype_counts[haplotype_counts$gene_names == gene,][['var1']]))){
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
        # samples_pop_quarter[samples_pop_quarter$var1 == Pop&
        #                       samples_pop_quarter$var2 == date,][['count']]
        
      }
    }
  }
  
  haplotype_counts %<>% mutate(gene_haplo = paste(gene_names, haplotype, sep = ": "))
  
  nhaplo = haplotype_counts %>% group_by(gene_names)%>% summarise(nhaplo = nlevels(as.factor(haplotype))) 
  
  max_ncolors = brewer.pal(max(nhaplo$nhaplo), 'Set3')
  
  colors = NULL
  
  for(n in 1:nrow(nhaplo)){
    colors = c(colors, max_ncolors[1:nhaplo[n,][['nhaplo']]])
  }
  
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
  
  
  drug_resistant_hap_list = list(aacigar_table = aacigar_table,
                                 phenotype_table = phenotype_table,
                                 haplotype_counts =haplotype_counts,
                                 haplo_freq_plot = haplo_freq_plot)
  
  return(drug_resistant_hap_list)
  
}

# get_polygenomic----

get_polygenomic = function(ampseq_object, strata, update_popsummary = T, na.rm = FALSE, filters = NULL){
  
  library(Hmisc)
  
  ampseq_loci_abd_table = ampseq_object@gt
  metadata = ampseq_object@metadata
  loci_performance = ampseq_object@loci_performance
  
  if(na.rm){
    ampseq_loci_abd_table = ampseq_loci_abd_table[!(is.na(metadata[[strata]]) | grepl('NA',metadata[[strata]])),]
    metadata = metadata[!(is.na(metadata[[strata]]) | grepl('NA',metadata[[strata]])),]
  }else if(length(metadata[is.na(metadata[[strata]]) | grepl('NA',metadata[[strata]]),][[strata]])>0){
    metadata[is.na(metadata[[strata]]) | grepl('NA',metadata[[strata]]),][[strata]] = 'missing data'
  }
  
  if(!is.null(filters)){
    ampseq_loci_abd_table = ampseq_loci_abd_table[grepl(filters,metadata[[strata]]),]
    metadata = metadata[grepl(filters,metadata[[strata]]),]
  }
  
  polygenomic = NULL
  
  for(sample in rownames(ampseq_loci_abd_table)){
    polygenomic = rbind(polygenomic, data.frame(NPolyLoci = sum(grepl("_",ampseq_loci_abd_table[sample, ])),
                                                Polyloci = paste(names(ampseq_loci_abd_table[sample, ])[which(grepl("_",ampseq_loci_abd_table[sample, ]))], collapse = "/"),
                                                alleles_at_loci = paste(ampseq_loci_abd_table[sample, ][which(grepl("_",ampseq_loci_abd_table[sample, ]))], collapse = "/"),
                                                nalleles_per_loci = paste(
                                                  sapply(
                                                    sapply(
                                                      ampseq_loci_abd_table[sample, ][which(grepl("_", ampseq_loci_abd_table[sample, ]))], function(x){strsplit(x, "_")}),
                                                    function(x) length(x)),
                                                  collapse = "/"),
                                                coi = ifelse(sum(grepl("_",ampseq_loci_abd_table[sample, ])) == 0,
                                                             1,
                                                             max(unlist(sapply(
                                                               sapply(
                                                                 ampseq_loci_abd_table[sample, ][which(grepl("_",ampseq_loci_abd_table[sample, ]))], function(x){strsplit(x, "_")}),
                                                               function(x) length(x)))))))
  }
  
  pop_summary = data.frame(
    pop = "Total",
    n = nrow(polygenomic),
    mean_coi = mean(polygenomic[["coi"]]),
    n_poly = nrow(polygenomic[polygenomic[["coi"]] > 1,]),
    prop_poly = binconf(x = nrow(polygenomic[polygenomic[["coi"]] > 1,]),
                        n = nrow(polygenomic),
                        method = "exact")[1],
    prop_poly_lower = binconf(x = nrow(polygenomic[polygenomic[["coi"]] > 1,]),
                              n = nrow(polygenomic),
                              method = "exact")[2],
    prop_poly_upper = binconf(x = nrow(polygenomic[polygenomic[["coi"]] > 1,]),
                              n = nrow(polygenomic),
                              method = "exact")[3])
  
  
  for(pop in levels(as.factor(metadata[[strata]]))){
    pop_summary = rbind(pop_summary,
                        data.frame(
                          pop = pop,
                          n = nrow(polygenomic[metadata[[strata]] == pop,]),
                          mean_coi = mean(polygenomic[metadata[[strata]] == pop,][["coi"]]),
                          n_poly = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["coi"]] > 1,]),
                          prop_poly = binconf(x = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["coi"]] > 1,]),
                                              n = nrow(polygenomic[metadata[[strata]] == pop,]),
                                              method = "exact")[1],
                          prop_poly_lower = binconf(x = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["coi"]] > 1,]),
                                                    n = nrow(polygenomic[metadata[[strata]] == pop,]),
                                                    method = "exact")[2],
                          prop_poly_upper = binconf(x = nrow(polygenomic[metadata[[strata]] == pop & polygenomic[["coi"]] > 1,]),
                                                    n = nrow(polygenomic[metadata[[strata]] == pop,]),
                                                    method = "exact")[3]))
  }
  
  # Defining the proportion of polygenomic infections detected by loci --
  
  loci_performance[["prop_poly_detected"]] = apply(ampseq_loci_abd_table, 2, function(x) sum(grepl("_",x))/length(x))
  
  if(update_popsummary){
    ampseq_object@metadata = cbind(ampseq_object@metadata, polygenomic)
    ampseq_object@loci_performance = loci_performance
    ampseq_object@pop_summary = if(is.null(ampseq_object@pop_summary)){pop_summary}else{cbind(ampseq_object@pop_summary, pop_summary)}
    return(ampseq_object)
  }else{
    return(pop_summary)
  }
  
}


# ampseq2loci----

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

# fs_checks----

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



# estimate_r_and_k----

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


# pairwise_hmmIBD----

pairwise_hmmIBD = function(loci_object, parallel = TRUE, w = 1, n = 1){
  library(parallel)
  library(doMC)
  library(svMisc)
  
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

# plot_relatedness_distribution----

plot_relatedness_distribution = function(pairwise_relatedness = pairwise_relatedness,
                                            metadata = ampseq@metadata,
                                            Population = 'Subnational_level2',
                                            fill_color = c("firebrick3", "firebrick1", "dodgerblue1", "dodgerblue3", "gold1", "gold3"),
                                            type_pop_comparison = 'within', # c('within', 'between', 'both')
                                            ncol = 5,
                                            pop_levels = NULL
){
  
  pairwise_relatedness_l = pairwise_relatedness
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population)], by.x = 'Yi', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-5], 'Yi_Population')
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population)], by.x = 'Yj', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-6], 'Yj_Population')
  
  pairwise_relatedness_l %<>% filter(!is.na(Yi_Population), !is.na(Yj_Population))
  
  pairwise_relatedness_l %<>% mutate(Pop_comparison = case_when(
    Yi_Population == Yj_Population ~ Yi_Population,
    Yi_Population < Yj_Population ~ paste(Yi_Population, Yj_Population, sep = '-'),
    Yi_Population > Yj_Population ~ paste(Yj_Population, Yi_Population, sep = '-')
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
    grepl("-",Pop_comparison) ~ "Between",
    !grepl("-",Pop_comparison) ~ "Within"
  ))
  
  # }
  
  pairwise_relatedness_l$Pop_comparison = factor(pairwise_relatedness_l$Pop_comparison,
                                                 levels = c(sort(unique(metadata[!is.na(metadata[[Population]]),][[Population]])),
                                                            apply(combn(unique(metadata[!is.na(metadata[[Population]]),][[Population]]),2), 2, function(x){paste(sort(x), collapse = '-')})))
  
  pairwise_relatedness_l$Type_of_comparison = factor(pairwise_relatedness_l$Type_of_comparison, levels =
                                                       c('Within', 'Between'))
  
  
  
  if(type_pop_comparison == 'within'){
    plot_pairwise_relatedness_distribution = pairwise_relatedness_l %>%
      filter(Type_of_comparison == 'Within')%>%
      ggplot(aes(x = rhat, fill = factor(Pop_comparison, levels = pop_levels))) +
      geom_histogram(position = "stack", alpha = .7)+
      geom_vline(xintercept = mean(pairwise_relatedness_l$rhat), linetype = 2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      facet_wrap(~factor(Pop_comparison, levels = pop_levels), ncol = ncol,
                 scales = "free_y")+
      labs(y = "Count",
           x = "Relatedness",
           fill = "Populations")+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.position = "none")
  }else if(type_pop_comparison == 'between'){
    plot_pairwise_relatedness_distribution = pairwise_relatedness_l %>%
      filter(Type_of_comparison == 'Between')%>%
      ggplot(aes(x = rhat, fill = factor(Pop_comparison, levels = pop_levels))) +
      geom_histogram(position = "stack", alpha = .7)+
      geom_vline(xintercept = mean(pairwise_relatedness_l$rhat), linetype = 2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      facet_wrap(~factor(Pop_comparison, levels = pop_levels), ncol = ncol,
                 scales = "free_y")+
      labs(y = "Count",
           x = "Relatedness",
           fill = "Populations")+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.position = "none")
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

# plot_frac_highly_related----

plot_frac_highly_related = function(pairwise_relatedness = pairwise_relatedness,
                                       metadata = ampseq@metadata,
                                       Population = 'Population',
                                       fill_color = c("dodgerblue3",  "firebrick3", "gold3", "gray50", "gray50", "gray50"),
                                       threshold = 0.99,
                                       type_pop_comparison = 'between',
                                       pop_levels = NULL){
  
  pairwise_relatedness_l = pairwise_relatedness
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population)], by.x = 'Yi', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-5], 'Yi_Population')
  
  pairwise_relatedness_l = merge(pairwise_relatedness_l, metadata[,c('Sample_id', Population)], by.x = 'Yj', by.y = 'Sample_id', all.x = TRUE)
  
  names(pairwise_relatedness_l) = c(names(pairwise_relatedness_l)[-6], 'Yj_Population')
  

  pairwise_relatedness_l %<>% filter(!is.na(Yi_Population), !is.na(Yj_Population))
  
  pairwise_relatedness_l %<>% mutate(Pop_comparison = case_when(
    Yi_Population == Yj_Population ~ Yi_Population,
    Yi_Population < Yj_Population ~ paste(Yi_Population, Yj_Population, sep = '-'),
    Yi_Population > Yj_Population ~ paste(Yj_Population, Yi_Population, sep = '-')
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
    grepl("-",Pop_comparison) ~ "Between",
    !grepl("-",Pop_comparison) ~ "Within"
  ))
  
  # }
  
  pairwise_relatedness_l$Pop_comparison = factor(pairwise_relatedness_l$Pop_comparison,
                                                 levels = c(sort(unique(metadata[!is.na(metadata[[Population]]),][[Population]])),
                                                            apply(combn(unique(metadata[!is.na(metadata[[Population]]),][[Population]]),2), 2, function(x){paste(sort(x), collapse = '-')})))
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
      ggplot(aes(x = factor(Pop_comparison, levels = pop_levels), y = prop, fill = factor(Pop_comparison, levels = pop_levels))) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }else if(type_pop_comparison == 'between'){
    
    plot_frac_highly_related = highly_related_table %>%
      filter(Type_of_comparison == 'Between') %>%
      ggplot(aes(x = factor(Pop_comparison, levels = pop_levels), y = prop, fill = factor(Pop_comparison, levels = pop_levels))) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }else if(type_pop_comparison == 'both'){
    
    plot_frac_highly_related = highly_related_table %>%
      ggplot(aes(x = factor(Pop_comparison, levels = pop_levels), y = prop, fill = factor(Pop_comparison, levels = pop_levels))) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }
  
  frac_highly_related = list(highly_related_table = highly_related_table,
                             plot = plot_frac_highly_related)
  
  return(frac_highly_related)
  
}

# plot_frac_highly_related_over_time----

plot_frac_highly_related_over_time = function(pairwise_relatedness = pairwise_relatedness,
                                                 metadata = ampseq@metadata,
                                                 Population = c('Population', 'Quarter_of_Collection'),
                                                 fill_color = c("dodgerblue3",  "firebrick3", "gold3", "gray50", "gray50", "gray50"),
                                                 threshold = 0.99,
                                                 type_pop_comparison = 'within',## c('within', 'between', 'both')
                                                 ncol = 4,
                                                 pop_levels = NULL){
  
  pairwise_relatedness_l = pairwise_relatedness
  
  # Giving wraning message, otherwise it is faster
  # pairwise_relatedness_l %<>% mutate(
  #   Population_Yi = metadata[grep(Yi, metadata[['samples']]),][[Population[1]]],
  #   Population_Yj = metadata[grep(Yj, metadata[['samples']]),][[Population[1]]],
  #   Date_Yi = metadata[grep(Yi, metadata[['samples']]),][[Population[2]]],
  #   Date_Yj = metadata[grep(Yj, metadata[['samples']]),][[Population[2]]])
  
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
  
  
  plot_IBD_correlation_matrix = pairwise_relatedness_l %>%
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
    scale_fill_gradient2(low = "white", high = "red", limit = c(0,.5), space = "Lab",
                         name=paste0('IBD >= ', threshold))+
    theme_bw()+
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, vjust = .5),
          axis.title.y = element_text(size = 12),
          strip.text = element_text(size = 12))
  
  
  pairwise_relatedness_l %<>% mutate(Pop_comparison = case_when(
    Population_Yi < Population_Yj ~ paste(Population_Yi, Population_Yj, sep = '_'),
    Population_Yi > Population_Yj ~ paste(Population_Yj, Population_Yi, sep = '_'),
    Population_Yi == Population_Yj ~ Population_Yi),
    Type_Pop_comparison = case_when(
      Population_Yi != Population_Yj ~ "Between",
      Population_Yi == Population_Yj ~ "Within"),
    Date_comparison = case_when(
      Date_Yi != Date_Yj ~ "Between",
      Date_Yi == Date_Yj ~ "Within")
  )
  
  pairwise_relatedness_l$Pop_comparison = factor(pairwise_relatedness_l$Pop_comparison,
                                                 levels = c(sort(unique(pairwise_relatedness_l$Pop_comparison)[!grepl("_", unique(pairwise_relatedness_l$Pop_comparison))]),
                                                            sort(unique(pairwise_relatedness_l$Pop_comparison)[grepl("_", unique(pairwise_relatedness_l$Pop_comparison))])))
  
  pairwise_relatedness_l$Type_Pop_comparison = factor(pairwise_relatedness_l$Type_Pop_comparison, levels =
                                                        c('Within', 'Between'))
  
  pairwise_relatedness_l %<>% filter(Date_comparison == 'Within')
  
  frac_highly_related = pairwise_relatedness_l %>%
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
      ggplot(aes(x = Date_Yi, y = prop, fill = factor(Pop_comparison, levels = pop_levels))) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      facet_wrap(Type_Pop_comparison~factor(Pop_comparison, levels = pop_levels))+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }else if(type_pop_comparison == 'within'){
    
    plot_frac_highly_related = frac_highly_related %>%
      filter(Type_Pop_comparison == 'Within')%>%
      ggplot(aes(x = Date_Yi, y = prop, fill = factor(Pop_comparison, levels = pop_levels))) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      facet_wrap(~factor(Pop_comparison, levels = pop_levels), ncol = ncol)+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
    
  }else if(type_pop_comparison == 'between'){
    
    plot_frac_highly_related = frac_highly_related %>%
      filter(Type_Pop_comparison == 'Between')%>%
      ggplot(aes(x = Date_Yi, y = prop, fill = factor(Pop_comparison, levels = pop_levels))) + 
      geom_col(alpha = .85)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
      scale_fill_manual(values = fill_color)+
      facet_wrap(~factor(Pop_comparison, levels = pop_levels), ncol = ncol)+
      theme_bw()+
      labs(y = paste0('Proportion of highly related samples, IBD >= ', threshold))+
      theme(axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = .5),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            strip.text = element_text(size = 12),
            legend.position = "none")
    
  }
  
  frac_highly_related = list(frac_highly_related = frac_highly_related,
                             plot_frac_highly_related = plot_frac_highly_related,
                             plot_IBD_correlation_matrix = plot_IBD_correlation_matrix)
  
  return(frac_highly_related)
  
}

# plot_network----

plot_network = function(pairwise_relatedness,  
                           threshold,
                           metadata,
                           sample_id,
                           group_by,
                           levels,
                           colors){
  
  if(sum(is.na(metadata[[group_by]])) > 0){
    metadata[is.na(metadata[[group_by]]),][[group_by]] = 'missing data'
  }
  
  pairwise_relatedness_l = pairwise_relatedness
  
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
                             vertex.size=4,
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



#log_scale_histogram----

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








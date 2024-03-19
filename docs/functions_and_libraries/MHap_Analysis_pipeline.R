#!/bin/r env

library(argparse)
library(stringr)
library(rmarkdown)
# Parse arguments ----

parser = ArgumentParser()

parser$add_argument("-wd", "--wd", 
                    help="Path to input to files and output files")

parser$add_argument("-fd", "--fd", default = 'NaN',
                    help="Path to function files")

parser$add_argument("-rd", "--rd", default = 'NaN',
                    help="Path to reference files")

parser$add_argument("-cigar_paths", "--cigar_paths",  default = 'NaN',
                    help="name of the cigar file")

parser$add_argument("-cigar_files", "--cigar_files",  default = 'NaN',
                    help="name of the cigar file")

parser$add_argument("-ampseqj", "--ampseq_jsonfile",  default = 'NaN',
                    help="name of the cigar file in json format")

parser$add_argument("-ampseqe", "--ampseq_excelfile",  default = 'NaN',
                    help="name of the cigar file in excel format")

parser$add_argument("-o", "--out",
                    help="Prefix of output files")

parser$add_argument("-sample_id_pattern", "--sample_id_pattern", default = ".",
                    help="")

parser$add_argument("-markers", "--markers", default = "NaN",
                    help="")

parser$add_argument("-min_abd", "--min_abd", default = 'NaN',
                    help="")

parser$add_argument("-min_ratio", "--min_ratio", default = 'NaN',
                    help="")

parser$add_argument("-PerformanceReport", "--PerformanceReport", default = FALSE,
                    help="")


parser$add_argument("-Drug_Surveillance_Report", "--Drug_Surveillance_Report", default = FALSE,
                    help="")


parser$add_argument("-Variants_of_Interest_Report", "--Variants_of_Interest_Report", default = FALSE,
                    help="")


parser$add_argument("-gff", "--ref_gff", default = 'NaN',
                    help="name of .gff file containing coordinates of genomic regions")

parser$add_argument("-fasta", "--ref_fasta", default = 'NaN',
                    help="name of .gff file containing coordinates of genomic regions")

parser$add_argument("-reference_alleles", "--reference_alleles", default = 'NaN',
                    help="name of .csv file containing sensitive alleles respect to a drug treatment")


parser$add_argument("-gene_names", "--gene_names", default = 'NaN',
                    help="gene name on the markers csv file")

parser$add_argument("-gene_ids", "--gene_ids", default = 'NaN',
                    help="gene ids on the .gff file")

parser$add_argument("-samprate", "--sample_ampl_rate", default = 'NaN',
                    help="Minimun proportion of amplified loci by sample")

parser$add_argument("-lamprate", "--locus_ampl_rate", default = 'NaN',
                    help="Minimun proportion of amplified samples by locus")

parser$add_argument("-t", "--nTasks", default = 1,
                    help="Number of Tasks arrays to split the convertion of the vcf to the rGenome format")

parser$add_argument("-tid", "--Task_id", default = 'NaN',
                    help="Tasks ID")

parser$add_argument("-ibd", "--ibd_thres", default = 'NaN',
                    help="Minimum IBD to define highly related samples")

parser$add_argument("-nchunks", "--nchunks", default = 500,
                    help="")

parser$add_argument("-metadata", "--metadata", default = 'NaN',
                    help="string with the file name of the metadata, Sample ID's should be labeled as Sample_id")

parser$add_argument("-join_by", "--join_by", default = 'NaN',
                    help="string that indicates the variable where Sample ID's are stored. Just if  Sample ID's are not labeled as Sample_id")

parser$add_argument("-Var1", "--Variable1", default = 'NaN',
                    help="name of .bed file containing core genomic regions to keep from the VCF file")

parser$add_argument("-Var2", "--Variable2", default = 'NaN',
                    help="name of .bed file containing core genomic regions to keep from the VCF file")

parser$add_argument("-Longitude", "--Longitude", default = 'NaN',
                    help="name of .bed file containing core genomic regions to keep from the VCF file")

parser$add_argument("-Latitude", "--Latitude", default = 'NaN',
                    help="name of .bed file containing core genomic regions to keep from the VCF file")

parser$add_argument("-na_var_rm", "--na_var_rm", default = FALSE,
                    help="name of .bed file containing no core genomic regions to exclude from the VCF file")

parser$add_argument("-na_hap_rm", "--na_hap_rm", default = FALSE,
                    help="name of .bed file containing no core genomic regions to exclude from the VCF file")

parser$add_argument("-drugs", "--drugs", default = 'NaN',
                    help="")

parser$add_argument("-var_filter", "--var_filter", default = NULL,
                    help="")

parser$add_argument("-include_all_drug_markers", "--include_all_drug_markers", default = FALSE,
                    help="")

parser$add_argument("-parallel", "--parallel", default = FALSE,
                    help="name of .bed file containing no core genomic regions to exclude from the VCF file")

parser$add_argument("-ibd_step", "--ibd_step",
                    help="name of .bed file containing no core genomic regions to exclude from the VCF file")

parser$add_argument("-ibd_ncol", "--ibd_ncol", default = 4,
                    help="name of .bed file containing no core genomic regions to exclude from the VCF file")

parser$add_argument("-pop_levels", "--pop_levels", default = NULL,
                    help="name of .bed file containing no core genomic regions to exclude from the VCF file")


# Defining and checking variables ----

print("starting to parse variables")
args = parser$parse_args()

# Working directory
wd = gsub("/$", "", args$wd)

print(paste0('wd: ', wd))
setwd(wd) 

# Tools or functions directory
fd = gsub("/$", "", args$fd)
print(paste0('fd: ', fd))

# Reference files directory
rd = gsub("/$", "", args$rd)
print(paste0('rd: ', rd))

# Starting vcf file
cigar_paths = args$cigar_paths
cigar_paths = ifelse(cigar_paths == 'NaN', NA, file.path(wd, cigar_paths))
print(paste0('cigar_paths: ', cigar_paths))

# Starting vcf file
cigar_files = args$cigar_files
cigar_files = ifelse(cigar_files == 'NaN', NA, file.path(wd, cigar_files))
print(paste0('cigar_files: ', cigar_files))

# Starting vcf file
ampseq_jsonfile = args$ampseq_jsonfile
ampseq_jsonfile = ifelse(ampseq_jsonfile == 'NaN', NA, file.path(wd, ampseq_jsonfile))
print(paste0('ampseq_jsonfile: ', ampseq_jsonfile))

# Starting vcf file
ampseq_excelfile = args$ampseq_excelfile
ampseq_excelfile = ifelse(ampseq_excelfile == 'NaN', NA, file.path(wd, args$ampseq_excelfile))
print(paste0('ampseq_excelfile: ', ampseq_excelfile))

# Pattern to differentiate between samples of interest and controls
sample_id_pattern = args$sample_id_pattern
print(paste0('sample_id_pattern: ', sample_id_pattern))

# csv table with markers information
markers = args$markers
markers = ifelse(markers == 'NaN', NA, file.path(rd, markers))
print(paste0('markers: ', markers))

# Minimum abundance to call an allele
min_abd = as.integer(args$min_abd)
print(paste0('min_abd: ', min_abd))

# Minimum ratio to call minor alleles

min_ratio = as.numeric(args$min_ratio)
print(paste0('min_ratio: ', min_ratio))

# sample_ampl_rate
sample_ampl_rate = as.numeric(args$sample_ampl_rate)
print(paste0('sample_ampl_rate: ', sample_ampl_rate))

# locus_ampl_rate
locus_ampl_rate  = as.numeric(args$locus_ampl_rate)
print(paste0('locus_ampl_rate: ', locus_ampl_rate))

PerformanceReport = as.logical(args$PerformanceReport)
print(paste0('PerformanceReport: ', PerformanceReport))

Drug_Surveillance_Report = as.logical(args$Drug_Surveillance_Report)
print(paste0('Drug_Surveillance_Report: ', Drug_Surveillance_Report))

Variants_of_Interest_Report = as.logical(args$Variants_of_Interest_Report)
print(paste0('Variants_of_Interest_Report: ', Variants_of_Interest_Report))

# Reference gff and fasta files, require to translate DNA cigar formats to aminoacid sequences

ref_gff = args$ref_gff
ref_gff = ifelse(ref_gff == 'NaN', NA, file.path(rd, ref_gff))
print(paste0('ref_gff: ', ref_gff))

ref_fasta = args$ref_fasta
ref_fasta = ifelse(ref_fasta == 'NaN', NA, file.path(rd, ref_fasta))
print(paste0('ref_fasta: ', ref_fasta))

reference_alleles = args$reference_alleles
reference_alleles = ifelse(reference_alleles == 'NaN', NA, file.path(rd, reference_alleles))
print(paste0('reference_alleles: ', reference_alleles))


# gene_names
gene_names = args$gene_names
if(gene_names == 'NaN'){
  gene_names = NA
}else{
  gene_names = strsplit(gene_names, ',')[[1]]
}
print('gene_names:')
print(gene_names)

# gene_ids
gene_ids = args$gene_ids
if(gene_ids == 'NaN'){
  gene_ids = NA
}else{
  gene_ids = strsplit(gene_ids, ',')[[1]]
}
print('gene_ids:')
print(gene_ids)

# metadata
metadata = as.character(args$metadata)
metadata = ifelse(metadata == 'NaN', NA, file.path(wd, metadata))
print(paste0('metadata: ', metadata))

# join metadata by
join_by = as.character(args$join_by)
join_by = ifelse(join_by == 'NaN', NA, join_by)
print(paste0('join_by: ', join_by))

# pop
Variable1 = args$Variable1
Variable1 = ifelse(Variable1 == 'NaN', NA, Variable1)
print(paste0('Variable1: ', Variable1))

# temporal_population
Variable2 = args$Variable2
Variable2 = ifelse(Variable2 == 'NaN', NA, Variable2)
print(paste0('Variable2: ', Variable2))

Longitude = args$Longitude
Longitude = ifelse(Longitude == 'NaN', NA, Longitude)
print(paste0('Longitude: ', Longitude))

Latitude = args$Latitude
Latitude = ifelse(Latitude == 'NaN', NA, Latitude)
print(paste0('Latitude: ', Latitude))

na_var_rm = as.logical(args$na_var_rm)
na_hap_rm = as.logical(args$na_hap_rm)
print(paste0('na_hap_rm: ', na_hap_rm))


drugs = args$drugs
if(drugs == 'NaN'){
  drugs = NA
}else{
  drugs = strsplit(drugs, ',')[[1]]
}
print('drugs: ')
print(drugs)

var_filter  = args$var_filter
if(var_filter == 'null'){
  var_filter = NULL
}else{
  var_filter = gsub('\\.', ' ', var_filter)
  var_filter = strsplit(var_filter, '/')[[1]]
}

print('var_filter: ')
print(var_filter)

include_all_drug_markers = as.logical(args$include_all_drug_markers)
print(paste0('include_all_drug_markers: ', include_all_drug_markers))


# ibd_thres
ibd_thres = args$ibd_thres
ibd_thres = ifelse(ibd_thres == 'NaN', NA, as.numeric(ibd_thres))
print(paste0('ibd_thres: ', ibd_thres))

# parallel
parallel = as.logical(args$parallel)
print(paste0('parallel: ', parallel))

# nTasks

nTasks = as.integer(args$nTasks)
print(paste0('nTasks: ', nTasks))

# Task_id
Task_id = as.integer(args$Task_id)
print(paste0('Task_id: ', Task_id))

# nchunks
nchunks = as.integer(args$nchunks)
print(paste0('nchunks: ', nchunks))

# ibd_step
ibd_step = as.character(args$ibd_step)
print(paste0('ibd_step: ', ibd_step))

# ibd_ncol
ibd_ncol = as.integer(args$ibd_ncol)
print(paste0('ibd_ncol: ', ibd_ncol))

# pop_levels
pop_levels = as.integer(args$pop_levels)
print(paste0('pop_levels: ', pop_levels))


# output pattern
output = args$out
print(paste0('output: ', output))

if(file.exists(file.path(wd, paste0(output, '.xlsx')))){
  cigar_paths = NA
  cigar_files = NA
  ampseq_excelfile = file.path(wd, paste0(output, '.xlsx'))
}


print("All variables checked")

# Check packages and functions----

print("Loading libraies and functions")
source(file.path(fd,'amplseq_required_libraries.R'))
source(file.path(fd,'amplseq_functions.R'))

# Upload Cigar object and convert to ampseq object----

print("starting to upload genotypic data")

if(!is.na(cigar_paths)|!is.na(cigar_files)){
  if(!is.na(cigar_paths)){
    
    print("Uploading genotipic data in cigar format from multiple paths")
    cigar_object = read_cigar_tables(paths = cigar_paths, sample_id_pattern = sample_id_pattern)
  }else if(!is.na(cigar_files)){
    
    if(as.logical(file.info(cigar_files)['isdir'])){
      print("Uploading genotipic data in cigar format from multiple csv files from a single path")
      list_of_cigar_files = file.path(cigar_files,list.files(cigar_files))
      cigar_object = read_cigar_tables(files = list_of_cigar_files, sample_id_pattern = sample_id_pattern)
    }else{
      print("Uploading genotipic data in cigar format from a single file")
      cigar_object = read_cigar_tables(files = cigar_files, sample_id_pattern = sample_id_pattern)
    }
  }
  
  markers = read.csv(markers)
  
  if(PerformanceReport){
    ampseq_object = cigar2ampseq(cigar_object, markers = markers, min_abd = min_abd, min_ratio = min_ratio, remove_controls = F)
 
  }else{
    ampseq_object = cigar2ampseq(cigar_object, markers = markers, min_abd = min_abd, min_ratio = min_ratio, remove_controls = T)
  }
  
  }else if(!is.na(ampseq_jsonfile)){
    print("Uploading genotipic data in ampseq format from json file")
  # In development
  
  }else if(!is.na(ampseq_excelfile)){
    print("Uploading genotipic data in ampseq format from excel file")
  ampseq_object = read_ampseq(file = ampseq_excelfile, format = 'excel')
  
}



# Cleaning, filtering and adding metadata to the  ampseq_object----

if(PerformanceReport){
  
  ampseq_object_filtered = ampseq_object
  
  # add metadata 
  
  if(!is.na(metadata)){
    # Upload metadata from an external source
    metadata = read.csv(metadata)
    
    # Merge the external metadata with our ampseq_object
    ampseq_object_filtered@metadata = left_join(ampseq_object_filtered@metadata,
                                                metadata,
                                                by = join_by)
    }
  
  # Controls
  
  # Coverage
  
  
  # if data was imported as cigar tables
  if(!is.na(cigar_paths)|!is.na(cigar_files)){
    ampseq_object_filtered = filter_samples(ampseq_object, v = ampseq_object@metadata$typeofSamp == 'Samples')
  }
  
  # if locus_ampl_rate was provided
  if(!is.na(locus_ampl_rate)){
    ampseq_object_filtered = locus_amplification_rate(ampseq_object_filtered, threshold = locus_ampl_rate)
  }
  
  # if sample_ampl_rate was provided
  if(!is.na(sample_ampl_rate)){
    ampseq_object_filtered = sample_amplification_rate(ampseq_object_filtered, threshold = sample_ampl_rate)
  }
  
  # if data was imported as cigar tables and no additional steps were provided
  if(!is.na(cigar_paths)|!is.na(cigar_files) & is.na(ibd_thres)){
    
    # save data in excel format
    # write_ampseq(ampseq_object = ampseq_object_filtered,
    #              format = 'excel',
    #              name = file.path(wd, paste0(output, '.xlsx')))
  }
  
  # save performance_report.RData
  
  
  }else{
  
  ampseq_object_filtered = ampseq_object
  
  # if locus_ampl_rate was provided
  if(!is.na(locus_ampl_rate) & "loci_performance" %in% slotNames(ampseq_object_filtered)){
    ampseq_object_filtered = locus_amplification_rate(ampseq_object_filtered, threshold = locus_ampl_rate)
  }
  
  # if sample_ampl_rate was provided
  if(!is.na(sample_ampl_rate) & "sample_ampl_rate" %in% colnames(ampseq_object_filtered@metadata)){
    ampseq_object_filtered = sample_amplification_rate(ampseq_object_filtered, threshold = sample_ampl_rate)
  }
  
  # add metadata 
  
  if(!is.na(metadata)){
    # Upload metadata from an external source
    print("Uploading metadata")
    metadata = read.csv(metadata)
    
    shared_variables = names(ampseq_object_filtered@metadata)[names(ampseq_object_filtered@metadata) %in% names(metadata)]
    shared_variables = shared_variables[!(shared_variables %in% c('Sample_id', join_by))]
    
    if(length(shared_variables) > 0){
      print(paste0('The following variables will be overwritten: ', paste(shared_variables, collapse = ', ')))
      
      for(variable in shared_variables){
        
        ampseq_object_filtered@metadata[[variable]] = NULL
        
      }
      
    }
    
    print("Adding metadata to ampseq object")
    # Merge the external metadata with our ampseq_object
    ampseq_object_filtered@metadata = left_join(ampseq_object_filtered@metadata,
                                                metadata,
                                                by = join_by)
  }
  
  # if data was imported as cigar tables and no additional steps were provided
  if(!is.na(cigar_paths)|!is.na(cigar_files) & is.na(ibd_thres) & !Drug_Surveillance_Report & !Variants_of_Interest_Report & !file.exists(file.path(wd, paste0(output, '.xlsx')))){
    
    # save data in excel format
    # write_ampseq(ampseq_object = ampseq_object_filtered,
    #              format = 'excel',
    #              name = file.path(wd, paste0(output, '.xlsx')))
  }
  
}


# Drug Resistance surveillance ----

if(Drug_Surveillance_Report){
  
  print('Starting drug surveillance report')
  
  variables = 'Sample_id'
  
  if(!is.na(Variable1)){
    variables = c(variables, Variable1)
  }
  
  if(!is.na(Variable2)){
    variables = c(variables, Variable2)
  }
  
  
  if(include_all_drug_markers){
    
    print('including drug markers with low amplification rate')
    ampseq_drug = ampseq_object_filtered
    
    ampseq_drug@gt = cbind(ampseq_drug@gt,
                                       ampseq_drug@discarded_loci$gt[rownames(ampseq_drug@discarded_loci$gt) %in%
                                                                                   rownames(ampseq_drug@gt),
                                                                                 grepl(paste0(gene_names, collapse = '|'),colnames(ampseq_drug@discarded_loci$gt))]
    )
    
    ampseq_drug@markers = rbind(ampseq_drug@markers,
                                ampseq_drug@discarded_loci$markers[grepl(paste0(gene_names, collapse = '|'),ampseq_drug@discarded_loci$markers$amplicon),])
    
  }else{
    
    ampseq_drug = ampseq_object_filtered
    
  }
  
  
  
  drug_resistant_haplotypes_plot = drug_resistant_haplotypes(ampseq_object = ampseq_drug,
                                                             reference_alleles = reference_alleles,
                                                             gene_names = gene_names,
                                                             gene_ids = gene_ids,
                                                             gff_file = ref_gff,
                                                             fasta_file = ref_fasta,
                                                             variables = variables,
                                                             Longitude = Longitude,
                                                             Latitude = Latitude,
                                                             na.var.rm = na_var_rm,
                                                             na.hap.rm = na_hap_rm,
                                                             drugs = drugs,
                                                             filters = var_filter)
  
  print('Generation of plots and tables for DRS report done')
  
  imagename = file.path(wd,paste0(output, '_DRS_Report.RData'))
  
  save(file = imagename, list = c('drug_resistant_haplotypes_plot'))
  
  system(paste0('cp ', file.path(fd, 'MHap_Analysis_DRS_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_DRS_Report.Rmd'))))
  
  # Assign variables based on command-line arguments
  render(file.path(wd, paste0(output, '_DRS_Report.Rmd')), params = list(
    RData_image = imagename),
  output_dir = wd)

  
  
  print("Leaving render script")
  
}

# Variants of interest ----

if(Variants_of_Interest_Report){
  
  variables = 'Sample_id'
  
  if(!is.na(Variable1)){
    variables = c(variables, Variable1)
  }
  
  if(!is.na(Variable2)){
    variables = c(variables, Variable2)
  }
  
  haplotypes_respect_to_reference_plot = haplotypes_respect_to_reference(ampseq_object_filtered,
                                                             gene_names = gene_names,
                                                             gene_ids = gene_ids,
                                                             gff_file = ref_gff,
                                                             fasta_file = ref_fasta,
                                                             variables = variables,
                                                             plot_haplo_freq = TRUE,
                                                             na.var.rm = na_var_rm,
                                                             filters = var_filter)
  
  imagename = file.path(wd,paste0(output, '_VoI_Report.RData'))
  
  save(file = imagename, list = c('haplotypes_respect_to_reference_plot'))
  
}


# Genetic Variants of interest---

# Genetic Relatedness ----

if(!is.na(ibd_thres)){
  
  # call hmmIBD and PCA functions from Rcpp
  sourceCpp(file.path(fd,'hmmloglikelihood.cpp'))
  
  
  if(nTasks > 1){
    
    if(ibd_step == 'pairwise'){
      
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
                                     pairwise_hmmIBD(ampseq_object_filtered, parallel = parallel, w = w, n = nChunks))
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
      
      ## Genetic connectivity----
      
      plot_relatedness_distribution_between = plot_relatedness_distribution(
        pairwise_relatedness = pairwise_relatedness,
        metadata = ampseq_object_filtered@metadata,
        Population = Variable1,
        fill_color = rep('gray50', 10),
        type_pop_comparison = 'between',
        ncol = ncol,
        pop_levels = pop_levels
      )
      
      
      plot_frac_highly_related_between = plot_frac_highly_related(
        pairwise_relatedness = pairwise_relatedness,
        metadata = ampseq_object_filtered@metadata,
        Population = Variable1,
        fill_color = rep('gray50', 10),
        threshold = 0.99,
        type_pop_comparison = 'between',
        pop_levels = pop_levels)
      
      imagename = paste0(output, '_ampseq.RData')
      
      save(file = imagename, list = c('ampseq_object_filtered',
                                      'pairwise_relatedness',
                                      'plot_relatedness_distribution_between', 
                                      'plot_frac_highly_related_between'))
      
      ## Population subdivision----
      
      # sourceCpp(file.path(fd,'Rcpp_functions.cpp'))
      #
      # evectors_IBD = IBD_evectors(ampseq_object = ampseq_object_filtered,
      #                             relatedness_table = pairwise_relatedness,
      #                             k = length(unique(c(pairwise_relatedness$Yi, pairwise_relatedness$Yi))),
      #                             Pop = Variable1, q = 2)
      # 
      # evectors_IBD$eigenvector %>% ggplot(aes(x = PC1, y = PC2, color = Population))+
      #   geom_point(alpha = .7, size = 2) +
      #   stat_ellipse(level = .6)+
      #   scale_color_manual(values = c('firebrick3', 'dodgerblue3', 'gold3'))+
      #   theme_bw()+
      #   labs(title = 'PvGTSeq',
      #        x = paste0('1st PC (', round(evectors_IBD$contrib[1],1), '%)'),
      #        y = paste0('2nd PC (', round(evectors_IBD$contrib[2],1), '%)'),
      #        color = 'Countries')
      
      
      ## Relatedness and transmission----
      
      
    }
    
  }else{
    
    pairwise_relatedness = NULL
    
    for(w in 1:nChunks){
      start = Sys.time()
      pairwise_relatedness = rbind(pairwise_relatedness,
                                   pairwise_hmmIBD(ampseq_object_filtered, parallel = parallel, w = w, n = nChunks))
      time_diff = Sys.time() - start
      
      print(paste0('step ', w, ' done in ', time_diff, ' secs'))
      
    }
    
    
    ## Genetic connectivity----
    
    plot_relatedness_distribution_between = plot_relatedness_distribution(
      pairwise_relatedness = pairwise_relatedness,
      metadata = ampseq_object_filtered@metadata,
      Population = geographic_population,
      fill_color = rep('gray50', 10),
      type_pop_comparison = 'between',
      ncol = ncol,
      pop_levels = pop_levels
    )
    
    
    plot_frac_highly_related_between = plot_frac_highly_related(
      pairwise_relatedness = pairwise_relatedness,
      metadata = ampseq_object_filtered@metadata,
      Population = geographic_population,
      fill_color = rep('gray50', 10),
      threshold = 0.99,
      type_pop_comparison = 'between',
      pop_levels = pop_levels)
    
    
    
    ## Population subdivision----
    
    # evectors_IBD = IBD_evectors(ampseq_object = ampseq_object_filtered,
    #                             relatedness_table = pairwise_relatedness,
    #                             k = length(unique(c(pairwise_relatedness$Yi, pairwise_relatedness$Yi))),
    #                             Pop = geographic_population, q = 2)
    # 
    # evectors_IBD$eigenvector %>% ggplot(aes(x = PC1, y = PC2, color = Population))+
    #   geom_point(alpha = .7, size = 2) +
    #   stat_ellipse(level = .6)+
    #   scale_color_manual(values = c('firebrick3', 'dodgerblue3', 'gold3'))+
    #   theme_bw()+
    #   labs(title = 'PvGTSeq',
    #        x = paste0('1st PC (', round(evectors_IBD$contrib[1],1), '%)'),
    #        y = paste0('2nd PC (', round(evectors_IBD$contrib[2],1), '%)'),
    #        color = 'Countries')
    
    
    ## Relatedness and transmission----
    
  }
  

}


# COI----

# Genetic diversity----



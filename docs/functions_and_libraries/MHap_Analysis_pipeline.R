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

parser$add_argument("-off_target_formula", "--off_target_formula", default = "dVSITES_ij>=0.3",
                    help="max density of variant sites per amplicon")

parser$add_argument("-flanking_INDEL_formula", "--flanking_INDEL_formula", default = "flanking_INDEL==TRUE&h_ij>=0.66",
                    help="Boolean for masking INDELs in flanking areas")

parser$add_argument("-PCR_errors_formula", "--PCR_errors_formula", default = "h_ij>=0.66&h_ijminor>=0.66",
                    help="max density of variant sites per amplicon")

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

parser$add_argument("-pairwise_relatedness_table", "--pairwise_relatedness_table", default = 'NaN',
                    help="string with the file name of the pairwise_relatedness_table")

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

parser$add_argument("-poly_quantile", "--poly_quantile", default = 'NaN',
                    help="Quantile to define polyclonal samples")

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




# off_target_formula filter
off_target_formula = as.character(args$off_target_formula)
off_target_formula = gsub('"',"",off_target_formula)

off_target_formula = gsub('&'," & ",off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('\\|'," \\| ",off_target_formula, ignore.case = TRUE)

if(grepl("\\w>\\d",off_target_formula)){
  patterns = str_extract_all(off_target_formula, "\\w>\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('>',' > ',pattern)
    off_target_formula = gsub(pattern,
                        replacement,
                        off_target_formula, ignore.case = TRUE)
  }
  
}

if(grepl("\\w<\\d",off_target_formula)){
  patterns = str_extract_all(off_target_formula, "\\w<\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('<',' < ',pattern)
    off_target_formula = gsub(pattern,
                        replacement,
                        off_target_formula, ignore.case = TRUE)
  }
  
}

off_target_formula = gsub('>=', " >= ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('<=', " <= ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('==', " == ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('!=', " != ", off_target_formula, ignore.case = TRUE)

off_target_formula = gsub('\\+', " \\+ ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('-', " - ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('\\*', " \\* ", off_target_formula, ignore.case = TRUE)
off_target_formula = gsub('/', " / ", off_target_formula, ignore.case = TRUE)


print(paste0('off_target_formula: ', off_target_formula))


# flanking_INDEL_formula filter
flanking_INDEL_formula = as.character(args$flanking_INDEL_formula)
flanking_INDEL_formula = gsub('"',"",flanking_INDEL_formula)

flanking_INDEL_formula = gsub('&'," & ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('\\|'," \\| ",flanking_INDEL_formula, ignore.case = TRUE)

if(grepl("\\w>\\d",flanking_INDEL_formula)){
  patterns = str_extract_all(flanking_INDEL_formula, "\\w>\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('>',' > ',pattern)
    flanking_INDEL_formula = gsub(pattern,
                        replacement,
                        flanking_INDEL_formula, ignore.case = TRUE)
  }
  
}

if(grepl("\\w<\\d",flanking_INDEL_formula)){
  patterns = str_extract_all(flanking_INDEL_formula, "\\w<\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('<',' < ',pattern)
    flanking_INDEL_formula = gsub(pattern,
                        replacement,
                        flanking_INDEL_formula, ignore.case = TRUE)
  }
  
}

flanking_INDEL_formula = gsub('>='," >= ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('<='," <= ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('=='," == ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('!='," != ",flanking_INDEL_formula, ignore.case = TRUE)

flanking_INDEL_formula = gsub('\\+'," \\+ ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('-'," - ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('\\*'," \\* ",flanking_INDEL_formula, ignore.case = TRUE)
flanking_INDEL_formula = gsub('/'," / ",flanking_INDEL_formula, ignore.case = TRUE)


print(paste0('flanking_INDEL_formula: ', flanking_INDEL_formula))


# PCR_errors_formula filter
PCR_errors_formula = as.character(args$PCR_errors_formula)
PCR_errors_formula = gsub('"',"",PCR_errors_formula)

PCR_errors_formula = gsub('&'," & ",PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('\\|'," \\| ",PCR_errors_formula, ignore.case = TRUE)

if(grepl("\\w>\\d",PCR_errors_formula)){
  patterns = str_extract_all(PCR_errors_formula, "\\w>\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('>',' > ',pattern)
    PCR_errors_formula = gsub(pattern,
                        replacement,
                        PCR_errors_formula, ignore.case = TRUE)
  }
  
}

if(grepl("\\w<\\d",PCR_errors_formula)){
  patterns = str_extract_all(PCR_errors_formula, "\\w<\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('<',' < ',pattern)
    PCR_errors_formula = gsub(pattern,
                        replacement,
                        PCR_errors_formula, ignore.case = TRUE)
  }
  
}

PCR_errors_formula = gsub('>='," >= ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('<='," <= ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('=='," == ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('!='," != ", PCR_errors_formula, ignore.case = TRUE)

PCR_errors_formula = gsub('\\+'," \\+ ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('-'," - ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('\\*'," \\* ", PCR_errors_formula, ignore.case = TRUE)
PCR_errors_formula = gsub('/'," / ", PCR_errors_formula, ignore.case = TRUE)


print(paste0('PCR_errors_formula: ', PCR_errors_formula))

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
metadata_file = as.character(args$metadata)
metadata_file = ifelse(metadata_file == 'NaN', NA, file.path(wd, metadata_file))
print(paste0('metadata: ', metadata_file))

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


# pairwise_relatedness_table
pairwise_relatedness_table = as.character(args$pairwise_relatedness_table)
pairwise_relatedness_table = ifelse(pairwise_relatedness_table == 'NaN', NA, file.path(wd, pairwise_relatedness_table))
print(paste0('pairwise_relatedness_table: ', pairwise_relatedness_table))

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

# poly_quantile

poly_quantile = args$poly_quantile
poly_quantile = ifelse(poly_quantile == 'NaN', NA, as.numeric(poly_quantile))
print(paste0('poly_quantile: ', poly_quantile))


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
    ampseq_object_abd1 = cigar2ampseq(cigar_object, markers = markers, min_abd = 1, min_ratio = min_ratio, remove_controls = F)
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
  
  # add metadata 
  
  if(!is.na(metadata_file)){
    # Upload metadata from an external source
    external_metadata = read.csv(metadata_file)
    
    # Merge the external metadata with our ampseq_object
    ampseq_object_abd1@metadata = left_join(ampseq_object_abd1@metadata,
                                            external_metadata,
                                                by = join_by)
    }
  
  ## Sample performance by different coverage----
  ### Overall sample performance by different coverage----
  
  ReadDepth_coverage = get_ReadDepth_coverage(ampseq_object_abd1, variable = NULL)
  
  sample_performance = ReadDepth_coverage$plot_read_depth_heatmap$data %>%
    mutate(Read_depth = case_when(
      is.na(Read_depth) ~ 0,
      !is.na(Read_depth) ~ Read_depth
    )) %>%
    summarise(amplified_amplicons1 = sum(Read_depth >= 1)/nrow(ampseq_object_abd1@markers),
              amplified_amplicons5 = sum(Read_depth >= 5)/nrow(ampseq_object_abd1@markers),
              amplified_amplicons10 = sum(Read_depth >= 10)/nrow(ampseq_object_abd1@markers),
              amplified_amplicons20 = sum(Read_depth >= 20)/nrow(ampseq_object_abd1@markers),
              amplified_amplicons50 = sum(Read_depth >= 50)/nrow(ampseq_object_abd1@markers),
              amplified_amplicons100 = sum(Read_depth >= 100)/nrow(ampseq_object_abd1@markers),
              .by = Sample_id) %>%
    pivot_longer(cols = starts_with('amplified_amplicons'), values_to = 'AmpRate', names_to = 'Threshold') %>%
    mutate(Threshold = as.integer(gsub('amplified_amplicons', '', Threshold)))
  
  if(!(min_abd %in% c(1,5,10,20,50,100))){
    
    sample_performance = rbind(sample_performance,
                               ReadDepth_coverage$plot_read_depth_heatmap$data %>%
                                 mutate(Read_depth = case_when(
                                   is.na(Read_depth) ~ 0,
                                   !is.na(Read_depth) ~ Read_depth
                                 )) %>%
                                 summarise(AmpRate = sum(Read_depth >= min_abd)/nrow(ampseq_object_abd1@markers),
                                           .by = Sample_id) %>%
                                 mutate(Threshold = min_abd) %>%
                                 select(Sample_id, Threshold, AmpRate))
    
  }
  
  plot_precentage_of_samples_over_min_abd = sample_performance %>%
    summarise(AmpRate5 = round(100*sum(AmpRate >= .05)/n(), 1),
              AmpRate10 = round(100*sum(AmpRate >= .10)/n(), 1),
              AmpRate15 = round(100*sum(AmpRate >= .15)/n(), 1),
              AmpRate20 = round(100*sum(AmpRate >= .20)/n(), 1),
              AmpRate25 = round(100*sum(AmpRate >= .25)/n(), 1),
              AmpRate30 = round(100*sum(AmpRate >= .30)/n(), 1),
              AmpRate35 = round(100*sum(AmpRate >= .35)/n(), 1),
              AmpRate40 = round(100*sum(AmpRate >= .40)/n(), 1),
              AmpRate45 = round(100*sum(AmpRate >= .45)/n(), 1),
              AmpRate50 = round(100*sum(AmpRate >= .50)/n(), 1),
              AmpRate55 = round(100*sum(AmpRate >= .55)/n(), 1),
              AmpRate60 = round(100*sum(AmpRate >= .60)/n(), 1),
              AmpRate65 = round(100*sum(AmpRate >= .65)/n(), 1),
              AmpRate70 = round(100*sum(AmpRate >= .70)/n(), 1),
              AmpRate75 = round(100*sum(AmpRate >= .75)/n(), 1),
              AmpRate80 = round(100*sum(AmpRate >= .80)/n(), 1),
              AmpRate85 = round(100*sum(AmpRate >= .85)/n(), 1),
              AmpRate90 = round(100*sum(AmpRate >= .90)/n(), 1),
              AmpRate95 = round(100*sum(AmpRate >= .95)/n(), 1),
              AmpRate100 = round(100*sum(AmpRate >= 1)/n(), 1),
              .by = c(Threshold)
    ) %>%
    pivot_longer(cols = paste0('AmpRate', seq(5, 100, 5)),
                 values_to = 'Percentage',
                 names_to = 'AmpRate') %>%
    mutate(AmpRate = as.numeric(gsub('AmpRate','', AmpRate)))%>%
    ggplot(aes(x = AmpRate, y = Percentage, color = as.factor(Threshold), group = as.factor(Threshold))) +
    geom_line() +
    geom_vline(xintercept = 100*sample_ampl_rate, linetype = 2) +
    theme_bw() +
    labs(x = '% of amplified loci (amplification rate)', y = '% of Samples', color = 'Min Coverage')
  
  ### Overall sample performance by different coverage per run----
  
  ReadDepth_coverage = get_ReadDepth_coverage(ampseq_object_abd1, variable = 'run')
  
  sample_performance = ReadDepth_coverage$plot_read_depth_heatmap$data %>%
    mutate(Read_depth = case_when(
      is.na(Read_depth) ~ 0,
      !is.na(Read_depth) ~ Read_depth
    )) %>%
    summarise(amplified_amplicons1 = sum(Read_depth >= 1)/nrow(ampseq_object_abd1@markers),
              amplified_amplicons5 = sum(Read_depth >= 5)/nrow(ampseq_object_abd1@markers),
              amplified_amplicons10 = sum(Read_depth >= 10)/nrow(ampseq_object_abd1@markers),
              amplified_amplicons20 = sum(Read_depth >= 20)/nrow(ampseq_object_abd1@markers),
              amplified_amplicons50 = sum(Read_depth >= 50)/nrow(ampseq_object_abd1@markers),
              amplified_amplicons100 = sum(Read_depth >= 100)/nrow(ampseq_object_abd1@markers),
              run = unique(var),
              .by = Sample_id) %>%
    pivot_longer(cols = starts_with('amplified_amplicons'), values_to = 'AmpRate', names_to = 'Threshold') %>%
    mutate(Threshold = as.integer(gsub('amplified_amplicons', '', Threshold)))
  
  if(!(min_abd %in% c(1,5,10,20,50,100))){
    
    sample_performance = rbind(sample_performance,
                               ReadDepth_coverage$plot_read_depth_heatmap$data %>%
                                 mutate(Read_depth = case_when(
                                   is.na(Read_depth) ~ 0,
                                   !is.na(Read_depth) ~ Read_depth
                                 )) %>%
                                 summarise(AmpRate = sum(Read_depth >= min_abd)/nrow(ampseq_object_abd1@markers),
                                           run = unique(var),
                                           .by = Sample_id) %>%
                                 mutate(Threshold = min_abd) %>%
                                 select(Sample_id, run, Threshold, AmpRate))
    
  }
  
  plot_precentage_of_samples_over_min_abd_byRun = sample_performance %>%
    summarise(AmpRate5 = round(100*sum(AmpRate >= .05)/n(), 1),
              AmpRate10 = round(100*sum(AmpRate >= .10)/n(), 1),
              AmpRate15 = round(100*sum(AmpRate >= .15)/n(), 1),
              AmpRate20 = round(100*sum(AmpRate >= .20)/n(), 1),
              AmpRate25 = round(100*sum(AmpRate >= .25)/n(), 1),
              AmpRate30 = round(100*sum(AmpRate >= .30)/n(), 1),
              AmpRate35 = round(100*sum(AmpRate >= .35)/n(), 1),
              AmpRate40 = round(100*sum(AmpRate >= .40)/n(), 1),
              AmpRate45 = round(100*sum(AmpRate >= .45)/n(), 1),
              AmpRate50 = round(100*sum(AmpRate >= .50)/n(), 1),
              AmpRate55 = round(100*sum(AmpRate >= .55)/n(), 1),
              AmpRate60 = round(100*sum(AmpRate >= .60)/n(), 1),
              AmpRate65 = round(100*sum(AmpRate >= .65)/n(), 1),
              AmpRate70 = round(100*sum(AmpRate >= .70)/n(), 1),
              AmpRate75 = round(100*sum(AmpRate >= .75)/n(), 1),
              AmpRate80 = round(100*sum(AmpRate >= .80)/n(), 1),
              AmpRate85 = round(100*sum(AmpRate >= .85)/n(), 1),
              AmpRate90 = round(100*sum(AmpRate >= .90)/n(), 1),
              AmpRate95 = round(100*sum(AmpRate >= .95)/n(), 1),
              AmpRate100 = round(100*sum(AmpRate >= 1)/n(), 1),
              .by = c(Threshold, run)
    ) %>%
    pivot_longer(cols = paste0('AmpRate', seq(5, 100, 5)),
                 values_to = 'Percentage',
                 names_to = 'AmpRate') %>%
    mutate(AmpRate = as.numeric(gsub('AmpRate','', AmpRate)))%>%
    ggplot(aes(x = AmpRate, y = Percentage, color = as.factor(Threshold), group = as.factor(Threshold))) +
    geom_line() +
    geom_vline(xintercept = 100*sample_ampl_rate, linetype = 2) +
    facet_wrap(run~., ncol = 3)+
    theme_bw() +
    labs(x = '% of amplified loci (amplification rate)', y = '% of Samples', color = 'Min Coverage')
  
  
  ### Overall sample performance by different coverage per Variable1----
  
  if(!is.na(Variable1)){
    
    ReadDepth_coverage = get_ReadDepth_coverage(ampseq_object_abd1, variable = Variable1)
    
    sample_performance = ReadDepth_coverage$plot_read_depth_heatmap$data %>%
      filter(!is.na(var))%>%
      mutate(Read_depth = case_when(
        is.na(Read_depth) ~ 0,
        !is.na(Read_depth) ~ Read_depth
      )) %>%
      summarise(amplified_amplicons1 = sum(Read_depth >= 1)/nrow(ampseq_object_abd1@markers),
                amplified_amplicons5 = sum(Read_depth >= 5)/nrow(ampseq_object_abd1@markers),
                amplified_amplicons10 = sum(Read_depth >= 10)/nrow(ampseq_object_abd1@markers),
                amplified_amplicons20 = sum(Read_depth >= 20)/nrow(ampseq_object_abd1@markers),
                amplified_amplicons50 = sum(Read_depth >= 50)/nrow(ampseq_object_abd1@markers),
                amplified_amplicons100 = sum(Read_depth >= 100)/nrow(ampseq_object_abd1@markers),
                Variable1 = unique(var),
                .by = Sample_id) %>%
      pivot_longer(cols = starts_with('amplified_amplicons'), values_to = 'AmpRate', names_to = 'Threshold') %>%
      mutate(Threshold = as.integer(gsub('amplified_amplicons', '', Threshold)))
    
    if(!(min_abd %in% c(1,5,10,20,50,100))){
      
      sample_performance = rbind(sample_performance,
                                 ReadDepth_coverage$plot_read_depth_heatmap$data %>%
                                   filter(!is.na(var))%>%
                                   mutate(Read_depth = case_when(
                                     is.na(Read_depth) ~ 0,
                                     !is.na(Read_depth) ~ Read_depth
                                   )) %>%
                                   summarise(AmpRate = sum(Read_depth >= min_abd)/nrow(ampseq_object_abd1@markers),
                                             Variable1 = unique(var),
                                             .by = Sample_id) %>%
                                   mutate(Threshold = min_abd) %>%
                                   select(Sample_id, Variable1, Threshold, AmpRate))
      
    }
    
    plot_precentage_of_samples_over_min_abd_byVariable1 = sample_performance %>%
      summarise(AmpRate5 = round(100*sum(AmpRate >= .05)/n(), 1),
                AmpRate10 = round(100*sum(AmpRate >= .10)/n(), 1),
                AmpRate15 = round(100*sum(AmpRate >= .15)/n(), 1),
                AmpRate20 = round(100*sum(AmpRate >= .20)/n(), 1),
                AmpRate25 = round(100*sum(AmpRate >= .25)/n(), 1),
                AmpRate30 = round(100*sum(AmpRate >= .30)/n(), 1),
                AmpRate35 = round(100*sum(AmpRate >= .35)/n(), 1),
                AmpRate40 = round(100*sum(AmpRate >= .40)/n(), 1),
                AmpRate45 = round(100*sum(AmpRate >= .45)/n(), 1),
                AmpRate50 = round(100*sum(AmpRate >= .50)/n(), 1),
                AmpRate55 = round(100*sum(AmpRate >= .55)/n(), 1),
                AmpRate60 = round(100*sum(AmpRate >= .60)/n(), 1),
                AmpRate65 = round(100*sum(AmpRate >= .65)/n(), 1),
                AmpRate70 = round(100*sum(AmpRate >= .70)/n(), 1),
                AmpRate75 = round(100*sum(AmpRate >= .75)/n(), 1),
                AmpRate80 = round(100*sum(AmpRate >= .80)/n(), 1),
                AmpRate85 = round(100*sum(AmpRate >= .85)/n(), 1),
                AmpRate90 = round(100*sum(AmpRate >= .90)/n(), 1),
                AmpRate95 = round(100*sum(AmpRate >= .95)/n(), 1),
                AmpRate100 = round(100*sum(AmpRate >= 1)/n(), 1),
                .by = c(Threshold, Variable1)
      ) %>%
      pivot_longer(cols = paste0('AmpRate', seq(5, 100, 5)),
                   values_to = 'Percentage',
                   names_to = 'AmpRate') %>%
      mutate(AmpRate = as.numeric(gsub('AmpRate','', AmpRate)))%>%
      ggplot(aes(x = AmpRate, y = Percentage, color = as.factor(Threshold), group = as.factor(Threshold))) +
      geom_line() +
      geom_vline(xintercept = 100*sample_ampl_rate, linetype = 2) +
      facet_wrap(Variable1~., ncol = 3)+
      theme_bw() +
      labs(x = '% of amplified loci (amplification rate)', y = '% of Samples', color = 'Min Coverage')
    
  }
  
}
  
## Identification of off target-products and PCR artifacts----

### Off-target products----
off_target_stats = frac_ofHet_pAlt_byAllele(ampseq_object)

off_target_formula_check = off_target_formula

print('Check filters to be applied for removing off-target products')

if(grepl("(h_ij|h_ijminor|p_ij|P_ij|H_ij|H_ijminor|nVSITES_ij|dVSITES_ij|nSNPs_ij|dSNPs_ij|nINDELs_ij|dINDELs_ij|flanking_INDEL)(<|>|!|=)+", off_target_formula_check)){
  stop("All mathematical and logical operators must be separated by blank spaces in off_target_formula_check")
}

# modify off_target_formula_check


if(grepl("flanking_INDEL ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "flanking_INDEL (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "flanking_INDEL (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("flanking_INDEL ", "off_target_stats[['flanking_INDEL']] ", off_target_formula_check)
  }else{
    stop("Filter flanking_INDEL is been called but there are spelling issues in this part of the off_target_formula_check")
  }
}


if(grepl("dINDELs_ij ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "dINDELs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "dINDELs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("dINDELs_ij ", "off_target_stats[['dINDELs_ij']] ", off_target_formula_check)
  }else{
    stop("Filter dINDELs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("nINDELs_ij ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "nINDELs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "nINDELs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("nINDELs_ij ", "off_target_stats[['nINDELs_ij']] ", off_target_formula_check)
  }else{
    stop("Filter nINDELs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("dSNPs_ij ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "dSNPs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "dSNPs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("dSNPs_ij ", "off_target_stats[['dSNPs_ij']] ", off_target_formula_check)
  }else{
    stop("Filter dSNPs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("nSNPs_ij ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "nSNPs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "nSNPs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("nSNPs_ij ", "off_target_stats[['nSNPs_ij']] ", off_target_formula_check)
  }else{
    stop("Filter nSNPs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("dVSITES_ij ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "dVSITES_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "dVSITES_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("dVSITES_ij ", "off_target_stats[['dVSITES_ij']] ", off_target_formula_check)
  }else{
    stop("Filter dVSITES_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("nVSITES_ij ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "nVSITES_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "nVSITES_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("nVSITES_ij ", "off_target_stats[['nVSITES_ij']] ", off_target_formula_check)
  }else{
    stop("Filter nVSITES_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("h_ij ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "h_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "h_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("h_ij ", "off_target_stats[['h_ij']] ", off_target_formula_check)
  }else{
    stop("Filter h_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("h_ijminor ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "h_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "h_ijminor (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("h_ijminor ", "off_target_stats[['h_ijminor']] ", off_target_formula_check)
  }else{
    stop("Filter h_ijminor is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("p_ij ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "p_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "p_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("p_ij ", "off_target_stats[['p_ij']] ", off_target_formula_check)
  }else{
    stop("Filter p_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("P_ij ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "P_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "P_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("P_ij ", "off_target_stats[['P_ij']] ", off_target_formula_check)
  }else{
    stop("Filter P_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("H_ij ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "H_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "H_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("H_ij ", "off_target_stats[['H_ij']] ", off_target_formula_check)
  }else{
    stop("Filter H_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("H_ijminor ", off_target_formula_check)){
  
  mask_filter = str_extract(off_target_formula_check, "H_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(off_target_formula_check, "H_ijminor (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    off_target_formula_check = gsub("H_ijminor ", "off_target_stats[['H_ijminor']] ", off_target_formula_check)
  }else{
    stop("Filter H_ijminor is been called but there are spelling issues in this part of the formula")
  }
}

mask_formula_check = str_split(off_target_formula_check, "&|\\|")[[1]]
mask_formula_check  = mask_formula_check[!grepl("off_target_stats", mask_formula_check)]


if(length(mask_formula_check) > 0){
  for(wrong_filter in mask_formula_check){
    print(paste0("Spelling error with filter ", wrong_filter))
  }
  stop("Execution halted, revise mask_filter argument.\nPossible filters are:\nh_ij, h_ijminor, p_ij, P_ij, H_ij, H_ijminor, nVSITES_ij, dVSITES_ij, nSNPs_ij, dSNPs_ij, nINDELs_ij, dINDELs_ij, flanking_INDEL")
}


n_off_target_alleles = off_target_stats[eval(parse(text = off_target_formula_check)),][['Allele']]


if(length(n_off_target_alleles) > 0){
  
  print(paste0(length(n_off_target_alleles), ' allele(s) matches the criteria to define off-target products'))
  
  if(PerformanceReport){
    
    plot_off_target_stats = off_target_stats %>%
      ggplot(aes(x= dVSITES_ij)) + 
      geom_vline(xintercept = as.numeric(gsub('dVSITES_ij (=|!|>|<)+ ',
                                              '',
                                              str_extract(off_target_formula,
                                                          "dVSITES_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"))),
                 linetype = 2) +
      geom_histogram(binwidth = 0.01) + 
      theme_bw()
    }
  
  
  ampseq_object@gt = mask_alt_alleles(ampseq_object, mask_formula = off_target_formula)
}else{
  
  print('No allele matches the criteria to define off-target products')
  
}


### Products with INDELs in their flanking region----

flanking_INDEL_stats = frac_ofHet_pAlt_byAllele(ampseq_object)

flanking_INDEL_formula_check = flanking_INDEL_formula

print('Check filters to be applied for removing off-target products')

if(grepl("(h_ij|h_ijminor|p_ij|P_ij|H_ij|H_ijminor|nVSITES_ij|dVSITES_ij|nSNPs_ij|dSNPs_ij|nINDELs_ij|dINDELs_ij|flanking_INDEL)(<|>|!|=)+", flanking_INDEL_formula_check)){
  stop("All mathematical and logical operators must be separated by blank spaces in flanking_INDEL_formula_check")
}

# modify flanking_INDEL_formula_check


if(grepl("flanking_INDEL ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "flanking_INDEL (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "flanking_INDEL (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("flanking_INDEL ", "flanking_INDEL_stats[['flanking_INDEL']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter flanking_INDEL is been called but there are spelling issues in this part of the flanking_INDEL_formula_check")
  }
}


if(grepl("dINDELs_ij ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "dINDELs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "dINDELs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("dINDELs_ij ", "flanking_INDEL_stats[['dINDELs_ij']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter dINDELs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("nINDELs_ij ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "nINDELs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "nINDELs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("nINDELs_ij ", "flanking_INDEL_stats[['nINDELs_ij']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter nINDELs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("dSNPs_ij ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "dSNPs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "dSNPs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("dSNPs_ij ", "flanking_INDEL_stats[['dSNPs_ij']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter dSNPs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("nSNPs_ij ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "nSNPs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "nSNPs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("nSNPs_ij ", "flanking_INDEL_stats[['nSNPs_ij']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter nSNPs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("dVSITES_ij ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "dVSITES_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "dVSITES_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("dVSITES_ij ", "flanking_INDEL_stats[['dVSITES_ij']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter dVSITES_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("nVSITES_ij ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "nVSITES_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "nVSITES_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("nVSITES_ij ", "flanking_INDEL_stats[['nVSITES_ij']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter nVSITES_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("h_ij ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "h_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "h_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("h_ij ", "flanking_INDEL_stats[['h_ij']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter h_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("h_ijminor ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "h_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "h_ijminor (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("h_ijminor ", "flanking_INDEL_stats[['h_ijminor']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter h_ijminor is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("p_ij ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "p_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "p_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("p_ij ", "flanking_INDEL_stats[['p_ij']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter p_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("P_ij ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "P_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "P_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("P_ij ", "flanking_INDEL_stats[['P_ij']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter P_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("H_ij ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "H_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "H_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("H_ij ", "flanking_INDEL_stats[['H_ij']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter H_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("H_ijminor ", flanking_INDEL_formula_check)){
  
  mask_filter = str_extract(flanking_INDEL_formula_check, "H_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(flanking_INDEL_formula_check, "H_ijminor (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    flanking_INDEL_formula_check = gsub("H_ijminor ", "flanking_INDEL_stats[['H_ijminor']] ", flanking_INDEL_formula_check)
  }else{
    stop("Filter H_ijminor is been called but there are spelling issues in this part of the formula")
  }
}

mask_formula_check = str_split(flanking_INDEL_formula_check, "&|\\|")[[1]]
mask_formula_check  = mask_formula_check[!grepl("flanking_INDEL_stats", mask_formula_check)]


if(length(mask_formula_check) > 0){
  for(wrong_filter in mask_formula_check){
    print(paste0("Spelling error with filter ", wrong_filter))
  }
  stop("Execution halted, revise mask_filter argument.\nPossible filters are:\nh_ij, h_ijminor, p_ij, P_ij, H_ij, H_ijminor, nVSITES_ij, dVSITES_ij, nSNPs_ij, dSNPs_ij, nINDELs_ij, dINDELs_ij, flanking_INDEL")
}

# Count number of alleles to be filtered

n_flanking_INDEL_alleles = flanking_INDEL_stats[eval(parse(text = flanking_INDEL_formula_check)),][['Allele']]


if(length(n_flanking_INDEL_alleles) > 0){
  
  print(paste0(length(n_flanking_INDEL_alleles), ' allele(s) matches the criteria to identify products with flanking INDELs'))
  
  if(PerformanceReport){
    
    h_ij_thres = as.numeric(gsub('h_ij (=|!|>|<)+ ',
                                 '',
                                 str_extract(flanking_INDEL_formula,
                                             "h_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")))
    
    h_ijminor_thres = as.numeric(gsub('h_ijminor (=|!|>|<)+ ',
                                      '',
                                      str_extract(flanking_INDEL_formula,
                                                  "h_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")))
    
    if(!is.na(h_ij_thres) & !is.na(h_ijminor_thres)){
      
      plot_flanking_INDEL_stats = flanking_INDEL_stats %>%
        mutate(h_ijminor_cat = case_when(
          h_ijminor < h_ijminor_thres ~ paste0('h_ijminor < ',h_ijminor_thres),
          h_ijminor >= h_ijminor_thres ~ paste0('h_ijminor >= ',h_ijminor_thres)
        ))%>%
        ggplot(aes(x = p_ij, 
                   y = h_ij,
                   color = h_ijminor))+
        geom_point()+
        geom_hline(yintercept = h_ij_thres,
                   linetype = 2) +
        theme_bw()+
        scale_color_continuous(type = 'viridis')+
        facet_grid(flanking_INDEL~h_ijminor_cat)+
        labs(x = 'Alternative allele frequency (p_ij)',
             y = 'h_ij (H_ij/P_ij)',
             color = 'h_ijminor')
      
    }else if(!is.na(h_ij_thres) & is.na(h_ijminor_thres)){
      
      plot_flanking_INDEL_stats = flanking_INDEL_stats %>%
        ggplot(aes(x = p_ij, 
                   y = h_ij,
                   color = h_ijminor))+
        geom_point()+
        geom_hline(yintercept = h_ij_thres,
                   linetype = 2) +
        theme_bw()+
        scale_color_continuous(type = 'viridis')+
        facet_grid(flanking_INDEL~.)+
        labs(x = 'Alternative allele frequency (p_ij)',
             y = 'h_ij (H_ij/P_ij)',
             color = 'h_ijminor')
      
    }else if(is.na(h_ij_thres) & !is.na(h_ijminor_thres)){
      
      plot_flanking_INDEL_stats = flanking_INDEL_stats %>%
        mutate(h_ijminor_cat = case_when(
          h_ijminor < h_ijminor_thres ~ paste0('h_ijminor < ',h_ijminor_thres),
          h_ijminor >= h_ijminor_thres ~ paste0('h_ijminor >= ',h_ijminor_thres)
        ))%>%
        ggplot(aes(x = p_ij, 
                   y = h_ij,
                   color = h_ijminor))+
        geom_point()+
        theme_bw()+
        scale_color_continuous(type = 'viridis')+
        facet_grid(flanking_INDEL~h_ijminor_cat)+
        labs(x = 'Alternative allele frequency (p_ij)',
             y = 'h_ij (H_ij/P_ij)',
             color = 'h_ijminor')
      
    }else if(is.na(h_ij_thres) & is.na(h_ijminor_thres)){
      
      plot_flanking_INDEL_stats = flanking_INDEL_stats %>%
        ggplot(aes(x = p_ij, 
                   y = h_ij,
                   color = h_ijminor))+
        geom_point()+
        theme_bw()+
        scale_color_continuous(type = 'viridis')+
        facet_grid(flanking_INDEL~.)+
        labs(x = 'Alternative allele frequency (p_ij)',
             y = 'h_ij (H_ij/P_ij)',
             color = 'h_ijminor')
      
    }
    
  }
  
  
  # mask flanking_INDEL
  ampseq_object@gt = mask_alt_alleles(ampseq_object, mask_formula = flanking_INDEL_formula)
}else{
  
  print('No allele matches the criteria to identify products with flanking INDELs')
  
}



### PCR Errors present as heterozygous----

PCR_errors_stats = frac_ofHet_pAlt_byAllele(ampseq_object)

PCR_errors_formula_check = PCR_errors_formula

print('Check filters to be applied for removing off-target products')

if(grepl("(h_ij|h_ijminor|p_ij|P_ij|H_ij|H_ijminor|nVSITES_ij|dVSITES_ij|nSNPs_ij|dSNPs_ij|nINDELs_ij|dINDELs_ij|flanking_INDEL)(<|>|!|=)+", PCR_errors_formula_check)){
  stop("All mathematical and logical operators must be separated by blank spaces in PCR_errors_formula_check")
}

# modify PCR_errors_formula_check


if(grepl("flanking_INDEL ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "flanking_INDEL (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "flanking_INDEL (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("flanking_INDEL ", "PCR_errors_stats[['flanking_INDEL']] ", PCR_errors_formula_check)
  }else{
    stop("Filter flanking_INDEL is been called but there are spelling issues in this part of the PCR_errors_formula_check")
  }
}


if(grepl("dINDELs_ij ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "dINDELs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "dINDELs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("dINDELs_ij ", "PCR_errors_stats[['dINDELs_ij']] ", PCR_errors_formula_check)
  }else{
    stop("Filter dINDELs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("nINDELs_ij ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "nINDELs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "nINDELs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("nINDELs_ij ", "PCR_errors_stats[['nINDELs_ij']] ", PCR_errors_formula_check)
  }else{
    stop("Filter nINDELs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("dSNPs_ij ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "dSNPs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "dSNPs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("dSNPs_ij ", "PCR_errors_stats[['dSNPs_ij']] ", PCR_errors_formula_check)
  }else{
    stop("Filter dSNPs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("nSNPs_ij ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "nSNPs_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "nSNPs_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("nSNPs_ij ", "PCR_errors_stats[['nSNPs_ij']] ", PCR_errors_formula_check)
  }else{
    stop("Filter nSNPs_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("dVSITES_ij ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "dVSITES_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "dVSITES_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("dVSITES_ij ", "PCR_errors_stats[['dVSITES_ij']] ", PCR_errors_formula_check)
  }else{
    stop("Filter dVSITES_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("nVSITES_ij ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "nVSITES_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "nVSITES_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("nVSITES_ij ", "PCR_errors_stats[['nVSITES_ij']] ", PCR_errors_formula_check)
  }else{
    stop("Filter nVSITES_ij is been called but there are spelling issues in this part of the formula")
  }
}


if(grepl("h_ij ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "h_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "h_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("h_ij ", "PCR_errors_stats[['h_ij']] ", PCR_errors_formula_check)
  }else{
    stop("Filter h_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("h_ijminor ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "h_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "h_ijminor (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("h_ijminor ", "PCR_errors_stats[['h_ijminor']] ", PCR_errors_formula_check)
  }else{
    stop("Filter h_ijminor is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("p_ij ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "p_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "p_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("p_ij ", "PCR_errors_stats[['p_ij']] ", PCR_errors_formula_check)
  }else{
    stop("Filter p_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("P_ij ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "P_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "P_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("P_ij ", "PCR_errors_stats[['P_ij']] ", PCR_errors_formula_check)
  }else{
    stop("Filter P_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("H_ij ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "H_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "H_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("H_ij ", "PCR_errors_stats[['H_ij']] ", PCR_errors_formula_check)
  }else{
    stop("Filter H_ij is been called but there are spelling issues in this part of the formula")
  }
}

if(grepl("H_ijminor ", PCR_errors_formula_check)){
  
  mask_filter = str_extract(PCR_errors_formula_check, "H_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")
  
  if(!is.na(mask_filter)){
    print(paste0('Filter ', str_extract(PCR_errors_formula_check, "H_ijminor (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)"), ' will be applied'))
    PCR_errors_formula_check = gsub("H_ijminor ", "PCR_errors_stats[['H_ijminor']] ", PCR_errors_formula_check)
  }else{
    stop("Filter H_ijminor is been called but there are spelling issues in this part of the formula")
  }
}

mask_formula_check = str_split(PCR_errors_formula_check, "&|\\|")[[1]]
mask_formula_check  = mask_formula_check[!grepl("PCR_errors_stats", mask_formula_check)]


if(length(mask_formula_check) > 0){
  for(wrong_filter in mask_formula_check){
    print(paste0("Spelling error with filter ", wrong_filter))
  }
  stop("Execution halted, revise mask_filter argument.\nPossible filters are:\nh_ij, h_ijminor, p_ij, P_ij, H_ij, H_ijminor, nVSITES_ij, dVSITES_ij, nSNPs_ij, dSNPs_ij, nINDELs_ij, dINDELs_ij, flanking_INDEL")
}

# Count number of alleles to be filtered

n_PCR_errors_alleles = PCR_errors_stats[eval(parse(text = PCR_errors_formula_check)),][['Allele']]


if(length(n_PCR_errors_alleles) > 0){
  
  print(paste0(length(n_PCR_errors_alleles), ' allele(s) matches the criteria to identify PCR_errors'))
  
  if(PerformanceReport){
    
    h_ij_thres = as.numeric(gsub('h_ij (=|!|>|<)+ ',
                                 '',
                                 str_extract(PCR_errors_formula,
                                             "h_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")))
    
    h_ijminor_thres = as.numeric(gsub('h_ijminor (=|!|>|<)+ ',
                                      '',
                                      str_extract(PCR_errors_formula,
                                                  "h_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+|TRUE|FALSE)")))
    
    if(!is.na(h_ij_thres) & !is.na(h_ijminor_thres)){
      
      plot_PCR_errors_stats = PCR_errors_stats %>%
        mutate(h_ijminor_cat = case_when(
          h_ijminor < h_ijminor_thres ~ paste0('h_ijminor < ',h_ijminor_thres),
          h_ijminor >= h_ijminor_thres ~ paste0('h_ijminor >= ',h_ijminor_thres)
        ))%>%
        ggplot(aes(x = p_ij, 
                   y = h_ij,
                   color = h_ijminor))+
        geom_point()+
        geom_hline(yintercept = h_ij_thres,
                   linetype = 2) +
        theme_bw()+
        scale_color_continuous(type = 'viridis')+
        facet_grid(.~h_ijminor_cat)+
        labs(x = 'Alternative allele frequency (p_ij)',
             y = 'h_ij (H_ij/P_ij)',
             color = 'h_ijminor')
      
    }else if(!is.na(h_ij_thres) & is.na(h_ijminor_thres)){
      
      plot_PCR_errors_stats = PCR_errors_stats %>%
        ggplot(aes(x = p_ij, 
                   y = h_ij,
                   color = h_ijminor))+
        geom_point()+
        geom_hline(yintercept = h_ij_thres,
                   linetype = 2) +
        theme_bw()+
        scale_color_continuous(type = 'viridis')+
        labs(x = 'Alternative allele frequency (p_ij)',
             y = 'h_ij (H_ij/P_ij)',
             color = 'h_ijminor')
      
    }else if(is.na(h_ij_thres) & !is.na(h_ijminor_thres)){
      
      plot_PCR_errors_stats = PCR_errors_stats %>%
        mutate(h_ijminor_cat = case_when(
          h_ijminor < h_ijminor_thres ~ paste0('h_ijminor < ',h_ijminor_thres),
          h_ijminor >= h_ijminor_thres ~ paste0('h_ijminor >= ',h_ijminor_thres)
        ))%>%
        ggplot(aes(x = p_ij, 
                   y = h_ij,
                   color = h_ijminor))+
        geom_point()+
        theme_bw()+
        scale_color_continuous(type = 'viridis')+
        facet_grid(.~h_ijminor_cat)+
        labs(x = 'Alternative allele frequency (p_ij)',
             y = 'h_ij (H_ij/P_ij)',
             color = 'h_ijminor')
      
    }else if(is.na(h_ij_thres) & is.na(h_ijminor_thres)){
      
      plot_PCR_errors_stats = PCR_errors_stats %>%
        ggplot(aes(x = p_ij, 
                   y = h_ij,
                   color = h_ijminor))+
        geom_point()+
        theme_bw()+
        scale_color_continuous(type = 'viridis')+
        labs(x = 'Alternative allele frequency (p_ij)',
             y = 'h_ij (H_ij/P_ij)',
             color = 'h_ijminor')
      
    }
    
  }
  
  
  
  # mask flanking_INDEL
  gt_masked = mask_alt_alleles(ampseq_object, mask_formula = PCR_errors_formula)
  
  print('gt_masked generated')
  
  ampseq_object@gt = gt_masked
  
}else{
  
  print('No allele matches the criteria to identify PCR_errors')
  
}


## Coverage by sample and amplicon----



if(!is.na(metadata_file)){
  # Merge the external metadata with our ampseq_object
  
  external_metadata = read.csv(metadata_file)
  
  ampseq_object@metadata = left_join(ampseq_object@metadata,
                                     external_metadata,
                                     by = join_by)
  
  print('external metadata added to ampseq object')
}


if(PerformanceReport == TRUE){
  
  print('Measuring coverage per sample per amplicon')
  
  # Read Depth coverage by the whole dataset or by Variable1
  
  if(is.na(Variable1)){
    ReadDepth_coverage = get_ReadDepth_coverage(ampseq_object, variable = NULL)
  }else{
    ReadDepth_coverage = get_ReadDepth_coverage(ampseq_object, variable = Variable1)
  }
  
  # Read Depth coverage by the whole dataset or by Run
  
  ReadDepth_coverage_by_run = get_ReadDepth_coverage(ampseq_object, variable = "run")
  
  # Read Depth coverage by the whole dataset or by Run for controls
  
  ampseq_object_controls = filter_samples(ampseq_object, v = ampseq_object@metadata$typeofSamp == 'Controls')
  
  ReadDepth_coverage_by_run_controls = get_ReadDepth_coverage(ampseq_object_controls, variable = "run")
  
}

# Remove controls

print('Removing Controls for further steps')

ampseq_object = filter_samples(ampseq_object, v = ampseq_object@metadata$typeofSamp == 'Samples')


print("Removing undesired categories based on var_filter")
if(!is.null(var_filter)){
  filters = strsplit(var_filter,';')
  for(temp_filter in 1:length(filters)){
    
    ampseq_object = filter_samples(ampseq_object,
                                   ampseq_object@metadata[[filters[[temp_filter]][1]]] %in% strsplit(filters[[temp_filter]][2],',')[[1]])
  }
}

# if locus_ampl_rate was provided

print('Measuring amplification rate by locus')

if(!is.na(locus_ampl_rate)){
  
  if(!is.na(Variable1)){
    ampseq_object = locus_amplification_rate(ampseq_object, threshold = locus_ampl_rate, strata = Variable1)
  }else{
    ampseq_object = locus_amplification_rate(ampseq_object, threshold = locus_ampl_rate)
  }
  
}else{
  
  if(!is.na(Variable1)){
    ampseq_object = locus_amplification_rate(ampseq_object, threshold = 0.65, strata = Variable1)
  }else{
    ampseq_object = locus_amplification_rate(ampseq_object, threshold = 0.65) 
  }
}

# if sample_ampl_rate was provided

print("measuring amplification rate by sample")

if(!is.na(sample_ampl_rate)){
  if(!is.na(Variable1)){
    ampseq_object = sample_amplification_rate(ampseq_object, threshold = sample_ampl_rate, strata = Variable1)
  }else{
    ampseq_object = sample_amplification_rate(ampseq_object, threshold = sample_ampl_rate)  
  }
  
}else{
  
  if(!is.na(Variable1)){
    ampseq_object = sample_amplification_rate(ampseq_object, threshold = 0.80, strata = Variable1)
  }else{
    ampseq_object = sample_amplification_rate(ampseq_object, threshold = 0.80)  
  }

}

## Printing performance report----

if(PerformanceReport == TRUE){
  
  all_loci_amplification_rate = ampseq_object@plots$all_loci_amplification_rate
  
  samples_amplification_rate = ampseq_object@plots$samples_amplification_rate
  
  cigar_table_unmasked_unfiltered = ampseq_object_abd1@gt 
  cigar_table_masked_filtered = ampseq_object@gt
  cigar_table_controls_masked_filtered = ampseq_object_controls@gt 
  print('Generation of plots and tables for Performance report done')
  
  
  Performance_Report_expected_outputs = c(
    'cigar_table_unmasked_unfiltered',
    'cigar_table_masked_filtered',
    'cigar_table_controls_masked_filtered',
    
    'plot_precentage_of_samples_over_min_abd',
    'plot_precentage_of_samples_over_min_abd_byRun',
    'plot_precentage_of_samples_over_min_abd_byVariable1',
    
    'n_off_target_alleles',
    'off_target_stats',
    'plot_off_target_stats',
    
    'n_flanking_INDEL_alleles',
    'flanking_INDEL_stats',
    'plot_flanking_INDEL_stats',
    
    'n_PCR_errors_alleles',
    'PCR_errors_stats',
    'plot_PCR_errors_stats',
    
    'ReadDepth_coverage',
    'ReadDepth_coverage_by_run',
    'ReadDepth_coverage_by_run_controls',
    'all_loci_amplification_rate',
    'samples_amplification_rate'
    )
  
  Performance_Report_outputs = Performance_Report_expected_outputs[Performance_Report_expected_outputs %in% ls()]
  
  imagename = file.path(wd,paste0(output, '_Performance_Report.RData'))
  
  save(file = imagename, list = Performance_Report_outputs)
  
  system(paste0('cp ', file.path(fd, 'MHap_Analysis_Performance_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_Performance_Report.Rmd'))))
  
  # Assign variables based on command-line arguments
  render(file.path(wd, paste0(output, '_Performance_Report.Rmd')), params = list(
    RData_image = imagename),
    output_dir = wd)
  
  
  
  print("Leaving render script")
  
}


# if data was imported as cigar tables and no additional steps were provided
# if(!is.na(cigar_paths)|!is.na(cigar_files) & is.na(ibd_thres)){
#   
#   # save data in excel format
#   # write_ampseq(ampseq_object = ampseq_object_filtered,
#   #              format = 'excel',
#   #              name = file.path(wd, paste0(output, '.xlsx')))
# }

# save performance_report.RData

# 
# 
# ampseq_object_filtered = ampseq_object
# 
# # if locus_ampl_rate was provided
# if(!is.na(locus_ampl_rate) & "loci_performance" %in% slotNames(ampseq_object_filtered)){
#   ampseq_object_filtered = locus_amplification_rate(ampseq_object_filtered, threshold = locus_ampl_rate)
# }
# 
# # if sample_ampl_rate was provided
# if(!is.na(sample_ampl_rate) & "sample_ampl_rate" %in% colnames(ampseq_object_filtered@metadata)){
#   ampseq_object_filtered = sample_amplification_rate(ampseq_object_filtered, threshold = sample_ampl_rate)
# }
# 
# # add metadata 
# 
# if(!is.na(metadata)){
#   # Upload metadata from an external source
#   print("Uploading metadata")
#   metadata = read.csv(metadata)
#   
#   shared_variables = names(ampseq_object_filtered@metadata)[names(ampseq_object_filtered@metadata) %in% names(metadata)]
#   shared_variables = shared_variables[!(shared_variables %in% c('Sample_id', join_by))]
#   
#   if(length(shared_variables) > 0){
#     print(paste0('The following variables will be overwritten: ', paste(shared_variables, collapse = ', ')))
#     
#     for(variable in shared_variables){
#       
#       ampseq_object_filtered@metadata[[variable]] = NULL
#       
#     }
#     
#   }
#   
#   print("Adding metadata to ampseq object")
#   # Merge the external metadata with our ampseq_object
#   ampseq_object_filtered@metadata = left_join(ampseq_object_filtered@metadata,
#                                               metadata,
#                                               by = join_by)
# }
# 
# # if data was imported as cigar tables and no additional steps were provided
# if(!is.na(cigar_paths)|!is.na(cigar_files) & is.na(ibd_thres) & !Drug_Surveillance_Report & !Variants_of_Interest_Report & !file.exists(file.path(wd, paste0(output, '.xlsx')))){
#   
#   # save data in excel format
#   # write_ampseq(ampseq_object = ampseq_object_filtered,
#   #              format = 'excel',
#   #              name = file.path(wd, paste0(output, '.xlsx')))
# }
  



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
    ampseq_drug = ampseq_object
    
    ampseq_drug@gt = cbind(ampseq_drug@gt,
                                       ampseq_drug@discarded_loci$gt[rownames(ampseq_drug@discarded_loci$gt) %in%
                                                                                   rownames(ampseq_drug@gt),
                                                                                 grepl(paste0(gene_names, collapse = '|'),colnames(ampseq_drug@discarded_loci$gt))]
    )
    
    ampseq_drug@markers = rbind(ampseq_drug@markers,
                                ampseq_drug@discarded_loci$markers[grepl(paste0(gene_names, collapse = '|'),ampseq_drug@discarded_loci$markers$amplicon),])
    
  }else{
    
    ampseq_drug = ampseq_object
    
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
                                                             filters = NULL)
  
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
  
  haplotypes_respect_to_reference_plot = haplotypes_respect_to_reference(ampseq_object,
                                                             gene_names = gene_names,
                                                             gene_ids = gene_ids,
                                                             gff_file = ref_gff,
                                                             fasta_file = ref_fasta,
                                                             variables = variables,
                                                             plot_haplo_freq = TRUE,
                                                             na.var.rm = na_var_rm,
                                                             filters = NULL)
  
  imagename = file.path(wd,paste0(output, '_VoI_Report.RData'))
  
  save(file = imagename, list = c('haplotypes_respect_to_reference_plot'))
  
}


# Genetic Variants of interest---

# Genetic Relatedness ----

if(!is.na(ibd_thres)){
  
  # call hmmIBD and PCA functions from Rcpp
  sourceCpp(file.path(fd,'hmmloglikelihood.cpp'))
  
  if(nTasks > 1 & !is.na(nTasks)){
    
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
                                     pairwise_hmmIBD(ampseq_object, parallel = parallel, w = w, n = nChunks))
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
        metadata = ampseq_object@metadata,
        Population = Variable1,
        fill_color = rep('gray50', 10),
        type_pop_comparison = 'between',
        ncol = 3,
        pop_levels = NULL
      )
      
      
      plot_frac_highly_related_between = plot_frac_highly_related(
        pairwise_relatedness = pairwise_relatedness,
        metadata = ampseq_object@metadata,
        Population = Variable1,
        fill_color = rep('gray50', 10),
        threshold = ibd_thres,
        type_pop_comparison = 'between',
        pop_levels = NULL)
      
      
      if(!is.na(Variable2)){
        plot_frac_highly_related_overtime_between = plot_frac_highly_related_over_time(
          pairwise_relatedness = pairwise_relatedness,
          metadata = ampseq_object@metadata,
          Population = c(Variable1, Variable2),
          fill_color = rep('gray50', 10),
          threshold = ibd_thres,
          type_pop_comparison = 'between',
          ncol = 3,
          pop_levels = NULL)
      }
      
      
      ## Population subdivision----
      
      evectors_IBD = IBD_evectors(ampseq_object = ampseq_object,
                                  relatedness_table = pairwise_relatedness,
                                  k = length(unique(ampseq_object@metadata$Sample_id)),
                                  Pop = Variable1, q = 2)
      
      names(evectors_IBD$eigenvector)[3] = 'Variable1'
      
      IBD_PCA = evectors_IBD$eigenvector %>% ggplot(aes(x = PC1, y = PC2, color = Variable1))+
        geom_point(alpha = .7, size = 2) +
        stat_ellipse(level = .6)+
        scale_color_manual(values = 
                             brewer.pal(n = nlevels(as.factor(ampseq_object@metadata[[Variable1]])), 
                                        name = 'Accent'))+
        theme_bw()+
        labs(x = paste0('1st PC (', round(evectors_IBD$contrib[1],1), '%)'),
             y = paste0('2nd PC (', round(evectors_IBD$contrib[2],1), '%)'),
             color = 'Countries')
      
      
      
      print('Generation of plots and tables for IBD and Connectivity report done')
      
      imagename = file.path(wd,paste0(output, '_IBD_Connectivity_Report.RData'))
      
      IBD_Connectivity_Report_expected_outputs = c('plot_relatedness_distribution_between',
                                                   'plot_frac_highly_related_between',
                                                   'plot_frac_highly_related_overtime_between',
                                                   'evectors_IBD',
                                                   'IBD_PCA',
                                                   'pairwise_relatedness',
                                                   'ibd_thres',
                                                   'ampseq_object',
                                                   'plot_network',
                                                   'create_ampseq',
                                                   'Variable1',
                                                   'Variable2')
      
      IBD_Connectivity_Report_outputs = IBD_Connectivity_Report_expected_outputs[IBD_Connectivity_Report_expected_outputs %in% ls()]
      
      save(file = imagename, list = IBD_Connectivity_Report_outputs)
      
      system(paste0('cp ', file.path(fd, 'MHap_Analysis_IBD_Connectivity_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_IBD_Connectivity_Report.Rmd'))))
      
      # Assign variables based on command-line arguments
      render(file.path(wd, paste0(output, '_IBD_Connectivity_Report.Rmd')), params = list(
        RData_image = imagename),
        output_dir = wd)
      
      print("Leaving render script")
      
      
      ## Relatedness and transmission----
      
      plot_relatedness_distribution_within = plot_relatedness_distribution(
        pairwise_relatedness = pairwise_relatedness,
        metadata = ampseq_object@metadata,
        Population = Variable1,
        fill_color = rep('gray50', 10),
        type_pop_comparison = 'within',
        ncol = 3,
        pop_levels = NULL
      )
      
      
      plot_frac_highly_related_within = plot_frac_highly_related(
        pairwise_relatedness = pairwise_relatedness,
        metadata = ampseq_object@metadata,
        Population = Variable1,
        fill_color = rep('gray50', 10),
        threshold = ibd_thres,
        type_pop_comparison = 'within',
        pop_levels = NULL)
      
      
      if(!is.na(Variable2)){
        plot_frac_highly_related_overtime_within = plot_frac_highly_related_over_time(
          pairwise_relatedness = pairwise_relatedness,
          metadata = ampseq_object@metadata,
          Population = c(Variable1, Variable2),
          fill_color = rep('gray50', 10),
          threshold = ibd_thres,
          type_pop_comparison = 'within',
          ncol = ibd_ncol,
          pop_levels = NULL)
      }
      
      print('Generation of plots and tables for IBD and Transmission report done')
      
      imagename = file.path(wd,paste0(output, '_IBD_Transmission_Report.RData'))
      
      IBD_Transmission_Report_expected_outputs = c('plot_relatedness_distribution_within',
                                                   'plot_frac_highly_related_within',
                                                   'plot_frac_highly_related_overtime_within',
                                                   'Variable2')
      
      IBD_Transmission_Report_outputs = IBD_Transmission_Report_expected_outputs[IBD_Transmission_Report_expected_outputs %in% ls()]
      
      save(file = imagename, list = IBD_Transmission_Report_outputs)
      
      system(paste0('cp ', file.path(fd, 'MHap_Analysis_IBD_Transmission_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_IBD_Transmission_Report.Rmd'))))
      
      # Assign variables based on command-line arguments
      render(file.path(wd, paste0(output, '_IBD_Transmission_Report.Rmd')), params = list(
        RData_image = imagename),
        output_dir = wd)
      
      print("Leaving render script")
      
      
    }
    
  }else if(nTasks == 1 | is.na(nTasks)){
    
    if(is.na(pairwise_relatedness_table)){
      pairwise_relatedness = NULL
      
      for(w in 1:nChunks){
        start = Sys.time()
        pairwise_relatedness = rbind(pairwise_relatedness,
                                     pairwise_hmmIBD(ampseq_object, parallel = parallel, w = w, n = nChunks))
        time_diff = Sys.time() - start
        
        print(paste0('step ', w, ' done in ', time_diff, ' secs'))
        
      }
    }else{
      
      pairwise_relatedness = read.csv(pairwise_relatedness_table)
      
    }
    
    
    ## Genetic connectivity----
    
    plot_relatedness_distribution_between = plot_relatedness_distribution(
      pairwise_relatedness = pairwise_relatedness,
      metadata = ampseq_object@metadata,
      Population = Variable1,
      fill_color = rep('gray50', 10),
      type_pop_comparison = 'between',
      ncol = 3,
      pop_levels = NULL
    )
    
    
    plot_frac_highly_related_between = plot_frac_highly_related(
      pairwise_relatedness = pairwise_relatedness,
      metadata = ampseq_object@metadata,
      Population = Variable1,
      fill_color = rep('gray50', 10),
      threshold = ibd_thres,
      type_pop_comparison = 'between',
      pop_levels = NULL)
    
    
    if(!is.na(Variable2)){
      plot_frac_highly_related_overtime_between = plot_frac_highly_related_over_time(
        pairwise_relatedness = pairwise_relatedness,
        metadata = ampseq_object@metadata,
        Population = c(Variable1, Variable2),
        fill_color = rep('gray50', 10),
        threshold = ibd_thres,
        type_pop_comparison = 'between',
        ncol = 3,
        pop_levels = NULL)
    }
    
    
    ## Population subdivision----
    
    evectors_IBD = IBD_evectors(ampseq_object = ampseq_object,
                                relatedness_table = pairwise_relatedness,
                                k = length(unique(ampseq_object@metadata$Sample_id)),
                                Pop = Variable1, q = 2)
    
    names(evectors_IBD$eigenvector)[3] = 'Variable1'

    IBD_PCA = evectors_IBD$eigenvector %>% ggplot(aes(x = PC1, y = PC2, color = Variable1))+
      geom_point(alpha = .7, size = 2) +
      stat_ellipse(level = .6)+
      scale_color_manual(values = 
                           brewer.pal(n = nlevels(as.factor(ampseq_object@metadata[[Variable1]])), 
                                      name = 'Accent'))+
      theme_bw()+
      labs(x = paste0('1st PC (', round(evectors_IBD$contrib[1],1), '%)'),
           y = paste0('2nd PC (', round(evectors_IBD$contrib[2],1), '%)'),
           color = 'Countries')
    
    
    
    print('Generation of plots and tables for IBD and Connectivity report done')
    
    imagename = file.path(wd,paste0(output, '_IBD_Connectivity_Report.RData'))
    
    IBD_Connectivity_Report_expected_outputs = c('plot_relatedness_distribution_between',
                                                 'plot_frac_highly_related_between',
                                                 'plot_frac_highly_related_overtime_between',
                                                 'evectors_IBD',
                                                 'IBD_PCA',
                                                 'pairwise_relatedness',
                                                 'ibd_thres',
                                                 'ampseq_object',
                                                 'plot_network',
                                                 'create_ampseq',
                                                 'Variable1',
                                                 'Variable2')
    
    IBD_Connectivity_Report_outputs = IBD_Connectivity_Report_expected_outputs[IBD_Connectivity_Report_expected_outputs %in% ls()]
    
    save(file = imagename, list = IBD_Connectivity_Report_outputs)
    
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_IBD_Connectivity_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_IBD_Connectivity_Report.Rmd'))))
    
    # Assign variables based on command-line arguments
    render(file.path(wd, paste0(output, '_IBD_Connectivity_Report.Rmd')), params = list(
      RData_image = imagename),
      output_dir = wd)
    
    print("Leaving render script")
    
    
    ## Relatedness and transmission----
    
    plot_relatedness_distribution_within = plot_relatedness_distribution(
      pairwise_relatedness = pairwise_relatedness,
      metadata = ampseq_object@metadata,
      Population = Variable1,
      fill_color = rep('gray50', 10),
      type_pop_comparison = 'within',
      ncol = 3,
      pop_levels = NULL
    )
    
    
    plot_frac_highly_related_within = plot_frac_highly_related(
      pairwise_relatedness = pairwise_relatedness,
      metadata = ampseq_object@metadata,
      Population = Variable1,
      fill_color = rep('gray50', 10),
      threshold = ibd_thres,
      type_pop_comparison = 'within',
      pop_levels = NULL)
    
    
    if(!is.na(Variable2)){
      plot_frac_highly_related_overtime_within = plot_frac_highly_related_over_time(
        pairwise_relatedness = pairwise_relatedness,
        metadata = ampseq_object@metadata,
        Population = c(Variable1, Variable2),
        fill_color = rep('gray50', 10),
        threshold = ibd_thres,
        type_pop_comparison = 'within',
        ncol = ibd_ncol,
        pop_levels = NULL)
    }
    
    print('Generation of plots and tables for IBD and Transmission report done')
    
    imagename = file.path(wd,paste0(output, '_IBD_Transmission_Report.RData'))
    
    IBD_Transmission_Report_expected_outputs = c('plot_relatedness_distribution_within',
                                                 'plot_frac_highly_related_within',
                                                 'plot_frac_highly_related_overtime_within',
                                                 'Variable2')
    
    IBD_Transmission_Report_outputs = IBD_Transmission_Report_expected_outputs[IBD_Transmission_Report_expected_outputs %in% ls()]
    
    save(file = imagename, list = IBD_Transmission_Report_outputs)
    
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_IBD_Transmission_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_IBD_Transmission_Report.Rmd'))))
    
    # Assign variables based on command-line arguments
    render(file.path(wd, paste0(output, '_IBD_Transmission_Report.Rmd')), params = list(
      RData_image = imagename),
      output_dir = wd)
    
    print("Leaving render script")
    
    
    
  }
  

}


# COI----

if(!is.na(poly_quantile)){
  
  
  if(!is.na(Variable1)){
    
    poly_by_Var1 = get_polygenomic(ampseq_object = ampseq_object, 
                                   strata = Variable1,
                                   update_popsummary = FALSE,
                                   na.rm = na_var_rm,
                                   filters = NULL,
                                   poly_quantile = poly_quantile)
    
    plot_poly_by_pop = poly_by_Var1$pop_summary %>% 
      ggplot(aes(x = factor(pop, 
                            levels = c(unique(poly_by_Var1$pop_summary$pop)[unique(poly_by_Var1$pop_summary$pop) != 'Total'], "Total")),
                 y = prop_poly,
                 fill = factor(pop, 
                               levels = c(unique(poly_by_Var1$pop_summary$pop)[unique(poly_by_Var1$pop_summary$pop) != 'Total'], "Total"))))+
      geom_col(alpha = .6) +
      geom_errorbar(aes(ymin = prop_poly_lower, ymax = prop_poly_upper), width = .2)+
      theme_bw() +
      labs(title = "Frequency of polyclonal infections",
           y = "Frecquency") +
      scale_fill_manual(values = c(brewer.pal(length(poly_by_Var1$pop_summary$pop) - 1, 'Accent'), "gray30"))+
      theme(axis.text = element_text(size = 12),
            axis.title = element_blank(),
            legend.position = "none")
    
  }else{
    
    poly_total = get_polygenomic(ampseq_object = ampseq_object, 
                                   strata = NULL,
                                   update_popsummary = FALSE,
                                   na.rm = na_var_rm,
                                   filters = NULL,
                                   poly_quantile = poly_quantile)
    
  }
  
  
  
  if(!is.na(Variable1) & !is.na(Variable2)){
    
    ampseq_object@metadata[['Var1_Var2']] = paste(ampseq_object@metadata[[Variable1]], ampseq_object@metadata[[Variable2]], sep = '_')
    
    poly_by_Var1_Var2 = get_polygenomic(ampseq_object = ampseq_object,
                                        strata = "Var1_Var2",
                                        update_popsummary = F,
                                        na.rm = TRUE,
                                        filters = NULL)
    
    
    plot_poly_by_Var1_Var2 = poly_by_Var1_Var2$pop_summary %>%
      filter(pop != 'Total')%>%
      mutate(
        Variable1 = stringr::str_split(pop, '_', simplify = TRUE)[,1],
        Variable2 = stringr::str_split(pop, '_', simplify = TRUE)[,2],
        prop_poly_lower = case_when(
          prop_poly == 0 ~ 0,
          prop_poly != 0 ~ prop_poly_lower),
        prop_poly_upper = case_when(
          prop_poly == 0 ~ 0,
          prop_poly != 0 ~ prop_poly_upper)
      )%>%
      ggplot(aes(x = Variable2,
                 y = prop_poly,
                 ymin = prop_poly_lower,
                 ymax = prop_poly_upper,
                 fill = Variable1))+
      geom_col()+
      geom_errorbar(width = .2)+
      facet_wrap(~Variable1, ncol = 5)+
      theme_bw()+
      scale_fill_brewer(palette = 'Accent')+
      labs(title = 'Temporal change of the proportion of polyclonal infections',
           y = "Polyclonal infections",
           x = Variable2)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
            legend.position =  "none")
    
  }
  

  print('Generation of plots and tables for COI report done')
  
  imagename = file.path(wd,paste0(output, '_COI_Report.RData'))
  
  COI_Report_expected_outputs = c('poly_total',
                                  'poly_by_Var1',
                                  'plot_poly_by_pop',
                                  'poly_by_Var1_Var2',
                                  'plot_poly_by_Var1_Var2',
                                  'Variable1',
                                  'Variable2')
  
  COI_Report_outputs = COI_Report_expected_outputs[COI_Report_expected_outputs %in% ls()]
  
  save(file = imagename, list = COI_Report_outputs)
  
  system(paste0('cp ', file.path(fd, 'MHap_Analysis_COI_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_COI_Report.Rmd'))))
  
  # Assign variables based on command-line arguments
  render(file.path(wd, paste0(output, '_COI_Report.Rmd')), params = list(
    RData_image = imagename),
    output_dir = wd)
  
  print("Leaving render script")
  
}



# Genetic diversity----



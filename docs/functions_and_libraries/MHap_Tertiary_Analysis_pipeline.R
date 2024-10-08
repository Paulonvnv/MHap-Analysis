#!/bin/r env

library(argparse)
library(stringr)
library(rmarkdown)
# Parse arguments ----

parser = ArgumentParser()

# PATHs to references and functions

parser$add_argument("-wd", "--wd", default = 'null',
                    help="Path to input and output files or folders")

parser$add_argument("-fd", "--fd", default = 'null',
                    help="Path to function files")

parser$add_argument("-rd", "--rd", default = 'null',
                    help="Path to reference files")

# Upload of raw data

parser$add_argument("-ampseqe", "--ampseq_excelfile",  default = 'null',
                    help="name of the ampseq file in excel format")

parser$add_argument("-ampseqc", "--ampseq_csvfolder",  default = 'null',
                    help="Name of the folder where all the slots of the ampseq object are stored in csv files")

# Paramters to export data and results

parser$add_argument("-o", "--out",
                    help="string that define the prefix to be used for naming all output files and headers in the reports.")

parser$add_argument("-ampseq_export_format", "--ampseq_export_format",  default = 'null',
                    help='String that specify the format to export the filtered and masked ampseq object, options are "xlsx", "csv", and "json"')


# Removing samples or loci with low performance

parser$add_argument("-samprate", "--sample_ampl_rate", default = 0.75,
                    help="Min proportion of amplified loci by a sample that is required for the sample to be kept (default 0.75).")

parser$add_argument("-lamprate", "--locus_ampl_rate", default = 0.75,
                    help="Min proportion of amplified samples by a locus that is required for the locus to be kept (default 0.75).")

# Adding metadata

parser$add_argument("-metadata", "--metadata", default = 'null',
                    help="Name of the metadata file in csv format. If not required write null, and a metadata slot will be created in the AmpSeq object, the created metadata table will contain the samples ids (Sample_id), the name of the run plate, the order in the PCR plate and the type of sample (sample of interest or controls).")

parser$add_argument("-join_by", "--join_by", default = 'null',
                    help="string that indicates the variable where Sample ID's are stored. Just if  Sample ID's are not labeled as Sample_id")

parser$add_argument("-Var1", "--Variable1", default = 'null',
                    help="Name of variable in metadata that we want to use to group samples and split the outputs or summary metrics. For some visual representations geographic information is better to be written in Variable1, while temporal information should be in Variable2.")

parser$add_argument("-Var2", "--Variable2", default = 'null',
                    help="Name of variables in metadata that we want to use to group samples and split the outputs or summary metrics. For some visual representations geographic information is better to be written in Variable1, while temporal information should be in Variable2.")

parser$add_argument("-Longitude", "--Longitude", default = 'null',
                    help="Name of variables in the metadata that indicates the geographic coordinates in decimal degrees WGS84 system. These variables are only required visualization of information in maps. The coordinates can be at any geographic level (household, village, district, region, country, etc.), and depending on the geographic level defined by Variable1, all observations in each category will be summarized using their centroid.")

parser$add_argument("-Latitude", "--Latitude", default = 'null',
                    help="Name of variables in the metadata that indicates the geographic coordinates in decimal degrees WGS84 system. These variables are only required visualization of information in maps. The coordinates can be at any geographic level (household, village, district, region, country, etc.), and depending on the geographic level defined by Variable1, all observations in each category will be summarized using their centroid.")

# Filtering desired or undesired populations or samples

parser$add_argument("-var_filter", "--var_filter", default = 'null',
                    help='String used to filter desired categories for further analysis. This argument was created because sometimes we only want to generate the report of a subset of the total number of categories in our dataset (select specific locations or time points). The syntax for this input is as follows: “Variable1_Name;Keep;Category1,Category2,Category3/Varibale2_Name;Remove;Category1,Category2” Where the slatch bar the variables of interest, Varable1_Name and Variable2_Name is the name of the picked variables from the metadata table, “;” separates the name of the variables, the action to be performed (keep or remove, capital or lowercase letters are accepted) and the categories within the variables, and “,” separates the categories within variables. If the name of the categories or the variables contains blank spaces use dots “.” represent them (The requirement of the dot (“.”) is because an error argparse function in R and Python).')

# Reports

parser$add_argument("-PerformanceReport", "--PerformanceReport", default = FALSE,
                    help="Boolean. If this input is true, the filtering step is executed after the generation of the report. This reports includes: 1) A heatmap of the read depth of each marker and each sample, this heatmap is facet based on Variable1. 2)Dot and jitter plots of read depth across different categories (Variable1). 3) Histograms of the amplification rate of the samples and the amplification rate of the loci.")


parser$add_argument("-Drug_Surveillance_Report", "--Drug_Surveillance_Report", default = FALSE,
                    help="Boolean. If this is set to true, other inputs are required, and they are going to be explained in next slides. The output of this reports includes: 1) Report card maps for each desired drug, where the list of drugs is specified by the argument drug. 2) Line plots and stacked bar plots for the different phenotypes for each desired drug. 3) Line plots and stacked bar plots for the different haplotypes for each gene target. 4) Table with the haplotype and phenotype of each sample for each targeted gene.")


parser$add_argument("-Variants_of_Interest_Report", "--Variants_of_Interest_Report", default = FALSE,
                    help="Boolean. This report is similar to the Drug_Surveillance_Report, however it is intended for gene targets whose phenotypes are unknown as in the case of Plasmodium vivax or vaccine candidates.")

parser$add_argument("-ibd", "--ibd_thres", default = 'null',
                    help="Numerical value that set Minimum IBD to define highly related samples. If this parameter is different than null, then two reports are going to be written.")

parser$add_argument("-poly_formula", "--poly_formula", default = "NHetLoci>=1&Fws<1",
                    help='Formula used to define a polyclonal sample. (default: "NHetLoci>=1&Fws<1".')


# Parameters for DSR

parser$add_argument("-gff", "--ref_gff", default = 'null',
                    help="name of .gff file containing coordinates of genomic regions")

parser$add_argument("-fasta", "--ref_fasta", default = 'null',
                    help="name of .fasta file of the genome of reference")

parser$add_argument("-amp_fasta", "--amplicon_fasta", default = 'null',
                    help="name of .fasta file of the inserts of the amplicons on the reference strain")

parser$add_argument("-reference_alleles", "--reference_alleles", default = 'null',
                    help="name of .csv file containing sensitive alleles respect to a drug treatment")

parser$add_argument("-hap_color_palette", "--hap_color_palette", default = 'auto',
                    help="Character indicating how colors will be assigned to haplotypes in drug plots. If 'auto' (defualt), the red scales will be applied based on the presence of mutations associated with resistance, instead if 'random', a random palette will be generated")

parser$add_argument("-gene_names", "--gene_names", default = 'null',
                    help="gene name on the markers csv file")

parser$add_argument("-gene_ids", "--gene_ids", default = 'null',
                    help="gene ids on the .gff file")

parser$add_argument("-drugs", "--drugs", default = 'null',
                    help="Vector which allows to define the drugs that we want to screen and generate plots for the presence of resistant mutations.")

parser$add_argument("-include_all_drug_markers", "--include_all_drug_markers", default = FALSE,
                    help="boolean that allows to incorporate amplicons that were  discarded during the filtering process for having low amplification rate.")

parser$add_argument("-na_var_rm", "--na_var_rm", default = FALSE,
                    help="boolean that removes samples that have incomplete metadata (Variable1 or Variable2 is missing). This is another way to remove controls or samples in the cigar tables that belongs to other project or study sites that we don’t want to include in the final report, however, if there are spelling errors in the sample ids then we won’t be able to detect them.")

parser$add_argument("-na_hap_rm", "--na_hap_rm", default = FALSE,
                    help="boolean that remove samples with incomplete haplotypes from the calculus of the haplotypes frequencies and the plots from the final reports.")


# Parameters for IBD and Conectivity report


## Other parameters for IBD

parser$add_argument("-pairwise_relatedness_table", "--pairwise_relatedness_table", default = 'null',
                    help="string with the file name of the pairwise_relatedness_table")

parser$add_argument("-nchunks", "--nchunks", default = 500,
                    help="Number of chunks to subdivide the pairwise comparisons to reduce the consumption of RAM memory")

parser$add_argument("-parallel", "--parallel", default = FALSE,
                    help="Boolean to allow parallelization of the estimation of IBD")

parser$add_argument("-ibd_ncol", "--ibd_ncol", default = 4,
                    help="Number of column to use to draw different panels in plots")

parser$add_argument("-pop_levels", "--pop_levels", default = 'null',
                    help="Order in which categories of Variable1 and Variable2 are going to be display in plots")

# Parameters for COI report

parser$add_argument("-poly_quantile", "--poly_quantile", default = 'null',
                    help="Quantile to define polyclonal samples if this parameter is included in poly_formula")

# Parameters through shinny
parser$add_argument("-selected_checkboxes", "--selected_checkboxes", default = 'null',
                    help="Data points selected with the Shinny App")

parser$add_argument("-env", "--environment", default = 'terra',
                    help="Define if the code will be run in a 'local' machine, in 'cromwell' or in 'terra' (default).")

# Defining and checking variables ----

print("starting to parse variables")
args = parser$parse_args()

## Defining environment----

## PATHs to references and functions----

code_environment = gsub("/$", "", args$environment)

print(paste0('code_environment: ', code_environment))

# Working directory
if(tolower(code_environment) == 'local'){
  wd = gsub("/$", "", args$wd)
  print(paste0('wd: ', wd))
  setwd(wd)
}

# Tools or functions directory
fd = gsub("/$", "", args$fd)
print(paste0('fd: ', fd))

# Reference files directory

if(tolower(code_environment) == 'local'){
  rd = gsub("/$", "", args$rd)
  print(paste0('rd: ', rd))
}

## Upload of raw data----

# Starting ampseq file in excel format
ampseq_excelfile = args$ampseq_excelfile
ampseq_excelfile = if(ampseq_excelfile == 'null'){NULL}else{
  if(tolower(code_environment) == 'local'){
   file.path(wd, ampseq_excelfile) # for local computer
  }else{
    ampseq_excelfile # for Terra or Cromwell
  }
  }
print(paste0('ampseq_excelfile: ', ampseq_excelfile))

# Starting ampseq file in csv format
ampseq_csvfolder = args$ampseq_csvfolder
ampseq_csvfolder = if(ampseq_csvfolder == 'null'){NULL}else{
  if(tolower(code_environment) == 'local'){
   file.path(wd, ampseq_csvfolder)
  }else{
    ampseq_csvfolder
  }
}
print(paste0('ampseq_csvfolder: ', ampseq_csvfolder))


## Parameters to export data and results----

# output pattern
output = args$out
print(paste0('output: ', output))

# format to export the ampseq object
ampseq_export_format = args$ampseq_export_format
ampseq_export_format = if(ampseq_export_format == 'null'){
  NULL
  }else{
    as.character(ampseq_export_format)}

## Removing samples or loci with low performance----

# sample_ampl_rate
sample_ampl_rate = as.numeric(args$sample_ampl_rate)
print(paste0('sample_ampl_rate: ', sample_ampl_rate))

# locus_ampl_rate
locus_ampl_rate  = as.numeric(args$locus_ampl_rate)
print(paste0('locus_ampl_rate: ', locus_ampl_rate))


#Obtain data points for plotting----
selected_checkboxes = args$selected_checkboxes

if(selected_checkboxes == 'null'){NULL}else{
  content = readLines(selected_checkboxes)
  
  checkboxes_names = gsub('__', '', content[grep("__$", content)])
  
  list_of_selected_checkboxes = NULL
  
  for(checkbox in 1:length(checkboxes_names)){
    
    if(checkbox < length(checkboxes_names)){
      list_of_selected_checkboxes[[checkboxes_names[checkbox]]] = 
        content[(which(grepl(checkboxes_names[checkbox], content)) + 1):(which(grepl(checkboxes_names[checkbox + 1], content)) - 1)]
    }else{
      list_of_selected_checkboxes[[checkboxes_names[checkbox]]] = 
        content[(which(grepl(checkboxes_names[checkbox], content)) + 1):length(content)]
      
      }
    
  }
  
}


## Adding metadata----

# metadata
metadata_file = as.character(args$metadata)
metadata_file = if(metadata_file == 'null'){
  NULL
  }else{
    if(tolower(code_environment) == 'local'){
      file.path(wd, metadata_file)
      }else{
        metadata_file
      }
  }
print(paste0('metadata: ', metadata_file))

# join metadata by
join_by = as.character(args$join_by)
join_by = if(join_by == 'null'){NULL}else{join_by}
print(paste0('join_by: ', join_by))

# pop
Variable1 = args$Variable1
Variable1 = if(Variable1 == 'null'){NULL}else{Variable1}
print(paste0('Variable1: ', Variable1))

# temporal_population
Variable2 = args$Variable2
Variable2 = if(Variable2 == 'null'){NULL}else{Variable2}
print(paste0('Variable2: ', Variable2))

Longitude = args$Longitude
Longitude = if(Longitude == 'null'){NULL}else{Longitude}
print(paste0('Longitude: ', Longitude))

Latitude = args$Latitude
Latitude = if(Latitude == 'null'){NULL}else{Latitude}
print(paste0('Latitude: ', Latitude))

## Filtering desired or undesired populations or samples----

if(!is.null(selected_checkboxes)){
  
  subnational_level = list_of_selected_checkboxes[[3]]
  period_of_collection = list_of_selected_checkboxes[[4]]
  
  var_filter = c(paste(c(Variable1, 
                         paste(gsub('\\.', ' ', subnational_level), 
                               collapse = ",")), 
                       collapse = ";"),
                 paste(c(Variable2, 
                         paste(period_of_collection, 
                               collapse = ",")), 
                       collapse =";"))
  
}else{
  
  var_filter  = args$var_filter
  if(var_filter == 'null'){
    var_filter = NULL
  }else{
    var_filter = gsub('\\.', ' ', var_filter)
    var_filter = strsplit(var_filter, '/')[[1]]
  }
  
}

print('var_filter: ')
print(var_filter)

## Reports----

# Performance Report
PerformanceReport = as.logical(args$PerformanceReport)
print(paste0('PerformanceReport: ', PerformanceReport))

# Drug surveillance report
Drug_Surveillance_Report = as.logical(args$Drug_Surveillance_Report)
print(paste0('Drug_Surveillance_Report: ', Drug_Surveillance_Report))

# Variants of interest report (in development)
Variants_of_Interest_Report = as.logical(args$Variants_of_Interest_Report)
print(paste0('Variants_of_Interest_Report: ', Variants_of_Interest_Report))

# Reports of IBD and Connectivity and IBD and transmission
ibd_thres = args$ibd_thres
ibd_thres = if(ibd_thres == 'null'){NULL}else{as.numeric(ibd_thres)}
print(paste0('ibd_thres: ', ibd_thres))

# Report for complexity of infection

# poly_formula filter
poly_formula = as.character(args$poly_formula)
poly_formula = if(poly_formula == 'null'){NULL}else{as.character(poly_formula)}

if(!is.null(poly_formula)){

  poly_formula = gsub('"',"",poly_formula)
  poly_formula = gsub('&'," & ",poly_formula, ignore.case = TRUE)
  poly_formula = gsub('\\|'," \\| ",poly_formula, ignore.case = TRUE)
  
  if(grepl("\\w>\\d",poly_formula)){
    patterns = str_extract_all(poly_formula, "\\w>\\d")[[1]]
    
    for(pattern in patterns){
      
      replacement = gsub('>',' > ',pattern)
      poly_formula = gsub(pattern,
                          replacement,
                          poly_formula, ignore.case = TRUE)
    }
    
  }
  
  if(grepl("\\w<\\d",poly_formula)){
    patterns = str_extract_all(poly_formula, "\\w<\\d")[[1]]
    
    for(pattern in patterns){
      
      replacement = gsub('<',' < ',pattern)
      poly_formula = gsub(pattern,
                          replacement,
                          poly_formula, ignore.case = TRUE)
    }
    
  }
  
  poly_formula = gsub('>='," >= ", poly_formula, ignore.case = TRUE)
  poly_formula = gsub('<='," <= ", poly_formula, ignore.case = TRUE)
  poly_formula = gsub('=='," == ", poly_formula, ignore.case = TRUE)
  poly_formula = gsub('!='," != ", poly_formula, ignore.case = TRUE)
  
  poly_formula = gsub('\\+'," \\+ ", poly_formula, ignore.case = TRUE)
  poly_formula = gsub('-'," - ", poly_formula, ignore.case = TRUE)
  poly_formula = gsub('\\*'," \\* ", poly_formula, ignore.case = TRUE)
  poly_formula = gsub('/'," / ", poly_formula, ignore.case = TRUE)
}

print(paste0('poly_formula: ', poly_formula))

### Parameters for DSR----

# Reference gff and fasta files, require to translate DNA cigar formats to aminoacid sequences

ref_gff = args$ref_gff
ref_gff = if(ref_gff == 'null'){NULL}else{
  if(tolower(code_environment) == 'local'){
    file.path(rd, ref_gff)  
  }else{
    ref_gff 
  }
  }
print(paste0('ref_gff: ', ref_gff))

ref_fasta = args$ref_fasta
ref_fasta = if(ref_fasta == 'null'){NULL}else{
  if(tolower(code_environment) == 'local'){
   file.path(rd, ref_fasta)}else{
    ref_fasta
   }
  }
print(paste0('ref_fasta: ', ref_fasta))

amplicon_fasta = args$amplicon_fasta
amplicon_fasta = if(amplicon_fasta == 'null'){NULL}else{
  if(tolower(code_environment) == 'local'){
   file.path(rd, amplicon_fasta)
  }else{amplicon_fasta}
  }
print(paste0('amplicon_fasta: ', amplicon_fasta))

reference_alleles = args$reference_alleles
reference_alleles = if(reference_alleles == 'null'){NULL}else{
  if(tolower(code_environment) == 'local'){
   file.path(rd, reference_alleles)
    }else{
      reference_alleles}
  }
print(paste0('reference_alleles: ', reference_alleles))

hap_color_palette = args$hap_color_palette
print(paste0('hap_color_palette: ', hap_color_palette))


# gene_names

if(!is.null(selected_checkboxes)){
  gene_names = list_of_selected_checkboxes[[1]]
}else{
  gene_names = args$gene_names
  if(gene_names == 'null'){
    gene_names = NULL
  }else{
    gene_names = strsplit(gene_names, ',')[[1]]
  }
}

print('gene_names:')
print(gene_names)

# gene_ids

if(!is.null(selected_checkboxes)){
  gene_ids = list_of_selected_checkboxes[[2]]
}else{
  gene_ids = args$gene_ids
  if(gene_ids == 'null'){
    gene_ids = NULL
  }else{
    gene_ids = strsplit(gene_ids, ',')[[1]]
  }
}

print('gene_ids:')
print(gene_ids)

drugs = args$drugs
if(drugs == 'null'){
  drugs = NULL
}else{
  drugs = strsplit(drugs, ',')[[1]]
}
print('drugs: ')
print(drugs)

include_all_drug_markers = as.logical(args$include_all_drug_markers)
print(paste0('include_all_drug_markers: ', include_all_drug_markers))

na_var_rm = as.logical(args$na_var_rm)
print(paste0('na_var_rm: ', na_var_rm))
      
na_hap_rm = as.logical(args$na_hap_rm)
print(paste0('na_hap_rm: ', na_hap_rm))


### Parameters for IBD and Conectivity report----

# pairwise_relatedness_table
pairwise_relatedness_table = as.character(args$pairwise_relatedness_table)
pairwise_relatedness_table = if(pairwise_relatedness_table == 'null'){
  NULL
  }else{
    if(tolower(code_environment) == 'local'){
     file.path(wd, pairwise_relatedness_table)
      }else{
        pairwise_relatedness_table
      }
  }

print(paste0('pairwise_relatedness_table: ', pairwise_relatedness_table))

# parallel
parallel = as.logical(args$parallel)
print(paste0('parallel: ', parallel))

# nchunks
nchunks = as.integer(args$nchunks)
print(paste0('nchunks: ', nchunks))

# ibd_ncol
ibd_ncol = as.integer(args$ibd_ncol)
print(paste0('ibd_ncol: ', ibd_ncol))

# pop_levels
pop_levels = args$pop_levels
pop_levels = if(pop_levels == 'null'){
  NULL
}else{
  as.character(pop_levels)}
print(paste0('pop_levels: ', pop_levels))

### Parameters for COI report----
# poly_quantile
poly_quantile = args$poly_quantile
poly_quantile = if(poly_quantile == 'null'){
  NULL
}else{
    as.numeric(poly_quantile)}
print(paste0('poly_quantile: ', poly_quantile))


print("All variables checked")

# Check packages and functions----

print("Loading libraies and functions")
source(file.path(fd,'amplseq_required_libraries.R'))
source(file.path(fd,'amplseq_functions.R'))

# Defining a random color palette for all plots with multiple categorical variables----

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Upload ampseq object----

print("starting to upload genotypic data")

if(!is.null(ampseq_excelfile)){
  print("Uploading genotypic data in ampseq format from excel file")
  ampseq_object = read_ampseq(file = ampseq_excelfile, format = 'excel')
  
  }else if(!is.null(ampseq_csvfolder)){
    print("Uploading genotypic data in ampseq format from excel file")
    ampseq_object = read_ampseq(file = ampseq_csvfolder, format = 'csv')
    
  }

# Filtering and adding metadata to the  ampseq_object----

if(PerformanceReport){
  
  # add metadata 
  
  if(!is.null(metadata_file)){
    # Upload metadata from an external source
    external_metadata = read.csv(metadata_file)
    
    # Check duplicate samples
    
    if(sum(duplicated(external_metadata$Sample_id)) > 0){
      
      print(paste0('There are ', sum(duplicated(external_metadata$Sample_id)), ' duplicated Ids that are going to be removed'))
      external_metadata %<>% filter(!duplicated(external_metadata$Sample_id))
      
    }
    
    
      shared_variables = names(ampseq_object@metadata)[names(ampseq_object@metadata) %in% names(external_metadata)]
      shared_variables = shared_variables[!(shared_variables %in% c('Sample_id', join_by))]

      if(length(shared_variables) > 0){
        print(paste0('The following variables will be overwritten: ', paste(shared_variables, collapse = ', ')))

        for(variable in shared_variables){

          ampseq_object@metadata[[variable]] = NULL

        }

      }

    print("Adding metadata to ampseq object")
      # Merge the external metadata with our ampseq_object
    ampseq_object@metadata = left_join(ampseq_object@metadata,
                                                external_metadata,
                                                  by = join_by)
    }
    
    
  
  ## Sample performance by different coverage----
  ### Overall sample performance by different coverage----
  
  ReadDepth_coverage = get_ReadDepth_coverage(ampseq_object, variable = NULL)
  
  min_abd = min(ReadDepth_coverage$plot_read_depth_heatmap$data$Read_depth, na.rm = T)
  
  sample_performance = ReadDepth_coverage$plot_read_depth_heatmap$data %>%
    mutate(Read_depth = case_when(
      is.na(Read_depth) ~ 0,
      !is.na(Read_depth) ~ Read_depth
    )) %>%
    summarise(amplified_amplicons1 = sum(Read_depth >= 1)/nrow(ampseq_object@markers),
              amplified_amplicons5 = sum(Read_depth >= 5)/nrow(ampseq_object@markers),
              amplified_amplicons10 = sum(Read_depth >= 10)/nrow(ampseq_object@markers),
              amplified_amplicons20 = sum(Read_depth >= 20)/nrow(ampseq_object@markers),
              amplified_amplicons50 = sum(Read_depth >= 50)/nrow(ampseq_object@markers),
              amplified_amplicons100 = sum(Read_depth >= 100)/nrow(ampseq_object@markers),
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
                                 summarise(AmpRate = sum(Read_depth >= min_abd)/nrow(ampseq_object@markers),
                                           .by = Sample_id) %>%
                                 mutate(Threshold = min_abd) %>%
                                 select(Sample_id, Threshold, AmpRate))
    
  }
  
  sample_performance_summary = sample_performance %>%
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
    )
  
  amplified_samples = sample_performance_summary[min_abd == sample_performance_summary$Threshold, paste0("AmpRate", sample_ampl_rate*100)]
  
  
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
  
  plot_precentage_of_samples_over_min_abd_reduced = sample_performance_summary[min_abd == sample_performance_summary$Threshold, ] %>%
    pivot_longer(cols = paste0('AmpRate', seq(5, 100, 5)),
                 values_to = 'Percentage',
                 names_to = 'AmpRate') %>%
    mutate(AmpRate = as.numeric(gsub('AmpRate','', AmpRate)))%>%
    ggplot(aes(x = AmpRate, y = Percentage, color = as.factor(Threshold), group = as.factor(Threshold))) +
    geom_line() +
    geom_vline(xintercept = 100*sample_ampl_rate, linetype = 2) +
    theme_bw() +
    labs(x = '% of amplified loci (amplification rate)', y = '% of Samples', color = 'Min Coverage')
  
  #WATERMARKED
  plot_precentage_of_samples_over_min_abd_AMPLIFICATION = sample_performance %>%
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
    mutate(AmpRate = as.numeric(gsub('AmpRate','', AmpRate))) %>%
    ggplot(aes(x = Percentage, y = AmpRate, color = as.factor(Threshold), group = as.factor(Threshold))) +
    geom_line() +
    geom_hline(yintercept = 100*sample_ampl_rate, linetype = 2) +
    theme_bw() +
    labs(x = '% of Samples', y = '% of amplified loci (amplification rate)', color = 'Min Coverage')
  
  
  ### Overall sample performance by different coverage per Run----
  
  ReadDepth_coverage = get_ReadDepth_coverage(ampseq_object, variable = 'Run')
  
  sample_performance = ReadDepth_coverage$plot_read_depth_heatmap$data %>%
    mutate(Read_depth = case_when(
      is.na(Read_depth) ~ 0,
      !is.na(Read_depth) ~ Read_depth
    )) %>%
    summarise(amplified_amplicons1 = sum(Read_depth >= 1)/nrow(ampseq_object@markers),
              amplified_amplicons5 = sum(Read_depth >= 5)/nrow(ampseq_object@markers),
              amplified_amplicons10 = sum(Read_depth >= 10)/nrow(ampseq_object@markers),
              amplified_amplicons20 = sum(Read_depth >= 20)/nrow(ampseq_object@markers),
              amplified_amplicons50 = sum(Read_depth >= 50)/nrow(ampseq_object@markers),
              amplified_amplicons100 = sum(Read_depth >= 100)/nrow(ampseq_object@markers),
              Run = unique(var),
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
                                 summarise(AmpRate = sum(Read_depth >= min_abd)/nrow(ampseq_object@markers),
                                           Run = unique(var),
                                           .by = Sample_id) %>%
                                 mutate(Threshold = min_abd) %>%
                                 select(Sample_id, Run, Threshold, AmpRate))
    
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
              .by = c(Threshold, Run)
    ) %>%
    pivot_longer(cols = paste0('AmpRate', seq(5, 100, 5)),
                 values_to = 'Percentage',
                 names_to = 'AmpRate') %>%
    mutate(AmpRate = as.numeric(gsub('AmpRate','', AmpRate)))%>%
    ggplot(aes(x = AmpRate, y = Percentage, color = as.factor(Threshold), group = as.factor(Threshold))) +
    geom_line() +
    geom_vline(xintercept = 100*sample_ampl_rate, linetype = 2) +
    facet_wrap(Run~., ncol = 3)+
    theme_bw() +
    labs(x = '% of amplified loci (amplification rate)', y = '% of Samples', color = 'Min Coverage')
  
  
  ### Overall sample performance by different coverage per Variable1----
  
  if(!is.null(Variable1)){
    
    ReadDepth_coverage = get_ReadDepth_coverage(ampseq_object, variable = Variable1)
    
    sample_performance = ReadDepth_coverage$plot_read_depth_heatmap$data %>%
      filter(!is.na(var))%>%
      mutate(Read_depth = case_when(
        is.na(Read_depth) ~ 0,
        !is.na(Read_depth) ~ Read_depth
      )) %>%
      summarise(amplified_amplicons1 = sum(Read_depth >= 1)/nrow(ampseq_object@markers),
                amplified_amplicons5 = sum(Read_depth >= 5)/nrow(ampseq_object@markers),
                amplified_amplicons10 = sum(Read_depth >= 10)/nrow(ampseq_object@markers),
                amplified_amplicons20 = sum(Read_depth >= 20)/nrow(ampseq_object@markers),
                amplified_amplicons50 = sum(Read_depth >= 50)/nrow(ampseq_object@markers),
                amplified_amplicons100 = sum(Read_depth >= 100)/nrow(ampseq_object@markers),
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
                                   summarise(AmpRate = sum(Read_depth >= min_abd)/nrow(ampseq_object@markers),
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
    
    data = sample_performance %>% select(1:2) %>% distinct() %>% count(!!sym(names(sample_performance)[2]))
    
    samples_by_location = ggplot(data, aes(x = !!sym(names(data)[1]), y = n, fill = !!sym(names(data)[1]))) +
      geom_bar(stat = "identity", color = "black", width = 0.7) + # Bar colors with spacing
      labs(x = "Location", y = "Number of Samples") + # Axis labels
      theme_minimal() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels if necessary
            panel.grid.major.x = element_blank(), # Remove vertical grid lines
            plot.title = element_text(hjust = 0.5)) + 
      scale_fill_discrete(name = "Location")
    
    
    
  }
  
}
  

## Coverage by sample and amplicon----


if(PerformanceReport == TRUE){
  
  print('Measuring coverage per sample per amplicon')
  
  # Read Depth coverage by the whole dataset or by Variable1
  
  if(is.null(Variable1)){
    ReadDepth_coverage = get_ReadDepth_coverage(ampseq_object, variable = NULL)
  }else{
    ReadDepth_coverage = get_ReadDepth_coverage(ampseq_object, variable = Variable1)
  }
  
  # Read Depth coverage by the whole dataset or by Run
  
  ReadDepth_coverage_by_run = get_ReadDepth_coverage(ampseq_object, variable = "Run")
  
  # Read Depth coverage by the whole dataset or by Run for controls
  
  # ampseq_object_controls = filter_samples(ampseq_object, v = ampseq_object@metadata$typeofSamp == 'Controls')
  # 
  # ReadDepth_coverage_by_run_controls = get_ReadDepth_coverage(ampseq_object = ampseq_object_controls, variable = "Run")
  # 
}

# Remove controls

print('Removing Controls for further steps')

ampseq_object = filter_samples(ampseq_object, v = ampseq_object@metadata$typeofSamp == 'Samples')

# if(PerformanceReport){
#   ampseq_object@controls[['gt']] = ampseq_object_controls@gt
#   ampseq_object@controls[['metadata']] = ampseq_object_controls@metadata
# }

print("Removing undesired categories based on var_filter")
if(!is.null(var_filter)){# ADAPT THIS SECTION IN THE SHINY APP
  filters = strsplit(var_filter,';')
  for(temp_filter in 1:length(filters)){
    
    #if(toupper(filters[[temp_filter]][2]) == 'KEEP'){
    
    ampseq_object = filter_samples(ampseq_object,
                                   ampseq_object@metadata[[filters[[temp_filter]][1]]] %in% strsplit(filters[[temp_filter]][2],',')[[1]])
    
    #}else if(toupper(filters[[temp_filter]][2]) == 'REMOVE'){
    
    #  ampseq_object = filter_samples(ampseq_object,
    #                                 !(ampseq_object@metadata[[filters[[temp_filter]][1]]] %in% strsplit(filters[[temp_filter]][3],',')[[1]]))
    
    #}
  }
}

# if locus_ampl_rate was provided

print('Measuring amplification rate by locus')

if(!is.null(locus_ampl_rate)){
  
  if(!is.null(Variable1)){
    ampseq_object = locus_amplification_rate(ampseq_object, threshold = locus_ampl_rate, strata = Variable1)
  }else{
    ampseq_object = locus_amplification_rate(ampseq_object, threshold = locus_ampl_rate)
  }
  
}else{
  
  if(!is.null(Variable1)){
    ampseq_object = locus_amplification_rate(ampseq_object, threshold = 0.65, strata = Variable1)
  }else{
    ampseq_object = locus_amplification_rate(ampseq_object, threshold = 0.65) 
  }
}

# if sample_ampl_rate was provided

print("measuring amplification rate by sample")

if(!is.null(sample_ampl_rate)){
  if(!is.null(Variable1)){
    ampseq_object = sample_amplification_rate(ampseq_object, threshold = sample_ampl_rate, strata = Variable1)
  }else{
    ampseq_object = sample_amplification_rate(ampseq_object, threshold = sample_ampl_rate)  
  }
  
}else{
  
  if(!is.null(Variable1)){
    ampseq_object = sample_amplification_rate(ampseq_object, threshold = 0.80, strata = Variable1)
  }else{
    ampseq_object = sample_amplification_rate(ampseq_object, threshold = 0.80)  
  }

}

## Printing performance report----

if(PerformanceReport == TRUE){
  
  all_loci_amplification_rate = ampseq_object@plots$all_loci_amplification_rate
  
  samples_amplification_rate = ampseq_object@plots$samples_amplification_rate
  
  cigar_table_unmasked_unfiltered = ampseq_object@gt 
  cigar_table_masked_filtered = ampseq_object@gt
  # cigar_table_controls_masked_filtered = ampseq_object_controls@gt
  
  metadata_kept_samples = ampseq_object@metadata
  metadata_removed_samples  = ampseq_object@discarded_samples$metadata
  
  
  print('Generation of plots and tables for Performance report done')
  
  
  Performance_Report_expected_outputs = c(
    'cigar_table_unmasked_unfiltered',
    'cigar_table_masked_filtered',
    # 'cigar_table_controls_masked_filtered',
    
    'metadata_kept_samples',
    'metadata_removed_samples',
    
    'plot_precentage_of_samples_over_min_abd',
    'plot_precentage_of_samples_over_min_abd_byRun',
    'plot_precentage_of_samples_over_min_abd_byVariable1',
    
    'ReadDepth_coverage',
    'ReadDepth_coverage_by_run',
    # 'ReadDepth_coverage_by_run_controls',
    'all_loci_amplification_rate',
    'samples_amplification_rate'
    )
  
  # Performance_Report_outputs = Performance_Report_expected_outputs[Performance_Report_expected_outputs %in% ls()]
  # 
  # if(tolower(code_environment) == 'local'){
  #   imagename = file.path(wd,paste0(output, '_Performance_Report.RData'))
  # 
  #   save(file = imagename, list = Performance_Report_outputs)
  # 
  #   system(paste0('cp ', file.path(fd, 'MHap_Analysis_Performance_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_Performance_Report.Rmd'))))
  # 
  #   # Assign variables based on command-line arguments
  #   render(file.path(wd, paste0(output, '_Performance_Report.Rmd')), params = list(
  #     RData_image = imagename),
  #     output_dir = wd)
  # }else{
  #   
  #   imagename = paste0(output, '_Performance_Report.RData')
  #   
  #   save(file = imagename, list = Performance_Report_outputs)
  #   
  #   system(paste0('cp ', file.path(fd, 'MHap_Analysis_Performance_Report_Template.Rmd'), ' ', paste0(output, '_Performance_Report.Rmd')))
  #   
  #   # Assign variables based on command-line arguments
  #   render(paste0(output, '_Performance_Report.Rmd'), params = list(
  #     RData_image = imagename),
  #     output_dir = 'Results')
  #   
  # }
  
  print("Leaving render script")
  
}

# Drug Resistance surveillance ----

if(Drug_Surveillance_Report){
  
  print('Starting drug surveillance report')
  
  variables = 'Sample_id'
  
  if(!is.null(Variable1)){
    variables = c(variables, Variable1)
  }
  
  if(!is.null(Variable2)){
    variables = c(variables, Variable2)
  }
  
  
  if(include_all_drug_markers){
    
    print('including drug markers with low amplification rate')
    ampseq_drug = ampseq_object
    
    ampseq_drug@gt = cbind(ampseq_drug@gt,
                           ampseq_drug@discarded_loci$gt[rownames(ampseq_drug@discarded_loci$gt) %in%
                                                           rownames(ampseq_drug@gt),
                                                         grepl(paste0(gene_names, collapse = '|'),
                                                               colnames(ampseq_drug@discarded_loci$gt))]
    )
    
    ampseq_drug@markers = rbind(ampseq_drug@markers,
                                ampseq_drug@discarded_loci$markers[grepl(paste0(gene_names, collapse = '|'),
                                                                         ampseq_drug@discarded_loci$markers$amplicon),])
    
  }else{
    
    ampseq_drug = ampseq_object
    
  }
  
  # Outputs for full report
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
                                                             filters = NULL,
                                                             hap_color_palette = hap_color_palette)
  
  minimal_aa_table = drug_resistant_haplotypes_plot$minimal_aa_table
  
  # Outputs for minimal report
  # tmp = drug_resistant_haplotypes_plot$aa_mutations
  # tmp = as.data.frame(tmp)
  # tmp$Sample_Id = row.names(tmp)
  # 
  # 
  # 
  # View(drug_resistant_haplotypes_plot$genotype_phenotype_table)
  # 
  # View(drug_resistant_haplotypes_plot$haplotype_freq_barplot$data)
  # View(drug_resistant_haplotypes_plot$drug_phenotype_barplot$data)
  # 
  # # Convert to long format
  # long_df <- tmp %>%
  #   pivot_longer(
  #     cols = -Sample_Id,            # Pivot all columns
  #     names_to = "Sample",            # New column for row names
  #     values_to = "Marker"            # New column for values
  #   )  %>%
  #   separate_rows(Marker, sep = " ")
  # 
  # # View the resulting long format table
  # #print(long_df)
  # 
  # drugR_reference_alleles_tmp = read.csv(reference_alleles)
  # drugR_reference_alleles_tmp$Description <- iconv(drugR_reference_alleles_tmp$Description, from = "latin1", to = "ASCII//TRANSLIT")
  # drugR_reference_alleles_tmp$Description <- gsub("\\^E$", "", drugR_reference_alleles_tmp$Description)
  # print(drugR_reference_alleles_tmp)
  # 
  # partial = long_df %>% 
  #   left_join(drugR_reference_alleles_tmp, by = c("Marker" = "Mutation")) %>%
  #   select(Sample_Id, Sample, Marker, Annotation) %>%
  #   filter(!is.na(Annotation)) %>% 
  #   select(-Sample_Id) %>%     # Remove the Sample_Id column
  #   group_by(Sample, Marker, Annotation) %>% # Group by the relevant columns
  #   summarise(Count = n(), .groups = 'drop') %>%  # Count repeated rows and add as a new column
  #   relocate(Annotation, .before = Sample) #%>%    # Move the Annotation column to the left
  # #rename(Drug = Annotation) 
  # 
  # #TABLE NAMES ADD
  # #Performance plot. Number that went in vs Number that went out. Changing wording.
  # 
  # #Drug Phenotype
  # #Locus
  # #Mutation/Marker
  # #Count
  # #Analyzed Samples/Total Samples/ Samples that yield intepretable data.
  # #Circle is the number of interpretable samples.
  # denominator = long_df %>% 
  #   left_join(drugR_reference_alleles_tmp, by = c("Marker" = "Mutation")) %>%
  #   select(Sample_Id, Sample, Marker, Annotation) %>%
  #   select(-Sample_Id) %>%     # Remove the Sample_Id column
  #   group_by(Sample) %>% # Group by the relevant columns
  #   summarise(Sample_Space = n(), .groups = 'drop') 
  # 
  # combined_df <- partial %>%
  #   left_join(denominator, by = "Sample")
  # 
  # # View the result
  # print(combined_df)
  # 
  # partial = as.data.frame(combined_df)
  # colnames(partial) = c("Resistance Association of gene", "Gene", "Marker", "Count", "Samples Genotyped")
  # 
  
  print('Generation of plots and tables for DRS report done')
  
  if(tolower(code_environment) == 'local'){
    imagename = file.path(wd,paste0(output, '_DRS_Report.RData'))
    save(file = imagename, list = c('drug_resistant_haplotypes_plot',
      'Longitude',
      'Latitude',
      'samples_by_location',
      'plot_precentage_of_samples_over_min_abd_reduced',
      'sample_ampl_rate',
      'min_abd',
      'amplified_samples',
      'minimal_aa_table'))
  
    # Full report
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_DRS_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_DRS_Report.Rmd'))))
    # Assign variables based on command-line arguments
    render(file.path(wd, paste0(output, '_DRS_Report.Rmd')), params = list(
      RData_image = imagename),
    output_dir = wd)
    # Minimal report
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_DRS_Minimal_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_DRS_Minimal_Report.Rmd'))))
    # Assign variables based on command-line arguments
    render(file.path(wd, paste0(output, '_DRS_Minimal_Report.Rmd')), params = list(
      RData_image = imagename),
      output_dir = wd)
  }else{
  
    imagename = paste0(output, '_DRS_Report.RData')
    save(file = imagename, list = c('drug_resistant_haplotypes_plot', 
                                    'Longitude', 
                                    'Latitude', 
                                    'samples_by_location', 
                                    'plot_precentage_of_samples_over_min_abd_reduced',
                                    'sample_ampl_rate',
                                    'min_abd',
                                    'amplified_samples',
                                    'minimal_aa_table'))
    
    # Full report
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_DRS_Report_Template.Rmd'), ' ', paste0(output, '_DRS_Report.Rmd')))
    # Assign variables based on command-line arguments
    render(paste0(output, '_DRS_Report.Rmd'), params = list(
      RData_image = imagename),
      output_dir = 'Results')
    # Minimal report
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_DRS_Minimal_Report_Template.Rmd'), ' ', paste0(output, '_DRS_Minimal_Report.Rmd')))
    # Assign variables based on command-line arguments
    render(paste0(output, '_DRS_Minimal_Report.Rmd'), params = list(
      RData_image = imagename),
      output_dir = 'Results')
  }
  print("Leaving render script")
  
}

# Variants of interest ----

if(Variants_of_Interest_Report){ # IN DEVELOPMENT
  
  # variables = 'Sample_id'
  # 
  # if(!is.null(Variable1)){
  #   variables = c(variables, Variable1)
  # }
  # 
  # if(!is.null(Variable2)){
  #   variables = c(variables, Variable2)
  # }
  # 
  # haplotypes_respect_to_reference_plot = haplotypes_respect_to_reference(ampseq_object,
  #                                                            gene_names = gene_names,
  #                                                            gene_ids = gene_ids,
  #                                                            gff_file = ref_gff,
  #                                                            fasta_file = ref_fasta,
  #                                                            variables = variables,
  #                                                            plot_haplo_freq = TRUE,
  #                                                            na.var.rm = na_var_rm,
  #                                                            filters = NULL)
  # 
  # imagename = file.path(wd,paste0(output, '_VoI_Report.RData'))
  # 
  # save(file = imagename, list = c('haplotypes_respect_to_reference_plot'))
  
}


# Genetic Variants of interest---

# Genetic Relatedness ----

if(!is.null(ibd_thres)){
  
  # call hmmIBD and PCA functions from Rcpp
  sourceCpp(file.path(fd,'hmmloglikelihood.cpp'))
  
  if(is.null(pairwise_relatedness_table)){
    pairwise_relatedness = NULL
    
    for(w in nChunks){
      start = Sys.time()
      pairwise_relatedness = rbind(pairwise_relatedness,
                                   pairwise_hmmIBD(ampseq_object, parallel = parallel, w = w, n = nChunks))
      time_diff = Sys.time() - start
      
      print(paste0('step ', w, ' done in ', time_diff, ' secs'))
      
    }
    
    if(tolower(code_environment) == 'local'){
      write.csv(pairwise_relatedness,
                file.path(wd, paste0(output, '_pairwise_ibd', '.csv')),
                quote = FALSE,
                row.names = FALSE)
    }else{

      write.csv(pairwise_relatedness,
                paste0(output, '_pairwise_ibd', '.csv'),
                quote = FALSE,
                row.names = FALSE)
    }
    
  }else{
    
    pairwise_relatedness = read.csv(pairwise_relatedness_table)
    
  }
  
  
  ## Genetic connectivity----
  
  plot_relatedness_distribution_between = plot_relatedness_distribution(
    pairwise_relatedness = pairwise_relatedness,
    metadata = ampseq_object@metadata,
    Population = Variable1,
    fill_color = rep('gray50', length(unique(ampseq_object@metadata[[Variable1]]))*(length(unique(ampseq_object@metadata[[Variable1]]))-1)/2),
    type_pop_comparison = 'between',
    ncol = 3,
    pop_levels = NULL
  )
  
  
  plot_frac_highly_related_between = plot_frac_highly_related(
    pairwise_relatedness = pairwise_relatedness,
    metadata = ampseq_object@metadata,
    Population = Variable1,
    fill_color = rep('gray50', length(unique(ampseq_object@metadata[[Variable1]]))*(length(unique(ampseq_object@metadata[[Variable1]]))-1)/2),
    threshold = ibd_thres,
    type_pop_comparison = 'between',
    pop_levels = NULL)
  
  
  if(!is.null(Variable2)){
    plot_frac_highly_related_overtime_between = plot_frac_highly_related_over_time(
      pairwise_relatedness = pairwise_relatedness,
      metadata = ampseq_object@metadata,
      Population = c(Variable1, Variable2),
      fill_color = rep('gray50', length(unique(ampseq_object@metadata[[Variable1]]))*(length(unique(ampseq_object@metadata[[Variable1]]))-1)/2),
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
  
  set.seed(1)
  
  IBD_PCA = evectors_IBD$eigenvector %>% ggplot(aes(x = PC1, y = PC2, color = Variable1))+
    geom_point(alpha = .7, size = 2) +
    stat_ellipse(level = .6)+
    scale_color_manual(values = sample(col_vector, nlevels(as.factor(ampseq_object@metadata[[Variable1]]))))+
    theme_bw()+
    labs(x = paste0('1st PCo (', round(evectors_IBD$contrib[1],1), '%)'),
         y = paste0('2nd PCo (', round(evectors_IBD$contrib[2],1), '%)'),
         color = 'Countries')
  
  
  
  print('Generation of plots and tables for IBD and Connectivity report done')
  
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
  
  if(tolower(code_environment) == 'local'){
    imagename = file.path(wd,paste0(output, '_IBD_Connectivity_Report.RData'))
  
    save(file = imagename, list = IBD_Connectivity_Report_outputs)
  
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_IBD_Connectivity_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_IBD_Connectivity_Report.Rmd'))))
  
    # Assign variables based on command-line arguments
    render(file.path(wd, paste0(output, '_IBD_Connectivity_Report.Rmd')), params = list(
      RData_image = imagename),
      output_dir = wd)
  }else{
  
    imagename = paste0(output, '_IBD_Connectivity_Report.RData')
    
    save(file = imagename, list = IBD_Connectivity_Report_outputs)
    
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_IBD_Connectivity_Report_Template.Rmd'), ' ', paste0(output, '_IBD_Connectivity_Report.Rmd')))
    
    # Assign variables based on command-line arguments
    render(paste0(output, '_IBD_Connectivity_Report.Rmd'), params = list(
      RData_image = imagename),
      output_dir = 'Results')
  }
  print("Leaving render script")
  
  
  ## Relatedness and transmission----
  
  plot_relatedness_distribution_within = plot_relatedness_distribution(
    pairwise_relatedness = pairwise_relatedness,
    metadata = ampseq_object@metadata,
    Population = Variable1,
    fill_color = rep('gray50', length(unique(ampseq_object@metadata[[Variable1]]))),
    type_pop_comparison = 'within',
    ncol = 3,
    pop_levels = NULL
  )
  
  
  plot_frac_highly_related_within = plot_frac_highly_related(
    pairwise_relatedness = pairwise_relatedness,
    metadata = ampseq_object@metadata,
    Population = Variable1,
    fill_color = rep('gray50', length(unique(ampseq_object@metadata[[Variable1]]))),
    threshold = ibd_thres,
    type_pop_comparison = 'within',
    pop_levels = NULL)
  
  
  if(!is.null(Variable2)){
    plot_frac_highly_related_overtime_within = plot_frac_highly_related_over_time(
      pairwise_relatedness = pairwise_relatedness,
      metadata = ampseq_object@metadata,
      Population = c(Variable1, Variable2),
      fill_color = rep('gray50', length(unique(ampseq_object@metadata[[Variable1]]))),
      threshold = ibd_thres,
      type_pop_comparison = 'within',
      ncol = 3,
      pop_levels = NULL)
  }
  
  print('Generation of plots and tables for IBD and Transmission report done')
  
  IBD_Transmission_Report_expected_outputs = c('plot_relatedness_distribution_within',
                                               'plot_frac_highly_related_within',
                                               'plot_frac_highly_related_overtime_within',
                                               'Variable2')
  
  IBD_Transmission_Report_outputs = IBD_Transmission_Report_expected_outputs[IBD_Transmission_Report_expected_outputs %in% ls()]
  
  if(tolower(code_environment) == 'local'){
    imagename = file.path(wd,paste0(output, '_IBD_Transmission_Report.RData'))
  
    save(file = imagename, list = IBD_Transmission_Report_outputs)
  
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_IBD_Transmission_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_IBD_Transmission_Report.Rmd'))))
  
    # Assign variables based on command-line arguments
    render(file.path(wd, paste0(output, '_IBD_Transmission_Report.Rmd')), params = list(
      RData_image = imagename),
      output_dir = wd)
  }else{
  
    imagename = paste0(output, '_IBD_Transmission_Report.RData')
    
    save(file = imagename, list = IBD_Transmission_Report_outputs)
    
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_IBD_Transmission_Report_Template.Rmd'), ' ', paste0(output, '_IBD_Transmission_Report.Rmd')))
    
    # Assign variables based on command-line arguments
    render(paste0(output, '_IBD_Transmission_Report.Rmd'), params = list(
      RData_image = imagename),
      output_dir = 'Results')
  }
  
  print("Leaving render script")
  
}

# COI----

if(!is.null(poly_formula)){
  
  print('Starting COI report')
  
  if(is.null(poly_quantile)){
    poly_quantile = 0.75
  }
  
  if(!is.null(Variable1)){
    
    print('Calculate COI metrics by Variable1')
    
    poly_by_Var1 = get_polygenomic(ampseq_object = ampseq_object, 
                                   strata = Variable1,
                                   update_popsummary = FALSE,
                                   na.rm = na_var_rm,
                                   filters = NULL,
                                   poly_quantile = poly_quantile,
                                   poly_formula = poly_formula)
    
    set.seed(1)
    
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
      scale_fill_manual(values = c(sample(col_vector, nlevels(as.factor(ampseq_object@metadata[[Variable1]]))), "gray30"))+
      theme(axis.text = element_text(size = 12),
            axis.title = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle = 90, vjust = 0.5))
    
  }else{
    
    print('Calculate COI metrics by the overall population')
    
    poly_total = get_polygenomic(ampseq_object = ampseq_object, 
                                 strata = NULL,
                                 update_popsummary = FALSE,
                                 na.rm = na_var_rm,
                                 filters = NULL,
                                 poly_quantile = poly_quantile,
                                 poly_formula = poly_formula)
    
  }
  
  
  
  if(!is.null(Variable1) & !is.null(Variable2)){
    
    print('Calculate COI metrics by Variable1 and Varibale2')
    
    ampseq_object@metadata[['Var1_Var2']] = paste(ampseq_object@metadata[[Variable1]], ampseq_object@metadata[[Variable2]], sep = '::')
    
    poly_by_Var1_Var2 = get_polygenomic(ampseq_object = ampseq_object,
                                        strata = "Var1_Var2",
                                        update_popsummary = F,
                                        na.rm = TRUE,
                                        filters = NULL,
                                        poly_quantile = poly_quantile, 
                                        poly_formula = poly_formula
                                        )
    
    set.seed(1)
    
    plot_poly_by_Var1_Var2 = poly_by_Var1_Var2$pop_summary %>%
      filter(pop != 'Total')%>%
      mutate(
        Variable1 = stringr::str_split(pop, '::', simplify = TRUE)[,1],
        Variable2 = stringr::str_split(pop, '::', simplify = TRUE)[,2],
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
      facet_wrap(~Variable1, ncol = 3)+
      theme_bw()+
      scale_fill_manual(values = sample(col_vector, nlevels(as.factor(ampseq_object@metadata[[Variable1]]))))+
      labs(title = 'Temporal change of the proportion of polyclonal infections',
           y = "Polyclonal infections",
           x = Variable2)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
            legend.position =  "none")
    
  }
  

  print('Generation of plots and tables for COI report done')
  
  COI_Report_expected_outputs = c('poly_total',
                                  'poly_by_Var1',
                                  'plot_poly_by_pop',
                                  'poly_by_Var1_Var2',
                                  'plot_poly_by_Var1_Var2',
                                  'Variable1',
                                  'Variable2')
  
  COI_Report_outputs = COI_Report_expected_outputs[COI_Report_expected_outputs %in% ls()]
  
  if(tolower(code_environment) == 'local'){
  
    imagename = file.path(wd,paste0(output, '_COI_Report.RData'))
  
    save(file = imagename, list = COI_Report_outputs)
  
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_COI_Report_Template.Rmd'), ' ', file.path(wd, paste0(output, '_COI_Report.Rmd'))))
  
    # Assign variables based on command-line arguments
    render(file.path(wd, paste0(output, '_COI_Report.Rmd')), params = list(
      RData_image = imagename),
      output_dir = wd)
  }else{
  
    imagename = paste0(output, '_COI_Report.RData')
    
    save(file = imagename, list = COI_Report_outputs)
    
    system(paste0('cp ', file.path(fd, 'MHap_Analysis_COI_Report_Template.Rmd'), ' ', paste0(output, '_COI_Report.Rmd')))
    
    # Assign variables based on command-line arguments
    render(paste0(output, '_COI_Report.Rmd'), params = list(
      RData_image = imagename),
      output_dir = 'Results')
  }
  
  print("Leaving render script")
  
}



# Genetic diversity----



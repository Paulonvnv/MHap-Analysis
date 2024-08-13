
# PATHs to references and functions----
wd = gsub("/$", "", "/Users/pam3650/Documents/Github/MHap-Analysis/docs/data/Pfal_example3/")
fd = gsub("/$", "", "/Users/pam3650/Documents/Github/MHap-Analysis/docs/functions_and_libraries/")
rd = gsub("/$", "", "/Users/pam3650/Documents/Github/MHap-Analysis/docs/reference/Pfal_3D7/")
## Upload of raw data----
cigar_paths = "null"

cigar_files = "cigar_tables"
asv_table_files = "asv_tables"
asv2cigar_files = "asv2cigar"
asv_seq_files = "asv_seqs"
zero_read_sample_list = "zeroReadSamples"

ampseq_jsonfile = "null"
ampseq_excelfile = "null"
ampseq_csvfolder = "null"

sample_id_pattern = "^SP"
markers = "markers.csv"

## Parameters to export data and results----

output = "MHap_test_excel"
ampseq_export_format = "excel"

## Removing or masking ASVs in the genotype table----

min_abd = as.integer(10)
min_ratio = as.numeric(0.1)

off_target_formula = as.character("dVSITES_ij>=0.3")
flanking_INDEL_formula = as.character("flanking_INDEL==TRUE&h_ij>=0.66")

homopolymer_length = as.integer(5)
SNV_in_homopolymer_formula = as.character("SNV_in_homopolymer==TRUE&h_ij>=0.66")
INDEL_in_homopolymer_formula = as.character("INDEL_in_homopolymer==TRUE&h_ij>=0.66")

bimera_formula = as.character("bimera==TRUE&h_ij>=0.66")
PCR_errors_formula = as.character("h_ij>=0.66&h_ijminor>=0.66&p_ij>=0.05")
sample_ampl_rate = as.numeric(0.75)
locus_ampl_rate  = as.numeric(0.75)

## Adding metadata----

metadata_file = as.character("colombia_metadata_cleaned.csv")
join_by = as.character("Sample_id")
Variable1 = "Subnational_level2"
Variable2 = "Quarter_of_Collection"
Longitude = "null"
Latitude = "null"
var_filter  = "Subnational_level2;Keep;Buenaventura,Quibdó,Guapi/Quarter_of_Collection;Keep;2021-Q1,2021-Q2,2021-Q3,2021-Q4,2022-Q1,2022-Q2"

## Reports----

PerformanceReport = as.logical("true")
Drug_Surveillance_Report = as.logical("false")
Variants_of_Interest_Report = as.logical("false")
ibd_thres = "null"
poly_formula = as.character("null")

### Parameters for DSR----

ref_gff = "PlasmoDB-59_Pfalciparum3D7.gff"
ref_fasta = "PlasmoDB-59_Pfalciparum3D7_Genome.fasta"
amplicon_fasta = "pf3d7_ref_updated_v3.fasta"
reference_alleles = "drugR_alleles.csv"
hap_color_palette = "random"
gene_names = "PfDHFR,PfMDR1,PfDHPS,PfKelch13,PF3D7_1447900"
gene_ids = "PF3D7_0417200,PF3D7_0523000,PF3D7_0810800,PF3D7_1343700,PF3D7_1447900"
drugs = "Artemisinin,Chloroquine,Pyrimethamine,Sulfadoxine,Lumefantrine,Mefloquine"
include_all_drug_markers = as.logical("true")
na_var_rm = as.logical("true")
na_hap_rm = as.logical("true")

### Parameters for IBD and Conectivity report----

#### Parameters for UGER----
nTasks = as.integer(1)
Task_id = as.integer(1)
ibd_step = as.character("null")

#### Other IBD parameters----

pairwise_relatedness_table = as.character("null")
parallel = as.logical("true")
nchunks = as.integer(500)
ibd_ncol = as.integer(4)
pop_levels = as.integer("null")

### Parameters for COI report----
# poly_quantile
poly_quantile = 0.75

# PATHs to references and functions----

print(paste0('wd: ', wd))
setwd(wd) 

# Tools or functions directory
print(paste0('fd: ', fd))

# Reference files directory
print(paste0('rd: ', rd))

## Upload of raw data----
# Starting cigar file
cigar_paths = if(cigar_paths == 'null'){NULL}else{file.path(wd, cigar_paths)}
print(paste0('cigar_paths: ', cigar_paths))

# Starting cigar file
cigar_files = if(cigar_files == 'null'){NULL}else{file.path(wd, cigar_files)}
print(paste0('cigar_files: ', cigar_files))

# Starting cigar file
asv_table_files = if(asv_table_files == 'null'){NULL}else{file.path(wd, asv_table_files)}
print(paste0('asv_table_files: ', cigar_files))

# Starting cigar file
asv2cigar_files = if(asv2cigar_files == 'null'){NULL}else{file.path(wd, asv2cigar_files)}
print(paste0('asv2cigar_files: ', asv2cigar_files))

# Starting cigar file
asv_seq_files = if(asv_seq_files == 'null'){NULL}else{file.path(wd, asv_seq_files)}
print(paste0('asv_seq_files: ', asv_seq_files))

# Starting cigar file
zero_read_sample_list = if(zero_read_sample_list == 'null'){NULL}else{file.path(wd, zero_read_sample_list)}
print(paste0('zero_read_sample_list: ', zero_read_sample_list))


# Starting ampseq file in json format
ampseq_jsonfile = if(ampseq_jsonfile == 'null'){NULL}else{file.path(wd, ampseq_jsonfile)}
print(paste0('ampseq_jsonfile: ', ampseq_jsonfile))

# Starting ampseq file in excel format
ampseq_excelfile = if(ampseq_excelfile == 'null'){NULL}else{file.path(wd, args$ampseq_excelfile)}
print(paste0('ampseq_excelfile: ', ampseq_excelfile))

# Starting ampseq file in csv format
ampseq_csvfolder = if(ampseq_csvfolder == 'null'){NULL}else{file.path(wd, args$ampseq_csvfolder)}
print(paste0('ampseq_csvfolder: ', ampseq_csvfolder))

# Pattern to differentiate between samples of interest and controls
print(paste0('sample_id_pattern: ', sample_id_pattern))

# csv table with markers information
markers = if(markers == 'null'){NULL}else{file.path(rd, markers)}
print(paste0('markers: ', markers))

## Parameters to export data and results----

# output pattern
print(paste0('output: ', output))

# format to export the ampseq object
ampseq_export_format = if(ampseq_export_format == 'null'){NULL}else{as.character(ampseq_export_format)}

## Removing or masking ASVs in the genotype table----

# Minimum abundance to call an allele
print(paste0('min_abd: ', min_abd))

# Minimum ratio to call minor alleles

print(paste0('min_ratio: ', min_ratio))

# off_target_formula filter
off_target_formula = if(off_target_formula == 'null'){NULL}else{as.character(off_target_formula)}

if(!is.null(off_target_formula)){
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
  
}

print(paste0('off_target_formula: ', off_target_formula))

# flanking_INDEL_formula filter
flanking_INDEL_formula = if(flanking_INDEL_formula == 'null'){NULL}else{as.character(flanking_INDEL_formula)}

if(!is.null(flanking_INDEL_formula)){
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
}

print(paste0('flanking_INDEL_formula: ', flanking_INDEL_formula))

# homopolymer_length
print(paste0('homopolymer_length: ', homopolymer_length))

# SNV_in_homopolymer_formula filter
SNV_in_homopolymer_formula = if(SNV_in_homopolymer_formula == 'null'){NULL}else{as.character(SNV_in_homopolymer_formula)}

if(!is.null(SNV_in_homopolymer_formula)){
  SNV_in_homopolymer_formula = gsub('"',"",SNV_in_homopolymer_formula)
  SNV_in_homopolymer_formula = gsub('&'," & ",SNV_in_homopolymer_formula, ignore.case = TRUE)
  SNV_in_homopolymer_formula = gsub('\\|'," \\| ",SNV_in_homopolymer_formula, ignore.case = TRUE)
  
  if(grepl("\\w>\\d",SNV_in_homopolymer_formula)){
    patterns = str_extract_all(SNV_in_homopolymer_formula, "\\w>\\d")[[1]]
    
    for(pattern in patterns){
      
      replacement = gsub('>',' > ',pattern)
      SNV_in_homopolymer_formula = gsub(pattern,
                                        replacement,
                                        SNV_in_homopolymer_formula, ignore.case = TRUE)
    }
    
  }
  
  if(grepl("\\w<\\d",SNV_in_homopolymer_formula)){
    patterns = str_extract_all(SNV_in_homopolymer_formula, "\\w<\\d")[[1]]
    
    for(pattern in patterns){
      
      replacement = gsub('<',' < ',pattern)
      SNV_in_homopolymer_formula = gsub(pattern,
                                        replacement,
                                        SNV_in_homopolymer_formula, ignore.case = TRUE)
    }
    
  }
  
  SNV_in_homopolymer_formula = gsub('>='," >= ",SNV_in_homopolymer_formula, ignore.case = TRUE)
  SNV_in_homopolymer_formula = gsub('<='," <= ",SNV_in_homopolymer_formula, ignore.case = TRUE)
  SNV_in_homopolymer_formula = gsub('=='," == ",SNV_in_homopolymer_formula, ignore.case = TRUE)
  SNV_in_homopolymer_formula = gsub('!='," != ",SNV_in_homopolymer_formula, ignore.case = TRUE)
  
  SNV_in_homopolymer_formula = gsub('\\+'," \\+ ",SNV_in_homopolymer_formula, ignore.case = TRUE)
  SNV_in_homopolymer_formula = gsub('-'," - ",SNV_in_homopolymer_formula, ignore.case = TRUE)
  SNV_in_homopolymer_formula = gsub('\\*'," \\* ",SNV_in_homopolymer_formula, ignore.case = TRUE)
  SNV_in_homopolymer_formula = gsub('/'," / ",SNV_in_homopolymer_formula, ignore.case = TRUE)
}

print(paste0('SNV_in_homopolymer_formula: ', SNV_in_homopolymer_formula))

# INDEL_in_homopolymer_formula filter
INDEL_in_homopolymer_formula = if(INDEL_in_homopolymer_formula == 'null'){NULL}else{as.character(INDEL_in_homopolymer_formula)}

if(!is.null(INDEL_in_homopolymer_formula)){
  INDEL_in_homopolymer_formula = gsub('"',"",INDEL_in_homopolymer_formula)
  INDEL_in_homopolymer_formula = gsub('&'," & ",INDEL_in_homopolymer_formula, ignore.case = TRUE)
  INDEL_in_homopolymer_formula = gsub('\\|'," \\| ",INDEL_in_homopolymer_formula, ignore.case = TRUE)
  
  if(grepl("\\w>\\d",INDEL_in_homopolymer_formula)){
    patterns = str_extract_all(INDEL_in_homopolymer_formula, "\\w>\\d")[[1]]
    
    for(pattern in patterns){
      
      replacement = gsub('>',' > ',pattern)
      INDEL_in_homopolymer_formula = gsub(pattern,
                                          replacement,
                                          INDEL_in_homopolymer_formula, ignore.case = TRUE)
    }
    
  }
  
  if(grepl("\\w<\\d",INDEL_in_homopolymer_formula)){
    patterns = str_extract_all(INDEL_in_homopolymer_formula, "\\w<\\d")[[1]]
    
    for(pattern in patterns){
      
      replacement = gsub('<',' < ',pattern)
      INDEL_in_homopolymer_formula = gsub(pattern,
                                          replacement,
                                          INDEL_in_homopolymer_formula, ignore.case = TRUE)
    }
    
  }
  
  INDEL_in_homopolymer_formula = gsub('>='," >= ",INDEL_in_homopolymer_formula, ignore.case = TRUE)
  INDEL_in_homopolymer_formula = gsub('<='," <= ",INDEL_in_homopolymer_formula, ignore.case = TRUE)
  INDEL_in_homopolymer_formula = gsub('=='," == ",INDEL_in_homopolymer_formula, ignore.case = TRUE)
  INDEL_in_homopolymer_formula = gsub('!='," != ",INDEL_in_homopolymer_formula, ignore.case = TRUE)
  
  INDEL_in_homopolymer_formula = gsub('\\+'," \\+ ",INDEL_in_homopolymer_formula, ignore.case = TRUE)
  INDEL_in_homopolymer_formula = gsub('-'," - ",INDEL_in_homopolymer_formula, ignore.case = TRUE)
  INDEL_in_homopolymer_formula = gsub('\\*'," \\* ",INDEL_in_homopolymer_formula, ignore.case = TRUE)
  INDEL_in_homopolymer_formula = gsub('/'," / ",INDEL_in_homopolymer_formula, ignore.case = TRUE)
}

print(paste0('INDEL_in_homopolymer_formula: ', INDEL_in_homopolymer_formula))


# bimera_formula filter
bimera_formula = if(bimera_formula == 'null'){NULL}else{as.character(bimera_formula)}

if(!is.null(bimera_formula)){
  bimera_formula = gsub('"',"",bimera_formula)
  bimera_formula = gsub('&'," & ",bimera_formula, ignore.case = TRUE)
  bimera_formula = gsub('\\|'," \\| ",bimera_formula, ignore.case = TRUE)
  
  if(grepl("\\w>\\d",bimera_formula)){
    patterns = str_extract_all(bimera_formula, "\\w>\\d")[[1]]
    
    for(pattern in patterns){
      
      replacement = gsub('>',' > ',pattern)
      bimera_formula = gsub(pattern,
                            replacement,
                            bimera_formula, ignore.case = TRUE)
    }
    
  }
  
  if(grepl("\\w<\\d",bimera_formula)){
    patterns = str_extract_all(bimera_formula, "\\w<\\d")[[1]]
    
    for(pattern in patterns){
      
      replacement = gsub('<',' < ',pattern)
      bimera_formula = gsub(pattern,
                            replacement,
                            bimera_formula, ignore.case = TRUE)
    }
    
  }
  
  bimera_formula = gsub('>='," >= ", bimera_formula, ignore.case = TRUE)
  bimera_formula = gsub('<='," <= ", bimera_formula, ignore.case = TRUE)
  bimera_formula = gsub('=='," == ", bimera_formula, ignore.case = TRUE)
  bimera_formula = gsub('!='," != ", bimera_formula, ignore.case = TRUE)
  
  bimera_formula = gsub('\\+'," \\+ ", bimera_formula, ignore.case = TRUE)
  bimera_formula = gsub('-'," - ", bimera_formula, ignore.case = TRUE)
  bimera_formula = gsub('\\*'," \\* ", bimera_formula, ignore.case = TRUE)
  bimera_formula = gsub('/'," / ", bimera_formula, ignore.case = TRUE)
}

print(paste0('bimera_formula: ', bimera_formula))

# PCR_errors_formula filter
PCR_errors_formula = if(PCR_errors_formula == 'null'){NULL}else{as.character(PCR_errors_formula)}

if(!is.null(PCR_errors_formula)){
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
}

print(paste0('PCR_errors_formula: ', PCR_errors_formula))

## Removing samples or loci with low performance----

# sample_ampl_rate
print(paste0('sample_ampl_rate: ', sample_ampl_rate))

# locus_ampl_rate
print(paste0('locus_ampl_rate: ', locus_ampl_rate))

## Adding metadata----

# metadata
metadata_file = if(metadata_file == 'null'){NULL}else{file.path(wd, metadata_file)}
print(paste0('metadata: ', metadata_file))

# join metadata by
join_by = if(join_by == 'null'){NULL}else{join_by}
print(paste0('join_by: ', join_by))

# pop
Variable1 = if(Variable1 == 'null'){NULL}else{Variable1}
print(paste0('Variable1: ', Variable1))

# temporal_population
Variable2 = if(Variable2 == 'null'){NULL}else{Variable2}
print(paste0('Variable2: ', Variable2))

Longitude = if(Longitude == 'null'){NULL}else{Longitude}
print(paste0('Longitude: ', Longitude))

Latitude = if(Latitude == 'null'){NULL}else{Latitude}
print(paste0('Latitude: ', Latitude))

## Filtering desired or undesired populations or samples----

if(var_filter == 'null'){
  var_filter = NULL
}else{
  var_filter = gsub('\\.', ' ', var_filter)
  var_filter = strsplit(var_filter, '/')[[1]]
}

print('var_filter: ')
print(var_filter)

## Reports----

# Performance Report
print(paste0('PerformanceReport: ', PerformanceReport))

# Drug surveillance report
print(paste0('Drug_Surveillance_Report: ', Drug_Surveillance_Report))

# Variants of interest report (in development)
print(paste0('Variants_of_Interest_Report: ', Variants_of_Interest_Report))

# Reports of IBD and Connectivity and IBD and transmission
ibd_thres = if(ibd_thres == 'null'){NULL}else{as.numeric(ibd_thres)}
print(paste0('ibd_thres: ', ibd_thres))

# Report for complexity of infection

# poly_formula filter
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

ref_gff = if(ref_gff == 'null'){NULL}else{file.path(rd, ref_gff)}
print(paste0('ref_gff: ', ref_gff))

ref_fasta = if(ref_fasta == 'null'){NULL}else{file.path(rd, ref_fasta)}
print(paste0('ref_fasta: ', ref_fasta))

amplicon_fasta = if(amplicon_fasta == 'null'){NULL}else{file.path(rd, amplicon_fasta)}
print(paste0('amplicon_fasta: ', amplicon_fasta))

reference_alleles = if(reference_alleles == 'null'){NULL}else{file.path(rd, reference_alleles)}
print(paste0('reference_alleles: ', reference_alleles))

print(paste0('hap_color_palette: ', hap_color_palette))

# gene_names
if(gene_names == 'null'){
  gene_names = NULL
}else{
  gene_names = strsplit(gene_names, ',')[[1]]
}
print('gene_names:')
print(gene_names)

# gene_ids
if(gene_ids == 'null'){
  gene_ids = NULL
}else{
  gene_ids = strsplit(gene_ids, ',')[[1]]
}
print('gene_ids:')
print(gene_ids)

if(drugs == 'null'){
  drugs = NULL
}else{
  drugs = strsplit(drugs, ',')[[1]]
}
print('drugs: ')
print(drugs)

print(paste0('include_all_drug_markers: ', include_all_drug_markers))

print(paste0('na_var_rm: ', na_var_rm))

print(paste0('na_hap_rm: ', na_hap_rm))


### Parameters for IBD and Conectivity report----

#### Parameters for UGER----
# nTasks
print(paste0('nTasks: ', nTasks))

# Task_id
print(paste0('Task_id: ', Task_id))

# ibd_step
print(paste0('ibd_step: ', ibd_step))

#### Other IBD parameters----

# pairwise_relatedness_table
pairwise_relatedness_table = if(pairwise_relatedness_table == 'null'){NULL}else{file.path(wd, pairwise_relatedness_table)}
print(paste0('pairwise_relatedness_table: ', pairwise_relatedness_table))

# parallel
print(paste0('parallel: ', parallel))

# nchunks
print(paste0('nchunks: ', nchunks))

# ibd_ncol
print(paste0('ibd_ncol: ', ibd_ncol))

# pop_levels
print(paste0('pop_levels: ', pop_levels))

### Parameters for COI report----
# poly_quantile
poly_quantile = if(poly_quantile == 'null'){NULL}else{as.numeric(poly_quantile)}
print(paste0('poly_quantile: ', poly_quantile))


# The next five lines are needed for UGER
# if(file.exists(file.path(wd, paste0(output, '.xlsx')))){
#   cigar_paths = NULL
#   cigar_files = NULL
#   ampseq_excelfile = file.path(wd, paste0(output, '.xlsx'))
# }

print("All variables checked")
args = list(NULL)


args$wd = '~/Documents/Github/MADHATTER_Senegal/Senegal_test3/cromwell-executions/ampseq/2f726715-08dc-48b7-80a4-fcf5c5f843db/call-t_005_asv_filtering/execution/Results/'
args$fd = '~/Documents/Github/malaria/WDL/Tasks/ASV_Filtering/malaria_asvfiltering/Code/'
args$rd = '~/Documents/Github/MADHATTER_Senegal/Senegal_test3/cromwell-executions/ampseq/2f726715-08dc-48b7-80a4-fcf5c5f843db/call-t_005_asv_filtering/execution/references/'
args$cigar_paths = 'null'
args$cigar_files = 'cigar_variants'
args$asv_table_files = 'asv_tables'
args$asv2cigar_files = 'asv2cigar'
args$asv_seq_files = 'asv_seq'
args$zero_read_sample_list = 'zeroReadSampleList'
args$ampseq_jsonfile = 'null'
args$ampseq_excelfile = 'null'
args$ampseq_csvfolder = 'null'
args$sample_id_pattern = '.'
args$markers = 'markersTable.csv'
args$out = 'out_prefix'
args$ampseq_export_format = 'excel'
args$min_abd = '10'
args$min_ratio = '0.1'
args$off_target_formula = 'dVSITES_ij>=0.2'
args$flanking_INDEL_formula = 'flanking_INDEL==TRUE'
args$homopolymer_length = '5'
args$SNV_in_homopolymer_formula = 'SNV_in_homopolymer==TRUE'
args$INDEL_in_homopolymer_formula = 'INDEL_in_homopolymer==TRUE'
args$bimera_formula = 'bimera==TRUE&h_ij>=0.66'
args$PCR_errors_formula = 'h_ij>=0.66&h_ijminor>=0.66&p_ij>=0.05'
args$sample_ampl_rate = '0.3'
args$locus_ampl_rate = '0.3'
args$metadata = 'metadata.csv'
args$join_by = 'Sample_id'
args$Variable1 = 'Country'
args$Variable2 = 'Country'
args$Longitude = 'Longitude'
args$Latitude = 'Latitude'
args$var_filter = 'null'
args$PerformanceReport = 'false'
args$Drug_Surveillance_Report = 'false'
args$Variants_of_Interest_Report = 'false'
args$ibd_thres = 'null'
args$poly_formula = 'null'
args$ref_gff = 'null'
args$ref_fasta = 'PkPfPmPoPv.fasta'
args$amplicon_fasta = 'reference.fasta'
args$reference_alleles = 'null'
args$hap_color_palette = 'null'
args$gene_names = 'null'
args$gene_ids = 'null'
args$drugs = 'null'
args$include_all_drug_markers = 'null'
args$na_var_rm = 'false'
args$na_hap_rm = 'false'
args$nTasks = '1'
args$Task_id = '1'
args$ibd_step = 'null'
args$pairwise_relatedness_table = 'null'
args$parallel = 'false'
args$nchunks = '1'
args$ibd_ncol = '4'
args$pop_levels = 'null'
args$poly_quantile = 'null'

print("starting to parse variables")


# 2.Paths to references and functions----

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

# 2.Upload of raw data----
# Upload cigar files, asv_tables and zero_read_files
cigar_paths = args$cigar_paths
cigar_paths = if(cigar_paths == 'null' || cigar_paths == ''){NULL}else{file.path(wd, cigar_paths)}
print(paste0('cigar_paths: ', cigar_paths))

cigar_files = args$cigar_files
cigar_files = if(cigar_files == 'null' || cigar_files == ''){NULL}else{file.path(wd, cigar_files)}
print(paste0('cigar_files: ', cigar_files))

asv_table_files = args$asv_table_files
asv_table_files = if(asv_table_files == 'null' || asv_table_files == ''){NULL}else{file.path(wd, asv_table_files)}
print(paste0('asv_table_files: ', asv_table_files))

asv2cigar_files = args$asv2cigar_files
asv2cigar_files = if(asv2cigar_files == 'null' || asv2cigar_files == ''){NULL}else{file.path(wd, asv2cigar_files)}
print(paste0('asv2cigar_files: ', asv2cigar_files))

asv_seq_files = args$asv_seq_files
asv_seq_files = if(asv_seq_files == 'null' || asv_seq_files == ''){NULL}else{file.path(wd, asv_seq_files)}
print(paste0('asv_seq_files: ', asv_seq_files))

zero_read_sample_list = args$zero_read_sample_list
zero_read_sample_list = if(zero_read_sample_list == 'null' || zero_read_sample_list == ''){NULL}else{file.path(wd, zero_read_sample_list)}
print(paste0('zero_read_sample_list: ', zero_read_sample_list))

# Starting ampseq file in json format
ampseq_jsonfile = args$ampseq_jsonfile
ampseq_jsonfile = if(ampseq_jsonfile == 'null' || ampseq_json_file == ''){NULL}else{file.path(wd, ampseq_jsonfile)}
print(paste0('ampseq_jsonfile: ', ampseq_jsonfile))

# Starting ampseq file in excel format
ampseq_excelfile = args$ampseq_excelfile
ampseq_excelfile = if(ampseq_excelfile == 'null' || ampseq_excelfile == ''){NULL}else{file.path(wd, args$ampseq_excelfile)}
print(paste0('ampseq_excelfile: ', ampseq_excelfile))

# Starting ampseq file in csv format
ampseq_csvfolder = args$ampseq_csvfolder
ampseq_csvfolder = if(ampseq_csvfolder == 'null' || ampseq_csvfolder == ''){NULL}else{file.path(wd, args$ampseq_csvfolder)}
print(paste0('ampseq_csvfolder: ', ampseq_csvfolder))

# Pattern to differentiate between samples of interest and controls
sample_id_pattern = args$sample_id_pattern
print(paste0('sample_id_pattern: ', sample_id_pattern))

# csv table with markers information
markers = args$markers
markers = if(markers == 'null' || markers == ''){NULL}else{file.path(rd, markers)}
print(paste0('markers: ', markers))

# 2.Parameters to export data and results----

# output pattern
output = args$out
print(paste0('output: ', output))

# format to export the ampseq object
ampseq_export_format = args$ampseq_export_format
ampseq_export_format = if(ampseq_export_format == 'null'){NULL}else{as.character(ampseq_export_format)}

# 2.Removing or masking ASVs in the genotype table----

# Minimum abundance to call an allele
min_abd = as.integer(args$min_abd)
print(paste0('min_abd: ', min_abd))

# Minimum ratio to call minor alleles
min_ratio = as.numeric(args$min_ratio)
print(paste0('min_ratio: ', min_ratio))

# 2.off_target_formula filter ----
off_target_formula = as.character(args$off_target_formula)
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

# 2.flanking_INDEL_formula filter ----
flanking_INDEL_formula = as.character(args$flanking_INDEL_formula)
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

# 2.homopolymer filter parameters ----
# homopolymer_length
homopolymer_length = as.integer(args$homopolymer_length)
print(paste0('homopolymer_length: ', homopolymer_length))

# SNV_in_homopolymer_formula filter
SNV_in_homopolymer_formula = as.character(args$SNV_in_homopolymer_formula)
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
INDEL_in_homopolymer_formula = as.character(args$INDEL_in_homopolymer_formula)
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

# 2.bimera_formula filter ----
bimera_formula = as.character(args$bimera_formula)
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

# 2.PCR_errors_formula filter ----
PCR_errors_formula = as.character(args$PCR_errors_formula)
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

# 2.Removing samples or loci with low performance----

# sample_ampl_rate
sample_ampl_rate = as.numeric(args$sample_ampl_rate)
print(paste0('sample_ampl_rate: ', sample_ampl_rate))

# locus_ampl_rate
locus_ampl_rate  = as.numeric(args$locus_ampl_rate)
print(paste0('locus_ampl_rate: ', locus_ampl_rate))

# 2.Adding metadata----

# metadata
metadata_file = as.character(args$metadata)
metadata_file = if(metadata_file == 'null'){NULL}else{file.path(wd, metadata_file)}
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

# 2.Filtering populations or samples----

var_filter  = args$var_filter
if(var_filter == 'null'){
  var_filter = NULL
}else{
  var_filter = gsub('\\.', ' ', var_filter)
  var_filter = strsplit(var_filter, '/')[[1]]
}

print('var_filter: ')
print(var_filter)

# 2.Define Reports----

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

# 2.Parameters for DSR----
# Reference gff and fasta files, require to translate DNA cigar formats to aminoacid sequences
ref_gff = args$ref_gff
ref_gff = if(ref_gff == 'null'){NULL}else{file.path(rd, ref_gff)}
print(paste0('ref_gff: ', ref_gff))

ref_fasta = args$ref_fasta
ref_fasta = if(ref_fasta == 'null'){NULL}else{file.path(rd, ref_fasta)}
print(paste0('ref_fasta: ', ref_fasta))

amplicon_fasta = args$amplicon_fasta
amplicon_fasta = if(amplicon_fasta == 'null'){NULL}else{file.path(rd, amplicon_fasta)}
print(paste0('amplicon_fasta: ', amplicon_fasta))

reference_alleles = args$reference_alleles
reference_alleles = if(reference_alleles == 'null'){NULL}else{file.path(rd, reference_alleles)}
print(paste0('reference_alleles: ', reference_alleles))

hap_color_palette = args$hap_color_palette
print(paste0('hap_color_palette: ', hap_color_palette))

# gene_names
gene_names = args$gene_names
if(gene_names == 'null'){
  gene_names = NULL
}else{
  gene_names = strsplit(gene_names, ',')[[1]]
}
print('gene_names:')
print(gene_names)

# gene_ids
gene_ids = args$gene_ids
if(gene_ids == 'null'){
  gene_ids = NULL
}else{
  gene_ids = strsplit(gene_ids, ',')[[1]]
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

# 2.Parameters for IBD and Connectivity report----

# Parameters for UGER
# nTasks
nTasks = as.integer(args$nTasks)
print(paste0('nTasks: ', nTasks))

# Task_id
Task_id = as.integer(args$Task_id)
print(paste0('Task_id: ', Task_id))

# ibd_step
ibd_step = as.character(args$ibd_step)
print(paste0('ibd_step: ', ibd_step))

# Other IBD parameters

# pairwise_relatedness_table
pairwise_relatedness_table = as.character(args$pairwise_relatedness_table)
pairwise_relatedness_table = if(pairwise_relatedness_table == 'null'){NULL}else{file.path(wd, pairwise_relatedness_table)}
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
pop_levels = as.integer(args$pop_levels)
print(paste0('pop_levels: ', pop_levels))

# 2.Parameters for COI report----
# poly_quantile
poly_quantile = args$poly_quantile
poly_quantile = if(poly_quantile == 'null'){NULL}else{as.numeric(poly_quantile)}
print(paste0('poly_quantile: ', poly_quantile))

# The next five lines are needed for UGER
# if(file.exists(file.path(wd, paste0(output, '.xlsx')))){
#   cigar_paths = NULL
#   cigar_files = NULL
#   ampseq_excelfile = file.path(wd, paste0(output, '.xlsx'))
# }

print("All variables checked")

# ampseq_object@asv_table$Amplicon
# 
# ampseq_object@asv_table$CIGAR_masked
# 
# test = cigar_strings2fasta(obj = ampseq_object@asv_table, ref_fasta = amplicon_fasta, cigar_string_col = 'CIGAR_masked', amplicon_col = 'Amplicon', format = 'DNAStringSet')


# # 
# if(length(slot(ampseq_object, 'asv_seqs')) != length(slot(ampseq_object, 'asv_seqs_masked'))){
#   stop('masked and unmasked fastas with different length')
# }

# if(length(slot(ampseq_object, 'asv_seqs')) != length(slot(ampseq_object, 'asv_seqs_masked'))){
#   stop('masked and unmasked fastas with different length')
# }
# 
# if(sum(names(slot(ampseq_object, 'asv_seqs')) != names(slot(ampseq_object, 'asv_seqs_masked'))) > 0){
#   stop('masked and unmasked fastas with different length')
# }









#!/bin/bash
#source /broad/software/scripts/useuse
#use R-4.4

function json_extract() {
  # Citation: https://stackoverflow.com/questions/1955505/parsing-json-with-unix-tools
  local key=$1
  local json=$2
  local string_regex='"([^"\]|\.)*"'
  local string_logial='([A-Za-z]+)*'
  local number_regex='-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][+-]?[0-9]+)?'
  local time_regex='[0-9][0-9]+:(([0-5][0-9])|60):(([0-5][0-9])|60)'
  local value_regex="${string_regex}|${number_regex}|${time_regex}|${string_logial}"
  local pair_regex="\"${key}\"[[:space:]]*:[[:space:]]*(${value_regex})"

  if [[ ${json} =~ ${pair_regex} ]]
  then
	echo $(sed 's/^"\|"$//g' <<< "${BASH_REMATCH[1]}")
  else
	return 1
  fi
	}

# Extract variables
json=$1

# PATHs to references and functions
wd=$(json_extract wd "$(cat ${json})")
wd=${wd%\"}
wd=${wd#\"}
echo "wd: "${wd}

fd=$(json_extract fd "$(cat ${json})")
fd=${fd%\"}
fd=${fd#\"}
echo "fd: "${fd}

rd=$(json_extract rd "$(cat ${json})")
rd=${rd%\"}
rd=${rd#\"}
echo "rd: "${rd}

# Upload of raw data

cigar_paths=$(json_extract cigar_paths "$(cat ${json})")
if [[ ${cigar_paths} != "" ]]
then 
cigar_paths=${cigar_paths%\"}
cigar_paths=${cigar_paths#\"}
else
cigar_paths="null"
fi
echo "cigar_paths: "${cigar_paths}

cigar_files=$(json_extract cigar_files "$(cat ${json})")
if [[ ${cigar_files} != "" ]]
then 
cigar_files=${cigar_files%\"}
cigar_files=${cigar_files#\"}
else
cigar_files="null"
fi
echo "cigar_files: "${cigar_files}


asv_table_files=$(json_extract asv_table_files "$(cat ${json})")
if [[ ${asv_table_files} != "" ]]
then 
asv_table_files=${asv_table_files%\"}
asv_table_files=${asv_table_files#\"}
else
asv_table_files="null"
fi
echo "asv_table_files: "${asv_table_files}


asv2cigar_files=$(json_extract asv2cigar_files "$(cat ${json})")
if [[ ${asv2cigar_files} != "" ]]
then 
asv2cigar_files=${asv2cigar_files%\"}
asv2cigar_files=${asv2cigar_files#\"}
else
asv2cigar_files="null"
fi
echo "asv2cigar_files: "${asv2cigar_files}


asv_seq_files=$(json_extract asv_seq_files "$(cat ${json})")
if [[ ${asv_seq_files} != "" ]]
then 
asv_seq_files=${asv_seq_files%\"}
asv_seq_files=${asv_seq_files#\"}
else
asv_seq_files="null"
fi
echo "asv_seq_files: "${asv_seq_files}


zero_read_sample_list=$(json_extract zero_read_sample_list "$(cat ${json})")
if [[ ${zero_read_sample_list} != "" ]]
then 
zero_read_sample_list=${zero_read_sample_list%\"}
zero_read_sample_list=${zero_read_sample_list#\"}
else
zero_read_sample_list="null"
fi
echo "zero_read_sample_list: "${zero_read_sample_list}


ampseq_jsonfile=$(json_extract ampseq_jsonfile "$(cat ${json})")
if [[ ${ampseq_jsonfile} != "" ]]
then 
ampseq_jsonfile=${ampseq_jsonfile%\"}
ampseq_jsonfile=${ampseq_jsonfile#\"}
else
ampseq_jsonfile="null"
fi
echo "ampseq_jsonfile: "${ampseq_jsonfile}


ampseq_excelfile=$(json_extract ampseq_excelfile "$(cat ${json})")
if [[ ${ampseq_excelfile} != "" ]]
then 
ampseq_excelfile=${ampseq_excelfile%\"}
ampseq_excelfile=${ampseq_excelfile#\"}
else
ampseq_excelfile="null"
fi
echo "ampseq_excelfile: "${ampseq_excelfile}


ampseq_csvfolder=$(json_extract ampseq_csvfolder "$(cat ${json})")
if [[ ${ampseq_csvfolder} != "" ]]
then 
ampseq_csvfolder=${ampseq_csvfolder%\"}
ampseq_csvfolder=${ampseq_csvfolder#\"}
else
ampseq_csvfolder="null"
fi
echo "ampseq_csvfolder: "${ampseq_csvfolder}

sample_id_pattern=$(json_extract sample_id_pattern "$(cat ${json})")
if [[ ${sample_id_pattern} != "" ]]
then 
sample_id_pattern=${sample_id_pattern%\"}
sample_id_pattern=${sample_id_pattern#\"}
else
sample_id_pattern="."
fi
echo "sample_id_pattern: "${sample_id_pattern}


markers=$(json_extract markers "$(cat ${json})")
if [[ ${markers} != "" ]]
then 
markers=${markers%\"}
markers=${markers#\"}
else
markers="null"
fi
echo "markers: "${markers}

# Parameters to export data and results

output=$(json_extract output "$(cat ${json})")
output=${output%\"}
output=${output#\"}
echo "output: "${output}

ampseq_export_format=$(json_extract ampseq_export_format "$(cat ${json})")
ampseq_export_format=${ampseq_export_format%\"}
ampseq_export_format=${ampseq_export_format#\"}
echo "ampseq_export_format: "${ampseq_export_format}

# Removing or masking ASVs in the genotype table

min_abd=$(json_extract min_abd "$(cat ${json})")
if [[ ${min_abd} == "" ]]
then 
min_abd=10
fi
echo "min_abd: "${min_abd}

min_ratio=$(json_extract min_ratio "$(cat ${json})")
if [[ ${min_ratio} == "" ]]
then 
min_ratio=0.1
fi
echo "min_ratio: "${min_ratio}


off_target_formula=$(json_extract off_target_formula "$(cat ${json})")
if [[ ${off_target_formula} == "" ]]
then 
off_target_formula="dVSITES_ij>=0.3"
fi
echo "off_target_formula: "${off_target_formula}

flanking_INDEL_formula=$(json_extract flanking_INDEL_formula "$(cat ${json})")
if [[ ${flanking_INDEL_formula} == "" ]]
then 
flanking_INDEL_formula="flanking_INDEL==TRUE&h_ij>=0.66"
fi
echo "flanking_INDEL_formula: "${flanking_INDEL_formula}

homopolymer_length=$(json_extract homopolymer_length "$(cat ${json})")
if [[ ${homopolymer_length} == "" ]]
then 
homopolymer_length=5
fi
echo "homopolymer_length: "${homopolymer_length}

SNV_in_homopolymer_formula=$(json_extract SNV_in_homopolymer_formula "$(cat ${json})")
if [[ ${SNV_in_homopolymer_formula} == "" ]]
then 
SNV_in_homopolymer_formula="SNV_in_homopolymer==TRUE&h_ij>=0.66"
fi
echo "SNV_in_homopolymer_formula: "${SNV_in_homopolymer_formula}

INDEL_in_homopolymer_formula=$(json_extract INDEL_in_homopolymer_formula "$(cat ${json})")
if [[ ${INDEL_in_homopolymer_formula} == "" ]]
then 
INDEL_in_homopolymer_formula="INDEL_in_homopolymer==TRUE&h_ij>=0.66"
fi
echo "INDEL_in_homopolymer_formula: "${INDEL_in_homopolymer_formula}

bimera_formula=$(json_extract bimera_formula "$(cat ${json})")
if [[ ${bimera_formula} == "" ]]
then 
bimera_formula="bimera==TRUE&h_ij>=0.66"
fi
echo "bimera_formula: "${bimera_formula}

PCR_errors_formula=$(json_extract PCR_errors_formula "$(cat ${json})")
if [[ ${PCR_errors_formula} == "" ]]
then 
PCR_errors_formula="h_ij>=0.66&h_ijminor>=0.66"
fi
echo "PCR_errors_formula: "${PCR_errors_formula}

# Removing samples or loci with low performance

sample_ampl_rate=$(json_extract sample_ampl_rate "$(cat ${json})")
if [[ ${sample_ampl_rate} == "" ]]
then 
sample_ampl_rate=0.75
fi
echo "sample_ampl_rate: "${sample_ampl_rate}

locus_ampl_rate=$(json_extract locus_ampl_rate "$(cat ${json})")
if [[ ${locus_ampl_rate} == "" ]]
then 
locus_ampl_rate=0.75
fi
echo "locus_ampl_rate: "${locus_ampl_rate}

# Adding metadata

metadata=$(json_extract metadata "$(cat ${json})")
if [[ ${metadata} != "" ]]
then 
metadata=${metadata%\"}
metadata=${metadata#\"}
else
metadata="null"
fi
echo "metadata: "${metadata}

join_by=$(json_extract join_by "$(cat ${json})")
if [[ ${join_by} != "" ]]
then 
join_by=${join_by%\"}
join_by=${join_by#\"}
else
join_by="null"
fi
echo "join_by: "${join_by}

Variable1=$(json_extract Variable1 "$(cat ${json})")
if [[ ${Variable1} != "" ]]
then 
Variable1=${Variable1%\"}
Variable1=${Variable1#\"}
else
Variable1="null"
fi
echo "Variable1: "${Variable1}

Variable2=$(json_extract Variable2 "$(cat ${json})")
if [[ ${Variable2} != "" ]]
then 
Variable2=${Variable2%\"}
Variable2=${Variable2#\"}
else
Variable2="null"
fi
echo "Variable2: "${Variable2}

Longitude=$(json_extract Longitude "$(cat ${json})")
if [[ ${Longitude} != "" ]]
then 
Longitude=${Longitude%\"}
Longitude=${Longitude#\"}
else
Longitude="null"
fi
echo "Longitude: "${Longitude}

Latitude=$(json_extract Latitude "$(cat ${json})")
if [[ ${Latitude} != "" ]]
then 
Latitude=${Latitude%\"}
Latitude=${Latitude#\"}
else
Latitude="null"
fi
echo "Latitude: "${Latitude}

# Filtering desired or undesired populations or samples

var_filter=$(json_extract var_filter "$(cat ${json})")
if [[ ${var_filter} != "" ]]
then 
var_filter=${var_filter%\"}
var_filter=${var_filter#\"}
else
var_filter="null"
fi
echo "var_filter: "${var_filter}

# Reports activation

PerformanceReport=$(json_extract PerformanceReport "$(cat ${json})")
if [[ ${PerformanceReport} == "" ]]
then 
PerformanceReport="false"
fi
echo "PerformanceReport: "${PerformanceReport}

Drug_Surveillance_Report=$(json_extract Drug_Surveillance_Report "$(cat ${json})")
if [[ ${Drug_Surveillance_Report} == "" ]]
then 
Drug_Surveillance_Report="false"
fi
echo "Drug_Surveillance_Report: "${Drug_Surveillance_Report}

Variants_of_Interest_Report=$(json_extract Variants_of_Interest_Report "$(cat ${json})")
if [[ ${Variants_of_Interest_Report} == "" ]]
then 
Variants_of_Interest_Report="false"
fi
echo "Variants_of_Interest_Report: "${Variants_of_Interest_Report}

ibd_thres=$(json_extract ibd_thres "$(cat ${json})")
if [[ ${ibd_thres} == "" ]]
then 
ibd_thres="null"
fi
echo "ibd_thres: "${ibd_thres}

poly_formula=$(json_extract poly_formula "$(cat ${json})")
if [[ ${poly_formula} == "" ]]
then 
poly_formula="null"
fi
echo "poly_formula: "${poly_formula}

#Parameters for DSR

ref_gff=$(json_extract ref_gff "$(cat ${json})")
if [[ ${ref_gff} != "" ]]
then 
ref_gff=${ref_gff%\"}
ref_gff=${ref_gff#\"}
else
ref_gff="null"
fi
echo "ref_gff: "${ref_gff}

ref_fasta=$(json_extract ref_fasta "$(cat ${json})")
if [[ ${ref_fasta} != "" ]]
then 
ref_fasta=${ref_fasta%\"}
ref_fasta=${ref_fasta#\"}
else
ref_fasta="null"
fi
echo "ref_fasta: "${ref_fasta}

amplicon_fasta=$(json_extract amplicon_fasta "$(cat ${json})")
if [[ ${amplicon_fasta} != "" ]]
then 
amplicon_fasta=${amplicon_fasta%\"}
amplicon_fasta=${amplicon_fasta#\"}
else
amplicon_fasta="null"
fi
echo "amplicon_fasta: "${amplicon_fasta}

reference_alleles=$(json_extract reference_alleles "$(cat ${json})")
if [[ ${reference_alleles} != "" ]]
then 
reference_alleles=${reference_alleles%\"}
reference_alleles=${reference_alleles#\"}
else
reference_alleles="null"
fi
echo "reference_alleles: "${reference_alleles}

hap_color_palette=$(json_extract hap_color_palette "$(cat ${json})")
if [[ ${hap_color_palette} != "" ]]
then 
hap_color_palette=${hap_color_palette%\"}
hap_color_palette=${hap_color_palette#\"}
else
hap_color_palette="auto"
fi
echo "hap_color_palette: "${hap_color_palette}

gene_names=$(json_extract gene_names "$(cat ${json})")
if [[ ${gene_names} != "" ]]
then 
gene_names=${gene_names%\"}
gene_names=${gene_names#\"}
else
gene_names="null"
fi
echo "gene_names: "${gene_names}

gene_ids=$(json_extract gene_ids "$(cat ${json})")
if [[ ${gene_ids} != "" ]]
then 
gene_ids=${gene_ids%\"}
gene_ids=${gene_ids#\"}
else
gene_ids="null"
fi
echo "gene_ids: "${gene_ids}

drugs=$(json_extract drugs "$(cat ${json})")
if [[ ${drugs} != "" ]]
then 
drugs=${drugs%\"}
drugs=${drugs#\"}
else
drugs="null"
fi
echo "drugs: "${drugs}

include_all_drug_markers=$(json_extract include_all_drug_markers "$(cat ${json})")
if [[ ${include_all_drug_markers} != "" ]]
then 
include_all_drug_markers=${include_all_drug_markers%\"}
include_all_drug_markers=${include_all_drug_markers#\"}
else
include_all_drug_markers="false"
fi
echo "include_all_drug_markers: "${include_all_drug_markers}

na_var_rm=$(json_extract na_var_rm "$(cat ${json})")
if [[ ${na_var_rm} != "" ]]
then 
na_var_rm=${na_var_rm%\"}
na_var_rm=${na_var_rm#\"}
else
na_var_rm="false"
fi
echo "na_var_rm: "${na_var_rm}

na_hap_rm=$(json_extract na_hap_rm "$(cat ${json})")
if [[ ${na_hap_rm} != "" ]]
then 
na_hap_rm=${na_hap_rm%\"}
na_hap_rm=${na_hap_rm#\"}
else
na_hap_rm="false"
fi
echo "na_hap_rm: "${na_hap_rm}

# Parameters for IBD and Conectivity report

nTasks=$(json_extract nTasks "$(cat ${json})")
if [[ ${nTasks} == "" ]]
then 
nTasks=1
fi
echo "nTasks: "${nTasks}

pairwise_relatedness_table=$(json_extract pairwise_relatedness_table "$(cat ${json})")
if [[ ${pairwise_relatedness_table} != "" ]]
then 
pairwise_relatedness_table=${pairwise_relatedness_table%\"}
pairwise_relatedness_table=${pairwise_relatedness_table#\"}
else
pairwise_relatedness_table="null"
fi
echo "pairwise_relatedness_table: "${pairwise_relatedness_table}

nchunks=$(json_extract nchunks "$(cat ${json})")
if [[ ${nchunks} == "" ]]
then 
nchunks=1
fi
echo "nchunks: "${nchunks}

parallel=$(json_extract parallel "$(cat ${json})")
if [[ ${parallel} == "" ]]
then 
parallel="false"
fi
echo "parallel: "${parallel}

ibd_ncol=$(json_extract ibd_ncol "$(cat ${json})")
if [[ ${ibd_ncol} == "" ]]
then 
ibd_ncol=4
fi
echo "ibd_ncol: "${ibd_ncol}

pop_levels=$(json_extract pop_levels "$(cat ${json})")
if [[ ${pop_levels} == "" ]]
then 
pop_levels="null"
fi
echo "pop_levels: "${pop_levels}

# Parameters for COI report

poly_quantile=$(json_extract poly_quantile "$(cat ${json})")
if [[ ${poly_quantile} == "" ]]
then 
poly_quantile="null"
fi
echo "poly_quantile: "${poly_quantile}

# Run MHap_Analysis_pipeline

Rscript ${fd}/MHap_Analysis_pipeline.R \
  -wd ${wd} \
  -fd ${fd} \
  -rd ${rd} \
  -cigar_paths ${cigar_paths} \
  -cigar_files ${cigar_files} \
  -asv_table_files ${asv_table_files} \
  -asv2cigar_files ${asv2cigar_files} \
  -asv_seq_files ${asv_seq_files} \
  -zero_read_sample_list ${zero_read_sample_list} \
  -ampseqj ${ampseq_jsonfile} \
  -ampseqe ${ampseq_excelfile} \
  -ampseqc ${ampseq_csvfolder} \
  -sample_id_pattern ${sample_id_pattern} \
  -markers ${markers} \
  -o ${output} \
  -ampseq_export_format ${ampseq_export_format} \
  -min_abd ${min_abd} \
  -min_ratio ${min_ratio} \
  -off_target_formula ${off_target_formula} \
  -flanking_INDEL_formula ${flanking_INDEL_formula} \
  -homopolymer_length ${homopolymer_length} \
  -SNV_in_homopolymer_formula ${SNV_in_homopolymer_formula} \
  -INDEL_in_homopolymer_formula ${INDEL_in_homopolymer_formula} \
  -bimera_formula ${bimera_formula} \
  -PCR_errors_formula ${PCR_errors_formula} \
  -samprate ${sample_ampl_rate} \
  -lamprate ${locus_ampl_rate} \
  -metadata ${metadata} \
  -join_by ${join_by} \
  -Var1 ${Variable1} \
  -Var2 ${Variable2} \
  -Longitude ${Longitude} \
  -Latitude ${Latitude} \
  -var_filter ${var_filter} \
  -PerformanceReport ${PerformanceReport} \
  -Drug_Surveillance_Report ${Drug_Surveillance_Report} \
  -Variants_of_Interest_Report ${Variants_of_Interest_Report} \
  -ibd ${ibd_thres} \
  -poly_formula ${poly_formula} \
  -gff ${ref_gff} \
  -fasta ${ref_fasta} \
  -amp_fasta ${amplicon_fasta} \
  -reference_alleles ${reference_alleles} \
  -hap_color_palette ${hap_color_palette} \
  -gene_names ${gene_names} \
  -gene_ids ${gene_ids} \
  -drugs ${drugs} \
  -include_all_drug_markers ${include_all_drug_markers} \
  -na_var_rm ${na_var_rm} \
  -na_hap_rm ${na_hap_rm} \
  -t ${nTasks} \
  -tid 1 \
  -pairwise_relatedness_table ${pairwise_relatedness_table} \
  -nchunks ${nchunks} \
  -parallel ${parallel} \
  -ibd_ncol ${ibd_ncol} \
  -pop_levels ${pop_levels} \
  -poly_quantile ${poly_quantile}


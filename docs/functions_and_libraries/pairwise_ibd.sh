#!/bin/bash
source /broad/software/scripts/useuse
use R-4.1

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

cigar_paths=$(json_extract cigar_paths "$(cat ${json})")

if [[ ${cigar_paths} != "" ]]
then 
cigar_paths=${cigar_paths%\"}
cigar_paths=${cigar_paths#\"}
else
cigar_paths="NaN"
fi

echo "cigar_paths: "${cigar_paths}

cigar_files=$(json_extract cigar_files "$(cat ${json})")

if [[ ${cigar_files} != "" ]]
then 
cigar_files=${cigar_files%\"}
cigar_files=${cigar_files#\"}
else
cigar_files="NaN"
fi

echo "cigar_files: "${cigar_files}

ampseq_jsonfile=$(json_extract ampseq_jsonfile "$(cat ${json})")

if [[ ${ampseq_jsonfile} != "" ]]
then 
ampseq_jsonfile=${ampseq_jsonfile%\"}
ampseq_jsonfile=${ampseq_jsonfile#\"}
else
ampseq_jsonfile="NaN"
fi

echo "ampseq_jsonfile: "${ampseq_jsonfile}


ampseq_excelfile=$(json_extract ampseq_excelfile "$(cat ${json})")

if [[ ${ampseq_excelfile} != "" ]]
then 
ampseq_excelfile=${ampseq_excelfile%\"}
ampseq_excelfile=${ampseq_excelfile#\"}
else
ampseq_excelfile="NaN"
fi

echo "ampseq_excelfile: "${ampseq_excelfile}


output=$(json_extract output "$(cat ${json})")
output=${output%\"}
output=${output#\"}

echo "output: "${output}


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
markers="NaN"
fi

echo "markers: "${markers}


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


PerformanceReport=$(json_extract PerformanceReport "$(cat ${json})")

if [[ ${PerformanceReport} == "" ]]
then 
PerformanceReport="false"
fi

echo "PerformanceReport: "${PerformanceReport}


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


ref_gff=$(json_extract ref_gff "$(cat ${json})")

if [[ ${ref_gff} != "" ]]
then 
ref_gff=${ref_gff%\"}
ref_gff=${ref_gff#\"}
else
ref_gff="NaN"
fi

echo "ref_gff: "${ref_gff}

ref_fasta=$(json_extract ref_fasta "$(cat ${json})")

if [[ ${ref_fasta} != "" ]]
then 
ref_fasta=${ref_fasta%\"}
ref_fasta=${ref_fasta#\"}
else
ref_fasta="NaN"
fi

echo "ref_fasta: "${ref_fasta}


ibd_thres=$(json_extract ibd_thres "$(cat ${json})")

if [[ ${ibd_thres} == "" ]]
then 
ibd_thres="NaN"
fi

echo "ibd_thres: "${ibd_thres}


nTasks=$(json_extract nTasks "$(cat ${json})")

if [[ ${nTasks} == "" ]]
then 
nTasks=1
fi

echo "nTasks: "${nTasks}

parallel=$(json_extract parallel "$(cat ${json})")

if [[ ${parallel} == "" ]]
then 
parallel="false"
fi

echo "parallel: "${parallel}


nchunks=$(json_extract nchunks "$(cat ${json})")

if [[ ${nchunks} == "" ]]
then 
nchunks=1
fi

echo "nchunks: "${nchunks}


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


metadata=$(json_extract metadata "$(cat ${json})")

if [[ ${metadata} != "" ]]
then 
metadata=${metadata%\"}
metadata=${metadata#\"}
else
metadata="NaN"
fi

echo "metadata: "${metadata}

join_by=$(json_extract join_by "$(cat ${json})")

if [[ ${join_by} != "" ]]
then 
join_by=${join_by%\"}
join_by=${join_by#\"}
else
join_by="NaN"
fi

echo "join_by: "${join_by}

geographic_population=$(json_extract geographic_population "$(cat ${json})")

if [[ ${geographic_population} != "" ]]
then 
geographic_population=${geographic_population%\"}
geographic_population=${geographic_population#\"}
else
geographic_population="NaN"
fi

echo "geographic_population: "${geographic_population}

temporal_population=$(json_extract temporal_population "$(cat ${json})")

if [[ ${temporal_population} != "" ]]
then 
temporal_population=${temporal_population%\"}
temporal_population=${temporal_population#\"}
else
temporal_population="NaN"
fi

echo "temporal_population: "${temporal_population}

ibd_step="pairwise"

# Run pre-filtering

Rscript ${fd}/MHap_Analysis_pipeline.R \
  -wd ${wd} \
  -fd ${fd} \
  -rd ${rd} \
#  -cigar_paths ${cigar_paths} \
#  -cigar_files ${cigar_files} \
#  -ampseqj ${ampseq_jsonfile} \
  -ampseqe ${ampseq_excelfile} \
  -o ${output} \
#  -sample_id_pattern ${sample_id_pattern} \
#  -markers ${markers} \
#  -min_abd ${min_abd} \
#  -min_ratio ${min_ratio} \
#  -PerformanceReport ${PerformanceReport} \
#  -gff ${ref_gff} \
#  -fasta ${ref_fasta} \
#  -samprate ${sample_ampl_rate} \
#  -lamprate ${locus_ampl_rate} \
#  -metadata ${metadata} \
#  -join_by ${join_by} \
  -gPop ${geographic_population} \
  -tPop ${temporal_population} \
  -t ${nTasks} \
  -tid ${SGE_TASK_ID} \
  -ibd ${ibd_thres} \
  -ibd_step ${ibd_step}
  -parallel ${parallel} \
  -ibd_ncol ${ibd_ncol} \
  -pop_levels ${pop_levels}

  
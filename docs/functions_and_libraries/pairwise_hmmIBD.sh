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

loci_tsvfolder=$(json_extract loci_tsvfolder "$(cat ${json})")

if [[ ${loci_tsvfolder} != "" ]]
then 
loci_tsvfolder=${loci_tsvfolder%\"}
loci_tsvfolder=${loci_tsvfolder#\"}
else
loci_tsvfolder="NaN"
fi

echo "loci_tsvfolder: "${loci_tsvfolder}

output=$(json_extract output "$(cat ${json})")
output=${output%\"}
output=${output#\"}

echo "output: "${output}


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

pairs=$(json_extract pairs "$(cat ${json})")
if [[ ${pairs} != "" ]]
then 
pairs=${pairs%\"}
pairs=${pairs#\"}
else
pairs="null"
fi
echo "pairs: "${pairs}


nChunks=$(json_extract nChunks "$(cat ${json})")

if [[ ${nChunks} == "" ]]
then 
nChunks=1
fi

echo "nChunks: "${nChunks}



max_k=$(json_extract max_k "$(cat ${json})")

if [[ ${max_k} == "" ]]
then 
max_k=200
fi

echo "max_k: "${max_k}


by=$(json_extract by "$(cat ${json})")

if [[ ${by} == "" ]]
then 
by="null"
fi

echo "by: "${by}


max_alleles=$(json_extract max_alleles "$(cat ${json})")

if [[ ${max_alleles} == "" ]]
then 
max_alleles="null"
fi

echo "max_alleles: "${max_alleles}


ibd_step="pairwise"

# Run pre-filtering

Rscript ${fd}/run_pairwise_hmmIBD.R \
  -wd ${wd} \
  -fd ${fd} \
  -locit ${loci_tsvfolder} \
  -o ${output} \
  -t ${nTasks} \
  -tid ${SGE_TASK_ID} \
  -ibd_step ${ibd_step} \
  -parallel ${parallel} \
  -pairs ${pairs} \
  -nChunks ${nChunks} \
  -max_k ${max_k} \
  -by ${by} \
  -max_alleles ${max_alleles}

  
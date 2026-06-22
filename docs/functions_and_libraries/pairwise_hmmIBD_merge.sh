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


output=$(json_extract output "$(cat ${json})")
output=${output%\"}
output=${output#\"}

echo "output: "${output}

ibd_step="merge"

# Run pre-filtering

Rscript ${fd}/run_pairwise_hmmIBD.R \
  -wd ${wd} \
  -fd ${fd} \
  -o ${output} \
  -ibd_step ${ibd_step}

  
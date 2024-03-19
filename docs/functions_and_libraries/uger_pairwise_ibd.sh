#!/bin/bash
source /broad/software/scripts/useuse
use UGER

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

vmem=$(json_extract vmem "$(cat ${json})")
cores=$(json_extract cores "$(cat ${json})")
h_rt1=$(json_extract h_rt1 "$(cat ${json})")
h_rt2=$(json_extract h_rt2 "$(cat ${json})")
h_rt3=$(json_extract h_rt3 "$(cat ${json})")
nTasks=$(json_extract nTasks "$(cat ${json})")


wd=$(json_extract wd "$(cat ${json})")

wd=${wd%\"}
wd=${wd#\"}


fd=$(json_extract fd "$(cat ${json})")
fd=${fd%\"}
fd=${fd#\"}

cd ${wd}

# Submit Task Array

if [${nTasks} == 1]; then
  qsub -l h_vmem=${vmem}G \
    -l h_rt=${h_rt2} \
    -j y \
    -o ${wd}/pairwise_ibd_output \
    ${fd}/pairwise_ibd.sh ${json}
else
  qsub -l h_vmem=${vmem}G \
    -l h_rt=${h_rt2} \
    -j y \
    -o ${wd}/pairwise_ibd_output \
    -t 1-${nTasks} \
    ${fd}/pairwise_ibd.sh ${json}
fi



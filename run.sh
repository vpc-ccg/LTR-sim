#!/bin/bash
set -e
SCRIPT_PATH="$(dirname `which $0`)"
SNAKEFILE=Snakefile
#Checking if tools exist
command -v snakemake && echo "... Exists :)" >&2 || { echo snakemake missing! >&2; FAILED=1; }

if [ -z $FAILED ]; then 
    if [[ "$@" == *"cluster-config"* ]]; then
        mkdir -p ~/.slurm-logs
        #
        snakemake  -s "${SCRIPT_PATH}/${SNAKEFILE}" --config SCPA=${SCRIPT_PATH}  "$@" --cluster "sbatch --parsable -J {cluster.name} -p \"long\" -c {cluster.c} --mem {cluster.mem}  -t \"23:30:00\" --output '$HOME/.slurm-logs/%j.out' --error '$HOME/.slurm-logs/%j.err'" --cluster-status $SCRIPT_PATH/slurm_status.py;
    else
        snakemake -s "${SCRIPT_PATH}/${SNAKEFILE}" --config SCPA=${SCRIPT_PATH}  "$@";
    fi
fi

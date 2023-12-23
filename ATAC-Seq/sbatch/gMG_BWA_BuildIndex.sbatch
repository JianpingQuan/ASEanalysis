#!/bin/bash

# Job name:
#SBATCH --job-name=Build_Index
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=2
#
# Memory per node:
#SBATCH --mem=300G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=8:00:00

# ============================================================

echo $(date) start.process

source "$env"

cd "$out"/$sample

$BWA index -p "$sample" "$dir"/"$sample"/"$fasta" > "$out"/$sample/"$sample".bwa.buildIndex.log 2>&1

echo $(date) done.main.process

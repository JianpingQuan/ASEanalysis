#!/bin/bash

# Job name:
#SBATCH --job-name=merge-vcfs
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
#SBATCH --mem=64G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=24:00:00

# ============================================================

source "$env"

java -Xmx"$mem" -XX:ParallelGCThreads=2 -jar $PICARD MergeVcfs \
  I=$input \
  O=$output \
  D=$DICT \
  tmp_dir=$tmp > $out/log/$prefix.mergeVCFs.log 2>&1

echo $(date) done.main.process

#!/bin/bash

# Job name:
#SBATCH --job-name=hardFilterSNP
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
#SBATCH --time=4:00:00

# ============================================================

source "$env"

java -Xmx"$mem" -XX:ParallelGCThreads=2 -Djava.io.tmpdir=$tmp \
  -jar $GATK VariantFiltration \
  -V $vcf \
  -O $output \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 50.0" --filter-name "QUAL50" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "FS > 100.0" --filter-name "FS100" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" > $out/log/snp.hardFilter.log 2>&1	

echo $(date) done.main.process

#!/bin/bash

# Job name:
#SBATCH --job-name=split-fastq-unknow
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
#SBATCH --mem=8G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=16:00:00

# ============================================================

echo "SLURM_JOB_ID: $SLURM_JOB_ID"

# check if there is unpaired lines
# ============================================================

unpline=`paste <(gunzip -c "$file1" | paste - - - - | cut -f 1 | sed 's/\/.*//' | awk '{print $1}') <(gunzip -c "$file2" | paste - - - - | cut -f 1 | sed 's/\/.*//' | awk '{print $1}') | awk '$1 != $2' | head | wc -l`

if [[ $unpline -gt 0 ]]
then
	echo $(date +"%m-%d-%Y-%T"):error: non-matching reads
	exit
fi
	
# get instrument info
# ============================================================

gunzip -c "$file1" | paste - - - - | cut -f 1 | sed 's/^@//' | awk -F ":" '{print "pe."$1"."$2"."$3"."$4" "$1" "$2" "$3" "$4}' | uniq | sort | uniq > "$dir"/file.info

# split fastq
# ============================================================

while read line           
do           
  # output fastq file
	ins=$(echo $line | cut -d " " -f 2)
	run=$(echo $line | cut -d " " -f 3)
	flc=$(echo $line | cut -d " " -f 4)
	lane=$(echo $line | cut -d " " -f 5)
	gunzip -c $file1 | grep --no-group-separator -A 3 ^"@$ins:$run:$flc:$lane" | gzip > $dir/"pe.$ins.$run.$flc.$lane"_1.fastq.gz
	gunzip -c $file2 | grep --no-group-separator -A 3 ^"@$ins:$run:$flc:$lane" | gzip > $dir/"pe.$ins.$run.$flc.$lane"_2.fastq.gz
done < "$dir"/file.info

chmod 440 $dir/*.fastq.gz

echo $(date +"%m-%d-%Y-%T"):info: done.main.process

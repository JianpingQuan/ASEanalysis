#! /bin/bash
# ============================================================

# pipeline script for processing fastq to gvcf
# ============================================================

# get command line options
# ============================================================

args=`getopt -o "r:,d:,t:,o:,s:" -l "resource:,dir:,tmp:,out:,sample:" -- "$@"`
>&2 echo "command arguments given: $args"
eval set -- "$args"

tmp=/tmp/

# parse arguments
# ============================================================

while true;
do
  case $1 in

    -r|--resource)
      resource=$2
      shift 2;;

    -d|--dir)
			dir=$2
			shift 2;; 

  	-t|--tmp)
			tmp=$2
			shift 2;;

    -o|--out)
			out=$2
			shift 2;;

    -s|--sample)
		  sample=$2
			shift 2;;

    --)
      shift
      break;;

  esac
done

# prepare directory
# ============================================================

id=$(date | md5sum | cut -d " " -f 1 | awk '{print "'$sample'"$0}')
if [[ -e $tmp/$id ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: $tmp/$id exists."
  exit 1
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: makedir $tmp/$id."
  mkdir $tmp/$id
	mkdir $tmp/$id/tmp
fi

if [[ -e $out/$sample ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: $out exists."
  exit 1
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: makedir $out/$sample."
  mkdir $out/$sample
	mkdir $out/$sample/log
fi

# check for file availability
# ============================================================

if [[ -e $resource/resource.env && -e $resource/sbatch/checkFastqSex.sbatch && -e $resource/sbatch/bwa.sbatch && -e $resource/sbatch/merge.sbatch && -e $resource/sbatch/markDuplicates.sbatch && -e $resource/sbatch/bqsr.sbatch && -e $resource/sbatch/applyBQSR.sbatch && -e $resource/sbatch/haplotypeCaller.sbatch  && -e $resource/sbatch/cov.sbatch ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: all slurm scripts found."
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: cannot find some slurm scripts."
  exit 1
fi

# check for directory structure
# ============================================================

if [[ -e $dir/$sample/file.info ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: found file.info in $dir/$sample."
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: cannot find file.info in $dir/$sample."
  exit 1
fi

# ==============================
# = 1. check sex of fastq file =
# ==============================

filepre=`head -n 1 $dir/$sample/file.info | cut -d " " -f 1`

# run sbatch
jobID=$(sbatch --export=env=$resource/resource.env,file=$dir/$sample/"$filepre"_1.fastq.gz,head=25000,tmp=$tmp/$id/,outfile=$out/$sample/sex.info --output=$out/$sample/log/checkFastqSex.out --error=$out/$sample/log/checkFastqSex.err $resource/sbatch/checkFastqSex.sbatch | cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success
reRunCounter=0

until [ `wc -l $out/$sample/log/checkFastqSex.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $out/$sample/log/checkFastqSex.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 1 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting sex check job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/checkFastqSex.fail.out
	sexTime=$(($reRunCounter))
	jobID=$(sbatch --export=env=$resource/resource.env,file=$dir/$sample/"$filepre"_1.fastq.gz,head=100000,outfile=$out/$sample/sex.info --time=$sexTime:30:00 --output=$out/$sample/log/checkFastqSex.out --error=$out/$sample/log/checkFastqSex.err $resource/sbatch/checkFastqSex.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 5m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $out/$sample/log/checkFastqSex.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" $out/$sample/log/checkFastqSex.out | wc -l | awk '{print $1}'` -gt 0 ]
then
	sex=`awk '{print $1}' $out/$sample/sex.info`
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed checkFastqSex. inferred sex = $sex."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/checkFastqSex.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for checkFastqSex."
	exit 1
fi

# ====================
# = 2. bwa alignment =
# ====================

# run sbatch
jobID=$(sbatch --export=env=$resource/resource.env,dir=$dir/$sample,tmp=$tmp/$id/,out=$out/$sample,sample=$sample --output=$out/$sample/log/bwa.out --error=$out/$sample/log/bwa.err $resource/sbatch/bwa.sbatch | cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success
reRunCounter=0

until [ `wc -l $out/$sample/log/bwa.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $out/$sample/log/bwa.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 2 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting bwa job."
	rm $tmp/$id/*.bam*
	rm $out/$sample/seq.sum
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/bwa.fail.out
	bwaTime=$(($reRunCounter * 12 + 60))
	bwaMem=$(($reRunCounter * 4 + 28))G
	jobID=$(sbatch --export=env=$resource/resource.env,dir=$dir/$sample,tmp=$tmp/$id/,out=$out/$sample,sample=$sample --time=$bwaTime:00:00 --mem=$bwaMem --output=$out/$sample/log/bwa.out --error=$out/$sample/log/bwa.err $resource/sbatch/bwa.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 5m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $out/$sample/log/bwa.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" $out/$sample/log/bwa.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed bwa."
  sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/bwa.out
else
	>&2 echo "$(date):error: error for bwa."
	exit 1
fi

# merge alignments
# run sbatch
jobID=$(sbatch --export=env=$resource/resource.env,dir=$dir/$sample,tmp=$tmp/$id,out=$out/$sample --output=$out/$sample/log/merge.out --error=$out/$sample/log/merge.err $resource/sbatch/merge.sbatch | cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success
reRunCounter=0

until [ `wc -l $out/$sample/log/merge.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $out/$sample/log/merge.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 2 ]
do
	echo "place" > $tmp/$id/sort.bam.txt
	rm $tmp/$id/sort.bam*
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/merge.fail.out
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting merge job."
	mergeTime=$(($reRunCounter * 4 + 16))
	mergeMem=$(($reRunCounter * 4 + 36))G
	jobID=$(sbatch --export=env=$resource/resource.env,dir=$dir/$sample,tmp=$tmp/$id,out=$out/$sample --time=$mergeTime:00:00 --mem=$mergeMem --output=$out/$sample/log/merge.out --error=$out/$sample/log/merge.err $resource/sbatch/merge.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 5m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $out/$sample/log/merge.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" $out/$sample/log/merge.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed merge alignments."
	rm $tmp/$id/*.map.bam
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/merge.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for merge."
	exit 1
fi

# ======================
# = 3. mark duplicates =
# ======================

jobID=$(sbatch --export=env=$resource/resource.env,tmp=$tmp/$id/tmp,input=$tmp/$id/sort.bam,output=$tmp/$id/rmdup.bam,out=$out/$sample,mem=32G,metric=$out/$sample/rmdup.metric.out --output=$out/$sample/log/markDuplicates.out --error=$out/$sample/log/markDuplicates.err $resource/sbatch/markDuplicates.sbatch | cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success
reRunCounter=0

until [ `wc -l $out/$sample/log/markDuplicates.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $out/$sample/log/markDuplicates.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 2 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting rmdup job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/markDuplicates.fail.out
	echo "place" > $tmp/$id/rmdup.txt
	rm $tmp/$id/rmdup.*
	rmdupTime=$(($reRunCounter * 4 + 28))
	rmdupMem=$(($reRunCounter * 32 + 96))G
	jobID=$(sbatch --export=env=$resource/resource.env,tmp=$tmp/$id/tmp,input=$tmp/$id/sort.bam,output=$tmp/$id/rmdup.bam,out=$out/$sample,mem=$rmdupMem,metric=$out/$sample/rmdup.metric.out --mem=$rmdupMem --time=$rmdupTime:00:00 --output=$out/$sample/log/markDuplicates.out --error=$out/$sample/log/markDuplicates.err $resource/sbatch/markDuplicates.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 5m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $out/$sample/log/markDuplicates.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" $out/$sample/log/markDuplicates.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed mark duplicates."
	rm $tmp/$id/sort.*
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/markDuplicates.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for mark duplicates."
	exit 1
fi

# =====================
# = 4. count coverage =
# =====================

jobID=$(sbatch --export=env=$resource/resource.env,input=$tmp/$id/rmdup.bam,output=$out/$sample/genome.cov.out,out=$out/$sample --output=$out/$sample/log/cov.out --error=$out/$sample/log/cov.err $resource/sbatch/cov.sbatch | cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success
reRunCounter=0

until [ `wc -l $out/$sample/log/cov.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $out/$sample/log/cov.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt 2 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting cov job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/cov.fail.out
	covTime=$(($reRunCounter * 4 + 16))
	covMem=$(($reRunCounter * 4 + 8))G
	jobID=$(sbatch --export=env=$resource/resource.env,input=$tmp/$id/rmdup.bam,output=$out/$sample/genome.cov.out,out=$out/$sample --time=$covTime:00:00 --mem=$covMem --output=$out/$sample/log/cov.out --error=$out/$sample/log/cov.err $resource/sbatch/cov.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 5m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $out/$sample/log/cov.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $out/$sample/log/cov.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed coverage count."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/cov.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for coverage count."
	exit 1
fi

# =================================
# = 5. base quality recalibration =
# =================================

jobID=$(sbatch --export=env=$resource/resource.env,tmp=$tmp/$id/tmp,input=$tmp/$id/rmdup.bam,output=$tmp/$id/recal.table,out=$out/$sample,mem=16G --output=$out/$sample/log/bqsr.out --error=$out/$sample/log/bqsr.err $resource/sbatch/bqsr.sbatch | cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success
reRunCounter=0

until [ `wc -l $out/$sample/log/bqsr.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $out/$sample/log/bqsr.out | wc -l | awk '{print $1}'` -gt 0 ] && [ `grep "ERROR" $out/$sample/log/bqsr.log | wc -l | awk '{print $1}'` -eq 0 ] || [ $reRunCounter -gt 2 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting bqsr job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/bqsr.fail.out
	bqsrTime=$(($reRunCounter * 4 + 16))
	bqsrMem=$(($reRunCounter * 8 + 16))G
	jobID=$(sbatch --export=env=$resource/resource.env,tmp=$tmp/$id/tmp,input=$tmp/$id/rmdup.bam,output=$tmp/$id/recal.table,out=$out/$sample,mem=$bqsrMem --mem=$bqsrMem --time=$bqsrTime:00:00 --output=$out/$sample/log/bqsr.out --error=$out/$sample/log/bqsr.err $resource/sbatch/bqsr.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 5m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $out/$sample/log/bqsr.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" $out/$sample/log/bqsr.out | wc -l | awk '{print $1}'` -gt 0 ] && [ `grep "ERROR" $out/$sample/log/bqsr.log | wc -l | awk '{print $1}'` -eq 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed bqsr."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/bqsr.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for bqsr."
	exit 1
fi

# apply bqsr
jobID=$(sbatch --export=env=$resource/resource.env,tmp=$tmp/$id/tmp,input=$tmp/$id/rmdup.bam,output=$tmp/$id/recal.bam,recal=$tmp/$id/recal.table,out=$out/$sample/,mem=16G --output=$out/$sample/log/applyBQSR.out --error=$out/$sample/log/applyBQSR.err $resource/sbatch/applyBQSR.sbatch | cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success
reRunCounter=0

until [ `wc -l $out/$sample/log/applyBQSR.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $out/$sample/log/applyBQSR.out | wc -l | awk '{print $1}'` -gt 0 ] && [ `grep "ERROR" $out/$sample/log/applyBQSR.log | wc -l | awk '{print $1}'` -eq 0 ] || [ $reRunCounter -gt 2 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting apply-bqsr job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/applyBQSR.fail.out
	bqsrTime=$(($reRunCounter * 4 + 16))
	jobID=$(sbatch --export=env=$resource/resource.env,tmp=$tmp/$id/tmp,input=$tmp/$id/rmdup.bam,output=$tmp/$id/recal.bam,recal=$tmp/$id/recal.table,out=$out/$sample/,mem=16G --time=$bqsrTime:00:00 --output=$out/$sample/log/applyBQSR.out --error=$out/$sample/log/applyBQSR.err $resource/sbatch/applyBQSR.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 5m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $out/$sample/log/applyBQSR.err | awk '{print $1}'` == 0 ] && [ `grep "done.main.process" $out/$sample/log/applyBQSR.out | wc -l | awk '{print $1}'` -gt 0 ] && [ `grep "ERROR" $out/$sample/log/applyBQSR.log | wc -l | awk '{print $1}'` -eq 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed applyBQSR."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/applyBQSR.out
	rm $tmp/$id/rmdup.* $tmp/$id/recal.table
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for applyBQSR."
	exit 1
fi

# =======================
# = 6. haplotype caller =
# =======================

# do a full genome ploidy = 2 and then only X chromosome ploidy = 1 for all samples
# regardless of sex
# because 1) sex call could change
# 2) PAR can be determined later
# 3) can use ploidy = 2 to test for heterozygosity

jobID=$(sbatch --export=env=$resource/resource.env,tmp=$tmp/$id/tmp,input=$tmp/$id/recal.bam,output1=$out/$sample/hc1.p2.g.vcf.gz,output2=$out/$sample/hc2.p2.g.vcf.gz,output3=$out/$sample/hc3.p2.g.vcf.gz,output4=$out/$sample/hc4.p2.g.vcf.gz,output5=$out/$sample/hc5.p1.g.vcf.gz,out=$out/$sample,mem=12G --output=$out/$sample/log/haplotypeCaller.out --error=$out/$sample/log/haplotypeCaller.err $resource/sbatch/haplotypeCaller.sbatch | cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

# check if run success
reRunCounter=0

until [ `wc -l $out/$sample/log/haplotypeCaller.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $out/$sample/log/haplotypeCaller.out | wc -l | awk '{print $1}'` -gt 0 ] && [ `cat $out/$sample/log/hc*.log | grep "ERROR" | wc -l | awk '{print $1}'` -eq 0 ] || [ $reRunCounter -gt 2 ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting haplotypeCaller job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/haplotypeCaller.fail.out
	rm $out/$sample/*.g.vcf.gz
	hcTime=$(($reRunCounter * 16 + 64))
	hcMem=$(($reRunCounter * 4 + 16))G
	hcMem2=$(($reRunCounter * 20 + 80))G
	jobID=$(sbatch --export=env=$resource/resource.env,tmp=$tmp/$id/tmp,input=$tmp/$id/recal.bam,output1=$out/$sample/hc1.p2.g.vcf.gz,output2=$out/$sample/hc2.p2.g.vcf.gz,output3=$out/$sample/hc3.p2.g.vcf.gz,output4=$out/$sample/hc4.p2.g.vcf.gz,output5=$out/$sample/hc5.p1.g.vcf.gz,out=$out/$sample,mem=$hcMem --mem=$hcMem2 --time=$hcTime:00:00 --output=$out/$sample/log/haplotypeCaller.out --error=$out/$sample/log/haplotypeCaller.err $resource/sbatch/haplotypeCaller.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 5m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l $out/$sample/log/haplotypeCaller.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" $out/$sample/log/haplotypeCaller.out | wc -l | awk '{print $1}'` -gt 0 ] && [ `cat $out/$sample/log/hc*.log | grep "ERROR" | wc -l | awk '{print $1}'` -eq 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed haplotype caller."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> $out/$sample/log/haplotypeCaller.out
	rm -r $tmp/$id
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):info: error for haplotype caller."
	exit 1
fi

#! /bin/bash
# ============================================================

# bash to call genotypes from combined gvcfs
# ============================================================

# get command line options
# ============================================================

args=`getopt -o "r:,d:,i:,t:,p:,u:" -l "resource:,dir:,input:,tmp:,pre:,rerun:" -- "$@"`
>&2 echo "command arguments given: $args"
eval set -- "$args"

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
    
    -i|--input)
      input=$2
      shift 2;;

    -t|--tmp)
      tmp=$2
      shift 2;;
      
    -p|--pre)
      pre=$2
      shift 2;;
    
    -u|--rerun)
      rerun=$2
      shift 2;;

    --)
      shift
      break;;
      
  esac
done

# check for file availability
# ============================================================

if [[ -e $resource/resource.env && -e $resource/sbatch/genotypeGVCFs.sbatch ]]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: all slurm scripts found."
else
  >&2 echo "$(date +"%m-%d-%Y-%T"):error: cannot find some slurm scripts."
  exit 1
fi

# run sbatch
# ============================================================

jobID=$(sbatch --export=env=/mnt/research/qgg/resource/swim/resource.env,dir=$dir,tmp=$tmp,prefix=$pre,input=$input,mem=64G --output="$dir"/log/genotypeGVCFs.$pre.out --error="$dir"/log/genotypeGVCFs.$pre.err $resource/sbatch/genotypeGVCFs.sbatch | cut -d " " -f 4)

sleep 1m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
  sleep 5m
  jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

reRunCounter=0
until [ `wc -l "$dir"/log/genotypeGVCFs.$pre.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" "$dir"/log/genotypeGVCFs.$pre.out | wc -l | awk '{print $1}'` -gt 0 ] || [ $reRunCounter -gt $rerun ]
do
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: resubmitting genotypeGVCFs job."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> "$dir"/log/genotypeGVCFs.$pre.fail.out
	runTime=$(($reRunCounter * 24 + 72))
	runMem=$(($reRunCounter * 16 + 80))G
  jobID=$(sbatch --export=env=/mnt/research/qgg/resource/swim/resource.env,dir=$dir,tmp=$tmp,prefix=$pre,input=$input,mem=$runMem --time=$runTime:00:00 --mem=$runMem --output="$dir"/log/genotypeGVCFs.$pre.out --error="$dir"/log/genotypeGVCFs.$pre.err $resource/sbatch/genotypeGVCFs.sbatch | cut -d " " -f 4)
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
	do
		sleep 5m
		jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
	done
	reRunCounter=$(($reRunCounter + 1))
done

if [ `wc -l "$dir"/log/genotypeGVCFs.$pre.err | awk '{print $1}'` -eq 0 ] && [ `grep "done.main.process" "$dir"/log/genotypeGVCFs.$pre.out | wc -l | awk '{print $1}'` -gt 0 ]
then
  >&2 echo "$(date +"%m-%d-%Y-%T"):info: completed genotypeGVCFs."
	sacct --units=G --format="User,state%15,Account%15,JobID%20,JobName%15,TimeLimit,Elapsed,TotalCPU%12,NCPUS,NTasks,NodeList%25,NNodes,ReqMem,AveRSS,MaxRSS,Submit,Start,End,Partition%48,WorkDir%50,ReqTRES,AllocTRES" -j $jobID >> "$dir"/log/combineGVCFs.$pre.out
else
	>&2 echo "$(date +"%m-%d-%Y-%T"):error: error for genotypeGVCFs."
	exit 1
fi

#!/bin/sh

#SBATCH --account=m22oc-s2448741
#SBATCH --partition=standard
#SBATCH --exclusive

#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --qos=short
#SBATCH --reservation=shortqos

#SBATCH --output=profile.out
#SBATCH --error=profile.err
#SBATCH --cpus-per-task=1


run_name=$1
echo $run_name

datetime=$(date +'%m-%d-%T')

reactor_results="reactorState_"$run_name"_"$datetime".out"
reactor_time="time_"$run_name"_"$datetime".out"
profiling="profiling_"$run_name"_"$datetime".out"

echo $reactor_results

make

./reactor "config_"$run_name".txt" $reactor_results > $reactor_time

gprof reactor gmon.out > $profiling

mv $reactor_results $reactor_time $profiling results

echo "Done"

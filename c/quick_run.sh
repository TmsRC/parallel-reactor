#!/bin/sh

#run_name=lean
run_name=$1
num_procs=$2
echo $run_name

datetime=$(date +'%m-%d-%T')

reactor_results="reactorState_"$run_name"_"$datetime".out"
reactor_time="time_"$run_name"_"$datetime".out"
#profiling="profiling_"$run_name"_"$datetime".out"

echo $reactor_results

make clean
make

mpirun -n $num_procs ./reactor "config_"$run_name".txt" $reactor_results > $reactor_time

#gprof reactor gmon.out > $profiling

#mv $reactor_results $reactor_time $profiling results
mv $reactor_results $reactor_time results

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

cd $SLURM_SUBMIT_DIR


module load -s intel-compilers-19
module load -s intel-vtune-19
module load -s intel-20.4/vtune
module load intel-mpi-19

run_name=$1
num_procs=$2
vtune_file=vtune-hotspots$3
echo $run_name

datetime=$(date +'%m-%d-%T')
hostname=.cirrus-login1

reactor_results="reactorState_"$run_name"_"$datetime".out"
reactor_time="time_"$run_name"_"$datetime".out"
profiling="profiling_"$run_name"_"$datetime".out"

echo $reactor_results

make clean
make

mpirun -n $num_procs amplxe-cl -collect hotspots -r $vtune_file ./reactor "config_"$run_name".txt" $reactor_results > $reactor_time

vtune -report gprof-cc -r $vtune_file$hostname > $profiling

mv $reactor_results $reactor_time $profiling results



echo "Done"

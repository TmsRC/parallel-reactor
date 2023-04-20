#!/bin/sh

#SBATCH --account=m22oc-s2448741
#SBATCH --partition=standard
#SBATCH --exclusive

#SBATCH --time=00:05:00
#SBATCH --nodes=4
#SBATCH --tasks-per-node=36
#SBATCH --qos=standard
# SBATCH --reservation=shortqos

#SBATCH --output=scaling.out
#SBATCH --error=scaling.err
#SBATCH --cpus-per-task=1

cd $SLURM_SUBMIT_DIR


module load -s intel-compilers-19
module load mpt

make clean
make

> proof_of_scale.out

srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=4 --ntasks=144 ./reactor config_scale.txt reactorState_scale.out >> proof_of_scale.out

echo "Done"

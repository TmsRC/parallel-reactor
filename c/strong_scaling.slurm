#!/bin/sh

#SBATCH --account=m22oc-xxxxxx
#SBATCH --partition=standard
#SBATCH --exclusive

#SBATCH --time=02:00:00
# SBATCH --nodes=2
# SBATCH --tasks-per-node=36
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

> strong_scaling.out

NPROCS=( 2 4 8 16 24 36 48 64 72)
NODES=( 1 1 1 1 1 1 2 2 2)

#NPROCS=( 128 144 )
#NODES=( 4 4 )

for i in {0..8}
do
    for j in {1..3}
    do
      srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=${NODES[i]} --ntasks=${NPROCS[i]} ./reactor config_lean.txt reactorState.out >> strong_scaling.out
    done
done

echo "Done"

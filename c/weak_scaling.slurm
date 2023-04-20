#!/bin/sh

#SBATCH --account=m22oc-xxxxxx
#SBATCH --partition=standard
#SBATCH --exclusive

#SBATCH --time=02:00:00
#SBATCH --nodes=1
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

> weak_scaling.out

NPROCS=( 2 4 8 16 24 36 48 64 72)
NODES=( 1 1 1 1 1 1 )
#PROBSIZE=( 100 200 400 800 1200 1800 2400 3200 3600 )
PROBSIZE=( 200 400 800 1600 2400 3600 )

#NPROCS=( 128 144 )
#NODES=( 4 4 )

#for i in {0..8}
for i in {0..5}
do
    for j in {1..3}
    do
      srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=${NODES[i]} --ntasks=${NPROCS[i]} ./reactor "config_"${PROBSIZE[i]}".txt" reactorState.out >> weak_scaling.out
    done
done

echo "Done"

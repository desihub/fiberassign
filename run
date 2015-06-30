#PBS -S /bin/bash
#PBS -N Assign
#PBS -l mppwidth=24,walltime=1:29:30
#PBS -q premium
#PBS -V

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=24
echo "# Running assign"

aprun -n 1 -N 1 -d 24 ./assign features.txt

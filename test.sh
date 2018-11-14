#!/bin/bash
#PBS -l walltime=00:15:00
#PBS -l nodes=1:ppn=20
#PBS -W group_list=newriver
#PBS -q normal_q
#PBS -A CMDA3634

cd $PBS_0_WORKDIR

module purge; module load gcc lm

gcc -O3 -o mandelbrot cudaMandelbrot.cu -lm

for i in `seq 1 20`
do
    ./cudaMandelbrot 4096 4096 $i
end

exit;

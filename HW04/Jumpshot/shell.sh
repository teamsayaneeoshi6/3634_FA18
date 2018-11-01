#!/bin/bash
#PBS -l walltime=00:15:00
#PBS -l nodes=1:ppn=30
#PBS -W group_list=newriver
#PBS -q normal_q
#PBS -A CMDA3634

cd 3634_FA18/HW04/Scaling

module purge; module load gcc mvapich2 jdk mpe2

export LDFLAG="-L$MPE_LIB -llmpe -lmpe -lm -mpilog -lpthread -lstdc++"

mpecc -o simpleRayTracer simpleRayTracer.h $LDFLAGS

for i in `seq 1 30`
do

    mpiexec -n $i ./simpleRayTracer 10

end

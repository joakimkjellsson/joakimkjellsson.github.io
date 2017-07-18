#!/bin/bash

module load easy
module load mpi
module load intelmpi
module load netcdf
time mpirun -np $SP_PROCS -machinefile $SP_HOSTFILE ./shallowwater.x 

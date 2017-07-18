#!/bin/bash

# Load the appropriate module of the batch system, t.ex.
module load easy

#Start the mpd daemons: one on every compute node
NR_NODES=$(cat ${SP_HOSTFILE} | wc -l)
echo $NR_NODES
mpdboot -n ${NR_NODES} -f ${SP_HOSTFILE}

# Now start the MPI application
mpiexec -n $1 ./shallowwater.x
# Maybe use another MPI program
#mpiexec -n <desired-number-of-mpi-processes> ./my_mpi_prog_2

# Done with MPI. Shutdown the mpd daemons.
mpdallexit
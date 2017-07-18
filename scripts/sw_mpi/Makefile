NC      = /pdc/vol/netcdf/4.0.1
MPI     = /software/mpi/openmpi/1.4.3/i1203
HDF     = /pdc/vol/hdf5/1.8.3
SZIP    = /pdc/vol/szip/2.1
ZLIB    = /pdc/vol/zlib/1.2.3
F95     = mpiifort
OPT     = -O2 -g -traceback -CB
LIB     = -L$(NC)/lib/ -L$(HDF)/lib/ -L$(SZIP)/lib/ -L$(ZLIB)/lib/
INC     = -I$(NC)/include/ -I$(HDF)/include/
MODS    = -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz

EXE     = shallowwater.x

default: shallowwater.x

clean:
	rm -f $(EXE)

shallowwater.x:
	$(F95) -o $(EXE) $(LIB) $(INC) $(OPT) src/shallowwater.F90 $(MODS)


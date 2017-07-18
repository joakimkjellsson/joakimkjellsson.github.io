##
##
## Post processing script for the shallow water model 
## used in the HPC1 course
##
## The model runs on MPI and each process writes its own
## netCDF file. This script patches the subdomains together. 
##
##
## Author: Joakim Kjellsson, October 2012
##
##

import Nio, os, sys
import numpy as np


##
## Prefix and suffix (none at the moment) for the files
##
prefix = '/nobackup/vagn2/x_joakj/data/hpc_project/shallow_10.nc.'
suffix = ''

print ' ---------------------------------------'
print ' Post process '
print ' ---------------------------------------'

##
## Create a list of all files
##
files = []
j = 1
file = (prefix+'%04d' % j)
while (os.path.isfile(file)):
   files.append(file)
   print ' Added file ',file
   j = j + 1
   file = (prefix+'%04d' % j)

##
## Create the master file using dimensions and variables 
## from the first file in the list
##
if (os.path.isfile(prefix+'master.nc')):
   os.system('rm '+prefix+'master.nc')
master = Nio.open_file(prefix+'master',mode='c',format='nc')

f = Nio.open_file(files[0],'r',format='nc')
for dim in f.dimensions:
   dimlength = f.dimensions[dim]
   if (dim == 'x'):
      dimlength = dimlength * len(files)
   master.create_dimension(dim,dimlength)

for var in f.variables:
   v = f.variables[var]
   kind = v.typecode()
   numDims = v.rank
   dimSizes = v.shape
   dimNames = v.dimensions
   master.create_variable(var,kind,dimNames)

   
##
## For each file and each variable, put the subdomain
## into the master file
##

for file in files:
   f = Nio.open_file(file,'r',format='nc')
   
   ji1 = f.variables['ji1'].get_value() # west point
   ji2 = f.variables['ji2'].get_value() # east point
   jj1 = f.variables['jj1'].get_value() # south point
   jj2 = f.variables['jj2'].get_value() # north point
   
   for var in f.variables:
      v = f.variables[var]
      numDims = v.rank
      print ' Writing ',var,' between ',ji1-1,' and ',ji2
      if (numDims == 1):
         if (var == 'vx'):
            master.variables[var][ji1-1:ji2] = f.variables[var][:]
         else:
            master.variables[var][:] = f.variables[var][:]
      elif (numDims == 2):
         master.variables[var][:,ji1-1:ji2] = f.variables[var][:,:]
      elif (numDims == 3):
         master.variables[var][:,:,ji1-1:ji2] = f.variables[var][:,:,:]
         

print ' ------------------------------------'
print ' The post process script successfully saved output to the master file '
print ' ------------------------------------'   
print master

master.close()

print ' ----------- END --------------------'

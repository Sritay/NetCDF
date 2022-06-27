# NetCDF utilities for LAMMPS
This folder contains utilities dealing with the NetCDF file format.

## Convert plaintext dump files to NetCDF files
File "nc_cdl_builder.cpp" consists of a C++ script which reads in
a plaintext LAMMPS dump file and turns it into cdlpart files,
which can then be used to build a NetCDF file using standard
NetCDF utilities, such as ncgen.

To compile: 
```
g++ nc_cdl_builder.cpp
```
To run:
```
./a.out
```
## Read LAMMPS NetCDF dump in parallel
File "pn_readdump.c" consists of a C script which reads a LAMMPS
NetCDF dump file in parallel using pNetCDF and MPI. This is useful 
where the filesystem supports parallel IO, such as lustre. The 
dump file must also be striped to multiple OSTs.
Parallelization is achieved by assigning different timesteps to 
different procs.

To compile: 
```
cc (or any C compiler with MPI wrapper)  pn_readdump.c -o pn_readdump  -lpnetcdf
```
(ensure pnetcdf library is loaded)

To run: 
```
mpirun/aprun/srun [MPI args] pn_readdump dump_file_name.nc
```

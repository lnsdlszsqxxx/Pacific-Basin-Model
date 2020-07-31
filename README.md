# How to install NCAR CSM Pacific Basin Model (Gent and Cane model) on my own MacBook Pro

## My laptop version
- MacBook Pro 10.15.5

## Download the Model
- Download source code from [here](http://www.cesm.ucar.edu/models/ocn-pac/)
- the user's manual is [here](http://www.cesm.ucar.edu/models/ocn-pac/pac3.e/doc0.html)

## Pre-installation
- a fortran compiler (I use gfortran)
- download a source file named [dgtsl.F](http://www.netlib.no/netlib/linpack/dgtsl.f), which solves a tridiagonal linear system and is invoked by integrate.F. Put it together with other source code files.
- [netCDF-C library](https://www.unidata.ucar.edu/downloads/netcdf/)
- [netCDF-Fortran library](https://www.unidata.ucar.edu/downloads/netcdf/)

The two netCDF libraries are used for accessing the nc format files through fortran. Install netCDF-C library first. The model needs a head file named netcdf.inc, which can be found after the netCDF-Fortran library has been installed.

To install the netCDF-C library, go to the unzipped folder and read the INSTALL.md for instructions. I used
```
$ ./configure --prefix=/home/username/local --disable-netcdf-4
$ make check install
```
where _/home/username/local_ is the location you want to install the library. I think I don't need to disable netcdf-4.

To install the netCDF-FORTRAN library, go to the unzipped folder and run
```
$ ./configure --prefix=/home/username/local
$ make check install
```

Check netcdf.inc is under /home/username/local/include. Change the INCLDIR in Macros.C90/Macros.J90 accordingly.

## Debug and modify the code
This model was run on NCAR CRAY server. There are some machine specific code and configurations, such as Mass Storage System (MSS), MSread/MSwrite, and keyword CRAY in the code. Modify the Macros.C90 or Macros.J90 file. Modify the Makefile. Also, debug the source code as well.

In my own model (the model I modify and run on my laptop), some folders are created to store the code and the log. To run it, run
```
~/pac3.e/csh -f test.01.nqs
```

In Macros file, '-L' gives the lib location, and '-l' gives the name of the lib file. For example, when '-l netcdff' is given, the compiler will look for a lib file named libnetcdff.a in the folder given by '-L'. 'a' stands for archive, which means that this is a static lib. Note that the lib file name is libnetcdff.a not libnetcdf.a, otherwise the compiler cann't find the netcdf related subroutines.

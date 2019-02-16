#
#  elcirc/selfe makefile
#  Mike Zulauf 10/15/2006
#

#
# the executable and source code file names (without extensions)
#
EXEC       = amd_zelfe1_5k6
MAIN       = elfe1_5k6
DSRC2      = dsrc2c
SFLUX      = sflux_9c
LAPACK     = lap

#
# specify different environment and/or preprocessing options
# (present environment options include XEON, AMD32, AMD64, and MAC)
#
ENV         = AMD32
#ENV         = XEON
USE_SFLUX   = yes
USE_NETCDF  = yes
USE_GOTM = yes
SELFE       = yes
#PREC_EVAP   = yes
#MM5         = yes

#
# if USE_SFLUX not set, then make sure USE_NETCDF is not set as well
#
ifndef USE_SFLUX
  USE_NETCDF   =
  PREC_EVAP   = 
endif

#
# options for different environments
# (compiler, linker, compilation flags, libraries, etc)
#

# Intel x86/Linux/32
ifeq ($(ENV),XEON)
  FC = ifort #ifc
  LN = $(FC)
#  NO_TR_15581 = yes
#  FFLAGS = -g -assume byterecl
  FFLAGS = -O3 -assume byterecl
  LFLAGS = -Bstatic
  FPP_FLAGS = -fpp
#  LIBS = -llapack -lblas -L/usr/lib/gcc-lib/i386-redhat-linux/3.2/ -lg2c
#  LIBS = $MKLPATH/libmkl_lapack.a $MKLPATH/libmkl_ia32.a $MKLPATH/libguide.a -lpthread
#  MKLPATH='/opt/intel/mkl70cluster/lib/32/'
  ifdef USE_NETCDF
    LIBS := $(LIBS) -L/usr/local/netcdf/lib -lnetcdf
    FFLAGS := $(FFLAGS) -I/usr/local/netcdf/include
  endif
  ifdef USE_GOTM
    LIBS := $(LIBS) -L/home/users/yinglong/GOTM/gotm-3.2.5/32bit/gotm-3.2.5/lib/IFORT/ -lturbulence_prod  -lutil_prod
    FFLAGS := $(FFLAGS) -I/home/users/yinglong/GOTM/gotm-3.2.5/32bit/gotm-3.2.5/modules/IFORT/
  endif
endif

# AMD Opteron/Linux/32
ifeq ($(ENV),AMD32)
  FC = ifort
  LN = $(FC)
#  NO_TR_15581 = yes
#  FFLAGS = -g -assume byterecl
  FFLAGS = -O3 -assume byterecl
  LFLAGS = 
  FPP_FLAGS = -fpp
  #LIBS = -L/usr/lib64/ -llapack
  ifdef USE_NETCDF
    LIBS := $(LIBS) -L/usr/local/netcdf/lib -lnetcdf
    FFLAGS := $(FFLAGS) -I/usr/local/netcdf/include
  endif
  ifdef USE_GOTM
#Locations of turbulence libraries and turburlence.mod
    LIBS := $(LIBS) -L/home/users/yinglong/GOTM/gotm-3.2.5/lib/IFORT/ -lturbulence_prod -lutil_prod
    FFLAGS := $(FFLAGS) -I/home/users/yinglong/GOTM/gotm-3.2.5/modules/IFORT/
  endif
endif

# AMD Opteron/Linux/64
ifeq ($(ENV),AMD64)
  FC = ifort
  LN = $(FC)
#  NO_TR_15581 = yes
  FFLAGS = -O3 -assume byterecl -mcmodel=large
  LFLAGS = -i_dynamic
  FPP_FLAGS = -fpp
  LIBS = -L/usr/lib64/ -llapack
  ifdef USE_NETCDF
    LIBS := $(LIBS) -L/usr/local/netcdf/lib -lnetcdf
    FFLAGS := $(FFLAGS) -I/usr/local/netcdf/include
  endif
  ifdef USE_GOTM
    LIBS := $(LIBS) -L/home/users/yinglong/GOTM/gotm-3.2.5/lib/IFORT/ -lturbulence_prod -lutil_prod
    FFLAGS := $(FFLAGS) -I/home/users/yinglong/GOTM/gotm-3.2.5/modules/IFORT/
  endif
endif

# Mac
#ifeq ($(ENV),MAC)
#  FC = g95
#  LN = $(FC)
##  NO_TR_15581 = yes
#  FFLAGS = -O -Wall
#  LFLAGS = 
#  FPP_FLAGS = -cpp
#  ifdef USE_NETCDF
#    LIBS = -L/sw/lib -lnetcdf
#    FFLAGS := $(FFLAGS) -I/sw/include
#  endif
#endif

#
# set preprocessing options (and combine with FFLAGS)
#

ifdef SELFE
  FPP_FLAGS := $(FPP_FLAGS) -DSELFE
endif

ifdef USE_NETCDF
  FPP_FLAGS := $(FPP_FLAGS) -DUSE_NETCDF
endif

ifdef USE_GOTM
  FPP_FLAGS := $(FPP_FLAGS) -DUSE_GOTM
endif

ifdef USE_SFLUX
  FPP_FLAGS := $(FPP_FLAGS) -DUSE_SFLUX
endif

ifdef PREC_EVAP
  FPP_FLAGS := $(FPP_FLAGS) -DPREC_EVAP
endif

ifdef MM5
  FPP_FLAGS := $(FPP_FLAGS) -DMM5
endif

ifdef NO_TR_15581
  FPP_FLAGS := $(FPP_FLAGS) -DNO_TR_15581
endif

FFLAGS := $(FPP_FLAGS) $(FFLAGS)

#
# comment out the following line if you don't want the build dependent
# upon the makefile, otherwise use the name of the makefile
#
#MAKEFILE = Makefile

#
# the object files
#
ifdef USE_SFLUX
  OBJS =  $(MAIN).o $(DSRC2).o $(SFLUX).o $(LAPACK).o
else
  OBJS =  $(MAIN).o $(DSRC2).o $(LAPACK).o
endif

#
# the actual build commands
#
$(EXEC): $(OBJS) $(MAKEFILE)
	$(LN) $(LFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

$(MAIN).o: $(MAIN).f90 $(MAKEFILE)
	$(FC) $(FFLAGS) -c $(MAIN).f90
$(DSRC2).o: $(DSRC2).f90 $(MAKEFILE)
	$(FC) $(FFLAGS) -c $(DSRC2).f90
$(SFLUX).o: $(SFLUX).f90 $(MAKEFILE)
	$(FC) $(FFLAGS) -c $(SFLUX).f90
$(LAPACK).o: $(LAPACK).f $(MAKEFILE)
	$(FC) $(FFLAGS) -c $(LAPACK).f

#
# how to clean _all_ up
#
clean:
	rm -f *.o *.mod $(EXEC)


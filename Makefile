PROG = IBIsCO

SRCS =    moduleparsing.f90 module_var.F90  \
          main.F90 WRITETP.f90 MAPS.f90  FTABLE.f90 WRITETRJ.f90 \
          AVERAGE.f90 ALLOCATEVAR.f90 SCALEBP.f90 MOMENTUM.f90 OUTPUT.F90 MOVE.F90 \
          UNIT.f90 COMVEL.f90 RDCONTROL.f90 \
          RDCOOR.f90 RDINTERACT.F90 PARSE.f90 SETLIS.f90 SHIFT.f90 \
          RDGAUSSIAN.f90 BONDTABLE.f90 ANGLETABLE.f90 LINKS.f90 \
          SCALEV.f90  \
          config.f90 WRITEPSF.f90\
          analysis.f90 SHAKE.f90 \
          ibi-preprocess.h RDVIRTUAL.f90 VIRTUAL_DEF.f90 \
          MAKE_LISTS.f90 DISTRIBUTE_VSFORCE.f90 NEW_FORCE.F90 NEW_LOOP.F90 \
          NONBONDED_FORCE.F90 UPDATE_NEIGHBOURLIST.f90 NEW_NEIGHBOUR_WITHLIST.f90 \
          BUILD_CONNECTIVITY.f90 BONDED_FORCE.F90 ALLOCATEVAR2.f90 \
          NEW_NEIGHBOUR_NOLIST.f90 \
          TIMING.F90 \
	  RDVIRTBONDS.f90 RDVIRTANGLES.f90

OBJS =    moduleparsing.o module_var.o \
          main.o WRITETP.o MAPS.o   FTABLE.o WRITETRJ.o \
          AVERAGE.o ALLOCATEVAR.o SCALEBP.o MOMENTUM.o OUTPUT.o MOVE.o \
          UNIT.o COMVEL.o  RDCONTROL.o RDCOOR.o \
          RDINTERACT.o PARSE.o SETLIS.o SHIFT.o RDGAUSSIAN.o \
          BONDTABLE.o ANGLETABLE.o LINKS.o SCALEV.o   \
          config.o WRITEPSF.o\
          analysis.o SHAKE.o RDVIRTUAL.o VIRTUAL_DEF.o \
          MAKE_LISTS.o DISTRIBUTE_VSFORCE.o NEW_FORCE.o NEW_LOOP.o \
          NONBONDED_FORCE.o UPDATE_NEIGHBOURLIST.o NEW_NEIGHBOUR_WITHLIST.o \
          BUILD_CONNECTIVITY.o BONDED_FORCE.o ALLOCATEVAR2.o \
          NEW_NEIGHBOUR_NOLIST.o \
          TIMING.o \
	  RDVIRTBONDS.o RDVIRTANGLES.o

.SUFFIXES: $(SUFFIXES) .f90 .F90

#F90 = pgf90
F90 = gfortran
#F90 = ifort
# ************ Profiling with gfortran *************
#LDFLAGS= -pg
#F90FLAGS  = -g -pg
#*****************************************
#Fast ifort options
#F90FLAGS = -openmp -openmp-report1 -O2
#LDFLAGS = -openmp -O2
#Debugging with gfortran
#F90FLAGS = -pg -fbounds-check -fbacktrace -O0
#LDFLAGS = -g -pg -fbounds-check -fbacktrace -O0
#F90FLAGS = -Mbounds -g pt=px-Bstatic
#Fast gfortran options
F90FLAGS = -O2 -fopenmp
LDFLAGS = -O2 -fopenmp
#F90FLAGS = -g -Wall -Wextra -Wconversion
#LIBS= -L/home/nicodemo/bin/lib -lmpack1  
# l'opzione -c compila e assembla ma non linka


all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

%.o: %.mod

#.f90.o:
#	$(F90) $(F90FLAGS) -c $<

#.F90.o:
#	$(F90) $(F90FLAGS) -c $<

%.o: %.f90 
	$(F90) $(F90FLAGS) -c $< -o $@
%.o: %.F90 
	$(F90) $(F90FLAGS) -c $< -o $@

clean:
	rm -f $(PROG) $(OBJS) *.mod

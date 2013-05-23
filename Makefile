PROG = IBIsCO-openmp-test

SRCS =    moduleparsing.f90 module_var.F90  module_RNEMD.f90  module_PAIR.f90   \
          main.F90 WRITETP.f90 MAPS.f90  FTABLE.f90 WRITETRJ.f90 \
          AVERAGE.f90 ALLOCATEVAR.f90 SCALEBP.f90 MOMENTUM.f90 OUTPUT.F90 MOVE2.f90 MOVE.f90 \
          FPANGLE.f90 FPTOR.f90 UNIT.f90 COMVEL.f90 RDCONTROL.f90 \
          RDCOOR.f90 RDINTERACT.f90 PARSE.f90 SETLIS.f90 SHIFT.f90 \
          RDGAUSSIAN.f90 BONDTABLE.f90 ANGLETABLE.f90 LINKS.f90 \
          SCALEV.f90  \
          RDVISC.f90 VISC_NEMD.f90  VISC_ATOM.f90  VISC_PROFILE.f90  \
          LINEAR.f90 PBCPOISEUILLE.f90 LOWEAND.f FDR.f FDR_STEP.f DPDVCMTP.f90 DPDVCMTP2.f90 \
          FDR_ADDITION.f90 FDR_STEP_ADDITION.f90 LA_ADDITION.f90 DPDV.f90 DPDV2.f90 DPDV22.f90 DPDVO.f90 DPDVO2.f90 DPDR.f90\
          PROFILEPPF.f90 PROFILERNEMD.f90 DPDVCMTPPV.f90 MOVERNEMD.f90 FPOUTPLANE.F90 \
	  NON_BOND_ARRAY_BEAD.f90 \
          config.f90 \
          NON_BOND_ARRAY_BEAD_SEC.f90 \
          analysis.f90 SHAKE.f90 \
          ibi-preprocess.h RDVIRTUAL.f90 VIRTUAL_DEF.f90 \
          MAKE_LISTS.f90 DISTRIBUTE_VSFORCE.f90 NEW_FORCE.f90 NEW_LOOP.f90 \
          NONBONDED_FORCE.F90 UPDATE_NEIGHBOURLIST.f90 NEW_NEIGHBOUR_WITHLIST.f90 \
          BUILD_CONNECTIVITY.f90 BONDED_FORCE.f90 ALLOCATEVAR2.f90


OBJS =    moduleparsing.o module_var.o  module_RNEMD.o \
          module_PAIR.o   main.o WRITETP.o MAPS.o   FTABLE.o WRITETRJ.o \
          AVERAGE.o ALLOCATEVAR.o SCALEBP.o MOMENTUM.o OUTPUT.o MOVE2.o MOVE.o \
          FPANGLE.o FPTOR.o UNIT.o COMVEL.o  RDCONTROL.o RDCOOR.o \
          RDINTERACT.o PARSE.o SETLIS.o SHIFT.o RDGAUSSIAN.o \
          BONDTABLE.o ANGLETABLE.o LINKS.o SCALEV.o   \
          RDVISC.o   VISC_NEMD.o  VISC_ATOM.o   VISC_PROFILE.o  \
          LINEAR.o PBCPOISEUILLE.o LOWEAND.o  FDR.o FDR_STEP.o DPDVCMTP.o DPDVCMTP2.o \
          FDR_ADDITION.o FDR_STEP_ADDITION.o LA_ADDITION.o DPDV.o DPDV2.o DPDV22.o DPDVO.o DPDVO2.o DPDR.o \
          PROFILEPPF.o PROFILERNEMD.o DPDVCMTPPV.o MOVERNEMD.o FPOUTPLANE.o \
          config.o \
          NON_BOND_ARRAY_BEAD_SEC.o \
          analysis.o SHAKE.o RDVIRTUAL.o VIRTUAL_DEF.o \
          MAKE_LISTS.o DISTRIBUTE_VSFORCE.o NEW_FORCE.o NEW_LOOP.o \
          NONBONDED_FORCE.o UPDATE_NEIGHBOURLIST.o NEW_NEIGHBOUR_WITHLIST.o \
          BUILD_CONNECTIVITY.o BONDED_FORCE.o ALLOCATEVAR2.o

.SUFFIXES: $(SUFFIXES) .f90 .F90

#F90 = pgf90
F90 = gfortran
#F90 = ifort
# ************ Profiling with gfortran *************
#LDFLAGS= -pg
#F90FLAGS  = -g -pg
#*****************************************
#F90FLAGS = -openmp -openmp-report1 -O2
#LDFLAGS = -openmp -O2
#F90FLAGS = -Mbounds -g pt=px-Bstatic
F90FLAGS = -g -fopenmp -fbacktrace -fbounds-check
LDFLAGS = -g -fopenmp -fbacktrace -fbounds-check
#F90FLAGS = -g -Wall -Wextra -Wconversion
#F90FLAGS  = -g -fbounds-check 
#F90FLAGS = -O3
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

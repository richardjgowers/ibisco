PROG = IBIsCO-openmp-test

SRCS =    moduleparsing.f90 module_var.F90  module_RNEMD.f90  module_PAIR.f90   \
          main.F90 LOOP.f90 LOOPDPD.f90  WRITETP.f90 MAPS.f90  FTABLE.f90 WRITETRJ.f90 \
          AVERAGE.f90 ALLOCATEVAR.f90 SCALEBP.f90 MOMENTUM.f90 OUTPUT.F90 MOVE2.f90 MOVE.f90 FPBOND.f90 \
          FPANGLE.f90 FPTOR.f90 FORCE.F90 UNIT.f90 COMVEL.f90 RDCONTROL.f90 \
          RDCOOR.f90 RDINTERACT.f90 PARSE.f90 SETLIS.f90 NON_BOND_ARRAY.f90 SHIFT.f90 \
          RDGAUSSIAN.f90 BONDTABLE.f90 ANGLETABLE.f90 NEIGHBOR_NOLIST.f90 NEIGHBOR_WITHLIST.f90 LINKS.f90 \
          SCALEV.f90  \
          RDVISC.f90 VISC_NEMD.f90  VISC_ATOM.f90  VISC_PROFILE.f90  \
          LINEAR.f90 PBCPOISEUILLE.f90    FORCE2.f90 LOWEAND.f FDR.f FDR_STEP.f DPDVCMTP.f90 DPDVCMTP2.f90 \
          FDR_ADDITION.f90 FDR_STEP_ADDITION.f90 LA_ADDITION.f90 DPDV.f90 DPDV2.f90 DPDV22.f90 DPDVO.f90 DPDVO2.f90 DPDR.f90\
          PROFILEPPF.f90 PROFILERNEMD.f90 DPDVCMTPPV.f90 MOVERNEMD.f90 FPOUTPLANE.F90 \
	      VIRT_NEIGHBOR_WITHLIST.f90 VIRT_LINKS.f90 VIRT_NEIGHBOR_WITHLIST_COM.f90\
	      VIRT_SEC_FORCE_COM.F90 NON_BOND_ARRAY_BEAD.f90 VIRT_LINKS_COM.f90 \
	      MTS.f90 LOOP_HYBR.f90 FORCE_HYBR.f90 FORCE_HYBR_MTS.f90 VIRT_SEC_FORCE_COM_MTS.F90 \
	      VIRT_FORCE_SEC_MTS.F90 LOOP_HYBR_MTS3_COM.f90 VIRT_FORCE_COM_MTS.F90 VIRT_FORCE_COM.F90 \
          config.f90 LOOP_HYBR_MTS3.f90 FORCE_HYBR_MTS_COM.f90\
          NON_BOND_ARRAY_BEAD_SEC.f90 NON_BOND_FIRST.f90 LOOP_HYBR_MTS4.f90 LOOP_HYBR_MTS5.f90 \
          FORCE_HYBR_PRE-MTS.f90 LOOP_HYBR_MTS4_COM.f90 FORCE_HYBR_PRE-MTS_COM.f90\
          LOOP_HYBR_MTS5_COM.f90 LOOP_HYBR_COM.f90 FORCE_HYBR_COM.f90 analysis.f90 SHAKE.f90 \
          HFPBOND.f90 HFPANGLE.f90 HFPTOR.f90 HAVERAGE.f90 ibi-preprocess.h RDVIRTUAL.f90 VIRTUAL_DEF.f90 \
          MAKE_LISTS.f90 ALLOCATEVAR2.f90

## Old VIRTUAL_SITE_COM.f90 VIRTUAL_SITE.f90


OBJS =    moduleparsing.o module_var.o  module_RNEMD.o \
          module_PAIR.o   main.o LOOP.o LOOPDPD.o WRITETP.o MAPS.o   FTABLE.o WRITETRJ.o \
          AVERAGE.o ALLOCATEVAR.o SCALEBP.o MOMENTUM.o OUTPUT.o MOVE2.o MOVE.o FPBOND.o \
          FPANGLE.o FPTOR.o FORCE.o UNIT.o COMVEL.o  RDCONTROL.o RDCOOR.o \
          RDINTERACT.o PARSE.o SETLIS.o NON_BOND_ARRAY.o SHIFT.o RDGAUSSIAN.o \
          BONDTABLE.o ANGLETABLE.o NEIGHBOR_NOLIST.o NEIGHBOR_WITHLIST.o  LINKS.o SCALEV.o   \
          RDVISC.o   VISC_NEMD.o  VISC_ATOM.o   VISC_PROFILE.o  \
          LINEAR.o PBCPOISEUILLE.o FORCE2.o LOWEAND.o  FDR.o FDR_STEP.o DPDVCMTP.o DPDVCMTP2.o \
          FDR_ADDITION.o FDR_STEP_ADDITION.o LA_ADDITION.o DPDV.o DPDV2.o DPDV22.o DPDVO.o DPDVO2.o DPDR.o \
          PROFILEPPF.o PROFILERNEMD.o DPDVCMTPPV.o MOVERNEMD.o FPOUTPLANE.o VIRT_NEIGHBOR_WITHLIST.o \
	      VIRT_LINKS.o VIRT_NEIGHBOR_WITHLIST_COM.o \
	      VIRT_SEC_FORCE_COM.o NON_BOND_ARRAY_BEAD.o VIRT_LINKS_COM.o \
          LOOP_HYBR.o FORCE_HYBR.o FORCE_HYBR_MTS.o VIRT_SEC_FORCE_COM_MTS.o VIRT_FORCE_SEC_MTS.o LOOP_HYBR_MTS3_COM.o \
          VIRT_FORCE_COM_MTS.o VIRT_FORCE_COM.o LOOP_HYBR_MTS4.o LOOP_HYBR_MTS5.o\
          config.o LOOP_HYBR_MTS3.o FORCE_HYBR_MTS_COM.o\
          NON_BOND_ARRAY_BEAD_SEC.o NON_BOND_FIRST.o FORCE_HYBR_PRE-MTS_COM.o\
          FORCE_HYBR_PRE-MTS.o LOOP_HYBR_MTS4_COM.o LOOP_HYBR_MTS5_COM.o LOOP_HYBR_COM.o FORCE_HYBR_COM.o \
          analysis.o SHAKE.o HFPBOND.o HFPANGLE.o HFPTOR.o HAVERAGE.o RDVIRTUAL.o VIRTUAL_DEF.o \
          MAKE_LISTS.o ALLOCATEVAR2.f90

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
F90FLAGS = -g -fopenmp -fbacktrace -fbounds-check -O2
LDFLAGS = -g -fopenmp -fbacktrace -fbounds-check -O2
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

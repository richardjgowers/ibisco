#ifndef DOXY_SKIP
#include "ibi-preprocess.h"
#endif
! Header page for doxygen

!> \mainpage
!!
!! \section IBIsCO Documentation!
!!
!! Each file in the source code has its own page, click around to explore the source code.
!!
!! Everything starts in:
!! main.F90
!! 
!! \subsection Inputs
!!
!! Read about the required input files in 
!!
!! RDCONTROL.f90 \n 
!! RDINTERACT.f90 \n
!! RDCOOR.f90 \n
!!
!! Read about hybrid simulation input files in
!!
!! RDVIRTUAL.f90 \n
!! RDVIRTBONDS.f90 \n
!! RDVIRTANGLES.f90 \n 
!!
!! \subsection Outputs
!!
!! Read about how the results are recorded
!!
!! OUTPUT.F90 \n
!! WRITETRJ.f90 \n
!!
!! Or check out the molecular dynamics loop in
!!
!! NEW_LOOP.F90

!> @file
!> @brief IBIsCO MD program
!!
!> @details This is the trunk of the program
!!
!! \section main_sec1 Reading program inputs
!! UNIT.f90 \n
!! RDCONTROL.f90 \n
!! ALLOCATEVAR.f90 \n
!! RDINTERACT.f90 \n
!! RDCOOR.f90 \n
!! ALLOCATEVAR2.f90 \n
!! RDVIRTUAL.f90 \n
!!
!! \section main_sec2 Preparing data structures
!! 
!! MAKE_LISTS.f90 \n
!! VIRTUAL_DEF.f90 \n
!! SHIFT.f90 \n
!! FTABLE.f90 \n
!! COMVEL.f90 \n
!! SCALEV.f90 \n
!! SETLIS.f90 \n
!! RDVIRTANGLES.f90 \n
!! BUILD_CONNECTIVITY.f90 \n
!! RDVIRTBONDS.f90 \n
!!
!! \section main_sec3 The MD loop
!!
!! NEW_LOOP.F90 \n

PROGRAM IBISCO

  USE MODULEPARSING

  USE VAR

  IMPLICIT NONE

  INTEGER :: I, A !< Used as counter variables

  WRITE(*,*) 'IBIsCO TIME! Revision 78'
  OPEN ( 115 , FILE = 's-md.out')
  OPEN ( 116 , FILE = 's-md.tp')
  OPEN ( 113 , FILE = 's-md.trj', form='UNFORMATTED', access='SEQUENTIAL')
#ifdef TIMING_ON
  OPEN ( 44  , FILE = 'timing.out')
#endif

  WRITE( 115, *)'IBIsCO Revision 78'
  OPEN (1, FILE='ERROR')
  WRITE(*,*)
  ISTOP = 0

  !DEFINE REDUCED UNITS
  CALL UNIT() !> \file UNIT.f90

  !Read control parameters
  NVIRTA = 0 !By default no virtual sites
  CALL RDCONTROL()
  IF (ISTOP == 1) STOP 'Failed in RDCONTROL'

  !     ALLOCATE SOME VARIABLES
  CALL ALLOCATEVAR()

  !READ PARAM FILE
  CALL RDINTERACT()
  IF (ISTOP == 1) STOP 'Failed in RDINTERACT'

  !Read coordinate file
  CALL  RDCOOR()
  IF (ISTOP == 1) STOP 'Failed in RDCOOR'

  CALL ALLOCATEVAR2()

  NCOARSE = 0
  IF(IBRDESCR .EQ. 0) THEN	
     !Read virtual file			
     CALL RDVIRTUAL()

     IF(ISTOP .eq. 1) STOP 'Failed at RDVIRTUAL'
     !Make linked lists of atoms or bead/VS
  END IF

  IF(IBRDESCR .EQ. 0) THEN
     CALL MAKE_LISTS()
     IF(ISTOP .eq. 1) STOP 'Failed at MAKE_LISTS'
  ELSE
     !Nonhybrid settings, single list containing all atoms
     NUMATOMS = NATOMS
     ALLOCATE(ATOM(NUMATOMS))
     DO I=1,NATOMS
        ATOM(I) = I
     END DO
  END IF

  IF(IBRDESCR .EQ. 0) THEN
     ! Calculate centers of VS, either COM or adopt an atom's coords
     CALL VIRTUAL_DEF()
  END IF

  !     SHIFT THE ATOMS INSIDE THE BOX
  CALL SHIFT ()

  MAXNAB_ATOM = 1000 !Max number of neighbours for a single atom
  ALLOCATE(LIST_ATOM(MAXNAB_ATOM,NUMATOMS),CELL_ATOM(NUMATOMS),LCLIST_ATOM(NUMATOMS))
  IF(IBRDESCR .eq. 0) THEN
     MAXNAB_BEAD = 1000
     ALLOCATE(LIST_BEAD(MAXNAB_BEAD,NCOARSE),CELL_BEAD(NCOARSE),LCLIST_BEAD(NCOARSE))
  END IF

  !MAKE TABLE FORCES
  CALL FTABLE ()

  !     GIVE THE VELOCITIES IF THE INITIAL TIME IS ZERO
  IF (RESTART == 1) THEN
     CALL COMVEL ()
     CALL SCALEV ()
  END IF

  !     CALCULATE THE INITIAL TEMPERATURE AND KINETIC ENERGY
  EK(1) = 0.0
  EK(2) = 0.0
  WRITE(*,*) NUMATOMS, NUMBEADS
  MKTEMP_ATOM = 2.0D0 / REAL(3.0D0 * (NUMATOMS - 1))
  MKTEMP_BEAD = 2.0D0 / REAL(3.0D0 * (NUMBEADS - 1))
  DO A = 1, NUMATOMS
     I = ATOM(A)
     EK(1) = EK(1) + MASS(ITYPE(I))*(VX(I)**2.0D0 + VY(I)**2.0D0 + VZ(I)**2.0D0)
     TOTMASS = TOTMASS + MASS(ITYPE(I))
  END DO

  DO A=1,NUMBEADS
     I = BEAD(A)
     EK(2) = EK(2) + MASS(ITYPE(I))*(VX(I)**2.0D0 + VY(I)**2.0D0 + VZ(I)**2.0D0)
     TOTMASS = TOTMASS + MASS(ITYPE(I))
  END DO

  EK(1) = 0.5D0 * EK(1)
  EK(2) = 0.5D0 * EK(2)
  TEMP = SUM(EK) * MKTEMP	
  TEMP_ATOM = EK(1) * MKTEMP_ATOM
  TEMP_BEAD = EK(2) * MKTEMP_BEAD

  !     CREATE LISTS OF BONDS, ANGLES AND TORSIONS
  CALL SETLIS
  IF (ISTOP == 1) STOP 'Failed in SETLIS'

  !Read virtangles
  IF(IBRDESCR .eq. 0) THEN
     CALL RDVIRTANGLES()
     IF(ISTOP .eq. 1) STOP 'Failed in RDVIRTANGLES'
  END IF

  CONNECTIONS = 0
  CONNECTED_TO = 0

  CALL BUILD_CONNECTIVITY(NUMATOMS,ATOM,NONBEXC_ATOM)
  IF(ISTOP .eq. 1) STOP 'Failed in BUILD_CONNECTIVITY'

  IF(IBRDESCR .eq. 0) THEN
     CALL BUILD_CONNECTIVITY(NCOARSE,BEAD,NONBEXC_BEAD)
     IF(ISTOP .eq. 1) STOP 'Failed in BUILD_CONNECTIVITY beads'
     CALL RDVIRTBONDS()
     IF(ISTOP .EQ. 1) STOP 'Failed in RDVIRTBONDS'
  END IF

#ifdef REPORT_EXCLUSIONS
  CALL REPORT_EXCLUSIONS() ! Reports all nonbonded exclusions
#endif

  CALL MAPS (MAP_ATOM,MAPSIZE_ATOM &
       , NCELLX_ATOM, NCELLY_ATOM, NCELLZ_ATOM) 
  CALL LINKS (HEAD_ATOM,MAXNUMCELL_ATOM,ATOM,NUMATOMS,CELL_ATOM &
       , NCELLX_ATOM, NCELLY_ATOM, NCELLZ_ATOM, LCLIST_ATOM)

  IF(IBRDESCR .eq. 0) THEN
     CALL MAPS (MAP_BEAD,MAPSIZE_BEAD &
          , NCELLX_BEAD, NCELLY_BEAD, NCELLZ_BEAD) 
     CALL LINKS (HEAD_BEAD,MAXNUMCELL_BEAD,BEAD,NCOARSE,CELL_BEAD &
          , NCELLX_BEAD, NCELLY_BEAD, NCELLZ_BEAD, LCLIST_BEAD)
  END IF

  CALL UPDATE_NEIGHBOURLIST()

  !     WRITE THE TOPOLOGY FILE 
  CALL WRITEPSF()

  !!*************************************************
  !*********** START OF MAIN LOOP*****************
  !*************************************************

  WRITE(*,*)
  WRITE(*,*) 'Beginning simulation'
  WRITE(*,*)

  CALL NEW_LOOP()


  !*************************************************
  !***********  END OF MAIN LOOP********************
  !*************************************************

  CLOSE (201)
  CLOSE (202)

  Write(*,*)'**************************************'
  Write(*,*)'********* END OF DYNAMICS ************'
  Write(*,*)'¤.¸.¤.¸.¤°´¯`´¯`°¤.¸.¤.¸.¤°´¯`´¯`°¤.¸'
  Write(*,*)'  (|_(|(|_(|       (|_(|     '
  Write(*,*)'  (="',':','" (="',':','")      (="',':','")    '
  Write(*,*)'  (,('')('')     (,('')('')    (,('')('')  '
  Write(*,*)'¤.¸.¤.¸.¤°´¯`´¯`°¤.¸.¤.¸.¤°´¯`´¯`°¤.¸'
  Write(*,*)'**************************************'

END PROGRAM IBISCO

!> @file
!> @brief Allocates some arrays
!! 
!> @details Called from main.f90
!

SUBROUTINE ALLOCATEVAR ()

  USE VAR

  IMPLICIT NONE

  INTEGER :: MAXBONDS = 4 !Maximum number of bonds an atom can have

  !	PROGRAM CONSTANTS
  FAC = DT / TAUT
  PFAC = BETA*DT*LIMAVP/ (TAUP)

  !	CALCULATE THE TEMPERATURE BY KINETIC ENERGY (TEMP = EK*MKTEMP)
  MKTEMP = 2.0D0 / REAL( 3.0D0 * (NATOMS-1))

  !Reset average totals to 0 to begin new averaging window
  AV_TOT_E = 0.0
  AV_POT_E = 0.0
  AV_KIN_E = 0.0
  AV_V_NB = 0.0
  AV_V_BOND = 0.0
  AV_V_ANGLE = 0.0
  AV_V_TORSION = 0.0
  AV_V_OOP = 0.0
  AV_TEMP = 0.0
  AV_PRES = 0.0
  AV_BOX = 0.0
  AV_DENS = 0.0

  !            ALLOCATE(RX(NATOMS)) !Allocated in ALLOCATEVAR2 instead
  !            ALLOCATE(RY(NATOMS))
  !            ALLOCATE(RZ(NATOMS))

  !Force and position variables
  !
  !F and R both include the virtual sites appended on the end (NITEMS = NATOMS + NVIRTA)
  !ITYPE also includes NVIRTA

  ALLOCATE(SX(NITEMS), SY(NITEMS), SZ(NITEMS) )
  ALLOCATE(VX(NATOMS), VY(NATOMS), VZ(NATOMS) )
  ALLOCATE(VTX(NATOMS), VTY(NATOMS), VTZ(NATOMS) )
  ALLOCATE(FX(NITEMS), FY(NITEMS), FZ(NITEMS) )
  ALLOCATE(RX(NITEMS), RY(NITEMS), RZ(NITEMS) )
  ALLOCATE(FXNB(NATOMS), FYNB(NATOMS), FZNB(NATOMS) )

  FX = 0.0D0
  FY = 0.0D0
  FZ = 0.0D0

  RX = 0.0D0
  RY = 0.0D0
  RZ = 0.0D0

  ALLOCATE(CONNECTIONS(NITEMS))
  CONNECTIONS = 0
  ALLOCATE(CONNECTED_TO(NITEMS,MAXCONNECTIONS))
  CONNECTED_TO = 0
  ALLOCATE(ITYPE(NITEMS))
  ALLOCATE(MOL(NATOMS))
  ALLOCATE(NBONDS(NITEMS))
  ALLOCATE(NIJK(NITEMS))
  ALLOCATE(NOANGLEIJK(NITEMS))
  ALLOCATE(NIJKL(NITEMS))
  ALLOCATE(FNIJKL(NITEMS))
  ALLOCATE(NOOPIJKL(NITEMS))
  ALLOCATE(JBOND(NITEMS,MAXBONDS))


  ALLOCATE(SP(LIMAVP))

  !###############################################################################
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !| Allocate some variables for Hybrid systems and MTS |
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IF (IBRDESCR .eq. 0) THEN

     if( MTS_CHECK .EQ. 0)then
        ! VECTOR FOR THE Multiple Time Step
        ALLOCATE(FXi0(NATOMS))
        ALLOCATE(FYi0(NATOMS))
        ALLOCATE(FZi0(NATOMS))
        ALLOCATE(FXi1(NATOMS))
        ALLOCATE(FYi1(NATOMS))
        ALLOCATE(FZi1(NATOMS))
        ALLOCATE(FXii(NATOMS))
        ALLOCATE(FYii(NATOMS))
        ALLOCATE(FZii(NATOMS))
        ALLOCATE(FXiii(NATOMS))
        ALLOCATE(FYiii(NATOMS))
        ALLOCATE(FZiii(NATOMS))
        if(type_mts .eq. 4)then
           ALLOCATE(FXiiii(NATOMS))
           ALLOCATE(FYiiii(NATOMS))
           ALLOCATE(FZiiii(NATOMS))
        elseif(type_mts .eq. 5)then
           ALLOCATE(FXiiii(NATOMS))
           ALLOCATE(FYiiii(NATOMS))
           ALLOCATE(FZiiii(NATOMS))
        endif
     end if
  END IF

!     FXi0,   FYi0,   FZi0   = components of the force at time i-1
!     FXi1,   FYi1,   FZi1   = components of the force at time i
!     FXii,   FYii,   FZii   = components of the force at time i+1
!     FXiii,  FYiii,  FZiii  = components of the force at time i+2
!     FXiiii, FYiiii, FZiiii = components of the force at time i+3     


  RETURN

END SUBROUTINE ALLOCATEVAR

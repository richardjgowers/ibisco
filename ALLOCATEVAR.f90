!     Modified by Nicodemo
!     2011
!
!     Vector Allocated for the Multiple Time Step
!
!     FXi0,   FYi0,   FZi0   = components of the force at time i-1
!     FXi1,   FYi1,   FZi1   = components of the force at time i
!     FXii,   FYii,   FZii   = components of the force at time i+1
!     FXiii,  FYiii,  FZiii  = components of the force at time i+2
!     FXiiii, FYiiii, FZiiii = components of the force at time i+3     
!
	!	*********************************************************************************************
SUBROUTINE ALLOCATEVAR ()

  USE VAR
  USE PAIR
  USE RNEMD 

  IMPLICIT NONE

  INTEGER :: MAXBONDS = 4 !Maximum number of bonds an atom can have

	!       ******************************************************************

  !	PROGRAM CONSTANTS
  FAC = DT / TAUT
  PFAC = BETA*DT*LIMAVP/ (TAUP)
  !       PFAC = BETA*DT/ (TAUP)
  DT2   = DT / 2.0D0
  DTSQ2 = DT * DT2

  !	CALCULATE THE TEMPERATURE BY KINETIC ENERGY (TEMP = EK*MKTEMP)
  MKTEMP = 2.0D0 / REAL( 3.0D0 * (NATOMS-1))

  IRAV = 0
  NRAV = 0
  REK = 0.0D0
  RVBOND = 0.0D0
  RVANGLE = 0.0D0
  RVTOR = 0.0D0
  RVNBOND = 0.0D0
  RV = 0.0D0
  RE = 0.0D0
  RTEMP = 0.0D0
  RPT11 = 0.0D0
  RPT22 = 0.0D0
  RPT33 = 0.0D0
  RPT12 = 0.0D0
  RPT13 = 0.0D0
  RPT23 = 0.0D0
  RP = 0.0D0
  RBOXX = 0.0D0
  RBOXY = 0.0D0
  RBOXZ = 0.0D0
  RVOL = 0.0D0
  RDENS = 0.0D0

!            ALLOCATE(RX(NATOMS)) !Allocated in ALLOCATEVAR2 instead
!            ALLOCATE(RY(NATOMS))
!            ALLOCATE(RZ(NATOMS))

!Force and position variables
!
!F and R both include the virtual sites appended on the end (NITEMS = NATOMS + NVIRTA)
!ITYPE also includes NVIRTA

  ALLOCATE(SX(NATOMS), SY(NATOMS), SZ(NATOMS) )
  ALLOCATE(VX(NATOMS), VY(NATOMS), VZ(NATOMS) )
  ALLOCATE(VTX(NATOMS), VTY(NATOMS), VTZ(NATOMS) )
  ALLOCATE(FX(NITEMS), FY(NITEMS), FZ(NITEMS) )
  ALLOCATE(RX(NITEMS), RY(NITEMS), RZ(NITEMS) )
  ALLOCATE(FXNB(NATOMS), FYNB(NATOMS), FZNB(NATOMS) )
  ALLOCATE(VFXNB(NATOMS), VFYNB(NATOMS), VFZNB(NATOMS) )
  VFXNB = 0
  VFYNB = 0
  VFZNB = 0
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
  ALLOCATE(NBONDS(NATOMS))
  ALLOCATE(NIJK(NATOMS))
  ALLOCATE(NOANGLEIJK(NATOMS))
  ALLOCATE(NIJKL(NATOMS))
  ALLOCATE(FNIJKL(NATOMS))
  ALLOCATE(NOOPIJKL(NATOMS))
  ALLOCATE(JBOND(NATOMS,MAXBONDS))



            
  IF (DPDINPUT.EQ.1) THEN
     ALLOCATE(VOX(NATOMS))
     ALLOCATE(VOY(NATOMS))
     ALLOCATE(VOZ(NATOMS))
  ENDIF

  IF (DPDINPUT.EQ.1.OR.LAINPUT.EQ.1) THEN
     ALLOCATE(DPDPOINT(NATOMS+1))
     ALLOCATE(DPDLIST(NATOMS*50))
  ENDIF
             
  ALLOCATE(STEK(LIMRAV))
  ALLOCATE(STE(LIMRAV))
  ALLOCATE(STEMP(LIMRAV))
  ALLOCATE(STVBOND(LIMRAV))
  ALLOCATE(STVNBOND(LIMRAV))
  ALLOCATE(STVTOR(LIMRAV))
  ALLOCATE(STVOOP(LIMRAV))
  ALLOCATE(STVANGLE(LIMRAV))
  ALLOCATE(STV(LIMRAV))
  ALLOCATE(STPT11(LIMRAV))
  ALLOCATE(STPT22(LIMRAV))
  ALLOCATE(STPT33(LIMRAV))
  ALLOCATE(STPT12(LIMRAV))
  ALLOCATE(STPT13(LIMRAV))
  ALLOCATE(STPT23(LIMRAV))
  ALLOCATE(STP(LIMRAV))
  ALLOCATE(SP(LIMAVP))

  ALLOCATE(SBOXX(LIMRAV))
  ALLOCATE(SBOXY(LIMRAV))
  ALLOCATE(SBOXZ(LIMRAV))
  ALLOCATE(SVOL(LIMRAV))
  ALLOCATE(SDENS(LIMRAV))

!###############################################################################
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| Allocate some variables for Hybrid systems and MTS |
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

IF (IBRDESCR .eq. 0) THEN
    ALLOCATE(STVBOND_At(LIMRAV))
    ALLOCATE(STVNBOND_At(LIMRAV))
    ALLOCATE(STVTOR_At(LIMRAV))
    ALLOCATE(STVOOP_At(LIMRAV))
    ALLOCATE(STVANGLE_At(LIMRAV))
    ALLOCATE(STV_At(LIMRAV))
    ALLOCATE(STVBOND_MIX(LIMRAV))
    ALLOCATE(STVNBOND_MIX(LIMRAV))
    ALLOCATE(STVTOR_MIX(LIMRAV))
    ALLOCATE(STVOOP_MIX(LIMRAV))
    ALLOCATE(STVANGLE_MIX(LIMRAV))
    ALLOCATE(STV_MIX(LIMRAV))
    ALLOCATE(STVBOND_CG(LIMRAV))
    ALLOCATE(STVNBOND_CG(LIMRAV))
    ALLOCATE(STVTOR_CG(LIMRAV))
    ALLOCATE(STVOOP_CG(LIMRAV))
    ALLOCATE(STVANGLE_CG(LIMRAV))
    ALLOCATE(STV_CG(LIMRAV))

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
          ALLOCATE(FXv(NATOMS))
          ALLOCATE(FYv(NATOMS))
          ALLOCATE(FZv(NATOMS))
       endif
    end if
END IF

!###############################################################################

! Allocate variables for DPD-thermostat

!         IF ( (STANDARD_DPDINPUT .EQ.  1) .OR. (TRANDPD_INPUT .EQ.1)) THEN
!                 ALLOCATE(VXO(NATOMS))
!                 ALLOCATE(VYO(NATOMS))
!                 ALLOCATE(VZO(NATOMS))
!                 ALLOCATE(DPD_POINT(NATOMS))
!                 ALLOCATE(DPD_LIST(MAXNAB))
!         ENDIF

!         IF (  STANDARD_DPDINPUT .EQ.  1)  THEN
!              ALLOCATE ( FDX(NATOMS) ) 
!              ALLOCATE ( FDY(NATOMS) )
!              ALLOCATE ( FDZ(NATOMS) )
!              ALLOCATE ( FRANDOMX(NATOMS))
!              ALLOCATE ( FRANDOMY(NATOMS))
!              ALLOCATE ( FRANDOMZ(NATOMS))
!         ENDIF
 
!           IF (TRANDPD_INPUT .EQ.  1)   THEN
!              ALLOCATE ( TRAN_FDX (NATOMS ) ) 
!              ALLOCATE ( TRAN_FDY (NATOMS ) )
!              ALLOCATE ( TRAN_FDZ (NATOMS ) )
!              ALLOCATE ( TRAN_FRANDOMX (NATOMS ) )
!              ALLOCATE ( TRAN_FRANDOMY (NATOMS ) )
!              ALLOCATE ( TRAN_FRANDOMZ (NATOMS ) )
!           ENDIF
           
      IF (VISCINPUT == 1) THEN
            ALLOCATE (AM(NATOMS))
!           ALLOCATE (VXPROF(NUMSLAB))
            ALLOCATE (VXMEAN_SLAB (NUMSLAB))
            ALLOCATE (VX_SLAB(NUMSLAB))
            ALLOCATE (NUMA_SLAB(NUMSLAB))
            ALLOCATE (SLAB_DENS(NUMSLAB))
            ALLOCATE (SLAB_TEMP(NUMSLAB))
            ALLOCATE (SLAB_MTEMP(NUMSLAB))
            ALLOCATE (SLAB_MDENS(NUMSLAB))
      ENDIF
       
      ALLOCATE (VPX(NATOMS))
      ALLOCATE (Z_POSITION(NATOMS))

      IF (PPF_INPUT .EQ. 1) THEN
            ALLOCATE (NUM(SLIDE))
            ALLOCATE (DENSITY(SLIDE))
            ALLOCATE (VXPRO(SLIDE))
            ALLOCATE (TEEMP(SLIDE))
            ALLOCATE (RDENSITY(SLIDE))
            ALLOCATE (RVXPRO(SLIDE))
            ALLOCATE (RTEEMP(SLIDE))
            ALLOCATE (PZX(SLIDE))
            ALLOCATE (RPZX(SLIDE))
            ALLOCATE (VVISCOSITY(SLIDE))
            ALLOCATE (RVVISCOSITY(SLIDE))
      ENDIF

      RETURN

    END SUBROUTINE ALLOCATEVAR

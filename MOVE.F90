#ifndef DOXY_SKIP
#include "ibi-preprocess.h"
#endif

!> @file
!> @brief Moves atoms in box according to LF algorithm
!> @details Called from NEW_LOOP.f90

SUBROUTINE MOVE ( )

  USE VAR
  USE OMP_LIB 
  IMPLICIT NONE
  
  INTEGER :: I, TI, A
  REAL(KIND=RKIND) :: LCFAC_ATOM, LCFAC_BEAD, INV_MASS, REG_MASS
  REAL(KIND=RKIND) :: VXI, VYI, VZI
  REAL(KIND=RKIND), DIMENSION(2) :: PT11K, PT22K, PT33K
  REAL(KIND=RKIND), DIMENSION(2) :: PT12K, PT23K, PT13K

  EK = 0.0D0

  PT11K = 0.0D0
  PT22K = 0.0D0
  PT33K = 0.0D0

  PT12K = 0.0D0
  PT23K = 0.0D0
  PT13K = 0.0D0

  IF ((ENSEMBLE == 1).OR.(ENSEMBLE == 2)) THEN !If NVT or NPT
     LCFAC_ATOM = SQRT(1.0D0+FAC*((TIN/TEMP_ATOM)-1.0D0))
     LCFAC_BEAD = SQRT(1.0D0+FAC*((TIN/TEMP_BEAD)-1.0D0)) 
  ELSE
     LCFAC_ATOM = 1.0D0
     LCFAC_BEAD = 1.0D0
  END IF

  DO A=1,NUMATOMS
     I = ATOM(A)
     TI = ITYPE(I)
     REG_MASS = MASS(TI)
     INV_MASS = INVMASS(TI)

     VXI = VX(I)
     VYI = VY(I)
     VZI = VZ(I)

     VX(I) = (VXI + DT * FX(I) * INV_MASS)*LCFAC_ATOM
     VY(I) = (VYI + DT * FY(I) * INV_MASS)*LCFAC_ATOM
     VZ(I) = (VZI + DT * FZ(I) * INV_MASS)*LCFAC_ATOM

     SX(I) = SX(I) + DT * VX(I)
     SY(I) = SY(I) + DT * VY(I)
     SZ(I) = SZ(I) + DT * VZ(I)

     VTX(I) = 0.5D0*(VX(I) + VXI)
     VTY(I) = 0.5D0*(VY(I) + VYI)
     VTZ(I) = 0.5D0*(VZ(I) + VZI)

     VXI = VTX(I)
     VYI = VTY(I)
     VZI = VTZ(I)
  
     PT11K(1) = PT11K(1) + REG_MASS *VXI * VXI
     PT22K(1) = PT22K(1) + REG_MASS *VYI * VYI  
     PT33K(1) = PT33K(1) + REG_MASS *VZI * VZI

     PT12K(1) = PT12K(1) + REG_MASS * VXI * VYI
     PT13K(1) = PT13K(1) + REG_MASS * VXI * VZI
     PT23K(1) = PT23K(1) + REG_MASS * VYI * VZI     
  END DO

  DO A=1,NUMBEADS
     I = BEAD(A)
     TI = ITYPE(I)
     REG_MASS = MASS(TI)
     INV_MASS = INVMASS(TI)

     VXI = VX(I)
     VYI = VY(I)
     VZI = VZ(I)

     VX(I) = (VXI + DT * FX(I) * INV_MASS)*LCFAC_BEAD
     VY(I) = (VYI + DT * FY(I) * INV_MASS)*LCFAC_BEAD
     VZ(I) = (VZI + DT * FZ(I) * INV_MASS)*LCFAC_BEAD

     SX(I) = SX(I) + DT * VX(I)
     SY(I) = SY(I) + DT * VY(I)
     SZ(I) = SZ(I) + DT * VZ(I)

     VTX(I) = 0.5D0*(VX(I) + VXI)
     VTY(I) = 0.5D0*(VY(I) + VYI)
     VTZ(I) = 0.5D0*(VZ(I) + VZI)

     VXI = VTX(I)
     VYI = VTY(I)
     VZI = VTZ(I)
  
     PT11K(2) = PT11K(2) + REG_MASS *VXI * VXI
     PT22K(2) = PT22K(2) + REG_MASS *VYI * VYI  
     PT33K(2) = PT33K(2) + REG_MASS *VZI * VZI

     PT12K(2) = PT12K(2) + REG_MASS * VXI * VYI
     PT13K(2) = PT13K(2) + REG_MASS * VXI * VZI
     PT23K(2) = PT23K(2) + REG_MASS * VYI * VZI  
  END DO

  PT11 = PT11 + SUM(PT11K)
  PT22 = PT22 + SUM(PT22K)
  PT33 = PT33 + SUM(PT33K)

  PT12 = PT12 + SUM(PT12K)
  PT13 = PT13 + SUM(PT13K)
  PT23 = PT23 + SUM(PT23K)

!  EK = 0.5D0 * (SUM(PT11K) + SUM(PT22K) + SUM(PT33K))
  EK(1) = 0.5 * (PT11K(1) + PT22K(1) + PT33K(1))
  EK(2) = 0.5 * (PT11K(2) + PT22K(2) + PT33K(2))

  TEMP = SUM(EK) * MKTEMP
  TEMP_ATOM = EK(1) * MKTEMP_ATOM
  TEMP_BEAD = EK(2) * MKTEMP_BEAD

  RETURN
END SUBROUTINE MOVE

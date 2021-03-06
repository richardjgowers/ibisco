!> @file
!> @brief Resets the net momentum
!> @details Called from NEW_LOOP.f90

SUBROUTINE MOMENTUM ()

  USE VAR
  IMPLICIT NONE

  INTEGER :: I, TI
  REAL*4, DIMENSION(3) :: SUM_XYZ
  REAL*4 :: REG_MASS, INV_MASS

  SUM_XYZ = 0.0

  DO I = 1, NATOMS
     TI = ITYPE(I)
     REG_MASS = MASS(TI)
     SUM_XYZ(:) = SUM_XYZ(:) + VXYZ(:,I) * REG_MASS
  END DO

  SUM_XYZ(:) = SUM_XYZ(:) / REAL(NATOMS)

  !	CHANGE THE VELOCITES SUCH THAT TOTAL MOMENTUM IN EACH DIRECTION IS SETED TO ZERO
  DO I = 1, NATOMS
     TI = ITYPE(I)
     INV_MASS = INVMASS(TI)
     VXYZ(:,I) = VXYZ(:,I) - SUM_XYZ(:) * INV_MASS
  END DO

  RETURN
END SUBROUTINE MOMENTUM

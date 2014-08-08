!> @file
!> @brief Thermostat
!> @details I have no idea if this works

SUBROUTINE SCALEV ()
  USE VAR

  IMPLICIT NONE

  REAL*8 :: FACTOR, T
  INTEGER :: I

  EK(1) = 0.0
  DO  I = 1, NATOMS
     EK(1) = EK(1) + MASS(ITYPE(I))*(VXYZ(1,I)**2.0 + VXYZ(2,I)**2.0 + VXYZ(3,I)**2.0)
  END DO

  EK(1) = 0.5 * EK(1)
  T = EK(1) *  MKTEMP
  FACTOR = TEMP_IN / T

  DO  I = 1, NATOMS
     VXYZ(1,I) = SQRT(FACTOR) * VXYZ(1,I)
     VXYZ(2,I) = SQRT(FACTOR) * VXYZ(2,I)
     VXYZ(3,I) = SQRT(FACTOR) * VXYZ(3,I)
  END DO

  RETURN
END SUBROUTINE SCALEV



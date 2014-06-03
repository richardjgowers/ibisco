!> @file
!> @brief Thermostat
!> @details I have no idea if this works

SUBROUTINE SCALEV ()
  USE VAR

  IMPLICIT NONE

  REAL*8 :: FACTOR, T
  INTEGER ::I

  EK(1) = 0.0D0
  DO  I = 1, NATOMS

     EK(1) = EK(1) + MASS(ITYPE(I))*(VX(I)**2.0D0 + VY(I)**2.0D0 + VZ(I)**2.0D0)

  END DO

  EK(1) = 0.5D0 * EK(1)
  T = EK(1) *  MKTEMP
  FACTOR = TEMP_IN / T

  DO  I = 1, NATOMS
     VX(I) = SQRT(FACTOR) * VX(I)
     VY(I) = SQRT(FACTOR) * VY(I)
     VZ(I) = SQRT(FACTOR) * VZ(I)
  END DO

  RETURN
END SUBROUTINE SCALEV



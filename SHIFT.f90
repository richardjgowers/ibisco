!> @file
!> @brief Moves all atoms to inside the primary unit cell

!
!
!	If SX is n times (where n is even) greater than half side of the box 
!	(in the x-direction in this case) SPX represent the "percentage" of the 
!	position in the	half part of the box where is located the image of the 
!	atom. Therefore the image is located at (SPX * BOX2(1)), because if n is even
!	INT(SPX) is always 0.
!
!	If SX is n times (where n is odd) greater than half side of the box 
!	(in the x-direction in this case) (SPX - 2.0*INT(SPX)) represent the 
!	"percentage" of the position in the half part of the box where is located 
!	the image of the atom. If n is odd [SPX > 1]. 
!
!
!
!
SUBROUTINE SHIFT()

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER :: I
  REAL*4, DIMENSION(3) ::SPXYZ

  ! INT(funzione parte intera) is an intrinsic function that reads a real, cuts the decimal 
  ! part of the number leaving only its integer part (the final format is real)

  !$OMP PARALLEL DO SCHEDULE(STATIC,1) DEFAULT(NONE) &
  !$OMP& SHARED(NATOMS,SXYZ,BOX2,RXYZ)&
  !$OMP& PRIVATE(I,SPXYZ)
  DO I = 1, NATOMS
     IF (ANY((SXYZ(:,I) > BOX2(:)) .OR. (SXYZ(:,I) < -BOX2(:)))) THEN
        RXYZ(:,I) = SXYZ(:,I) / BOX2(:)
        SPXYZ(:) = RXYZ(:,I) - 2.0 * INT(RXYZ(:,I) / 2.0)
        RXYZ(:,I) = (SPXYZ(:) - 2.0 * INT(SPXYZ(:))) * BOX2(:)
     ELSE
        RXYZ(:,I) = SXYZ(:,I)
     END IF
  END DO
  !$OMP END PARALLEL DO

  RETURN
END SUBROUTINE SHIFT

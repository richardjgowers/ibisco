!> @file
!> @brief Moves all atoms to inside the primary unit cell

!
!
!	If SX is n times (where n is even) greater than half side of the box 
!	(in the x-direction in this case) SPX represent the "percentage" of the 
!	position in the	half part of the box where is located the image of the 
!	atom. Therefore the image is located at (SPX * BOXX2), because if n is even
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
  !      *******************************************************************

  ! INT(funzione parte intera) is an intrinsic function that reads a real, cuts the decimal 
  ! part of the number leaving only its integer part (the final format is real)

  !$OMP PARALLEL DO SCHEDULE(STATIC,1) DEFAULT(NONE) &
  !$OMP& SHARED(NATOMS,SXYZ,BOXX2,BOXY2,BOXZ2,RXYZ)&
  !$OMP& PRIVATE(I,SPXYZ)
  DO I = 1, NATOMS
     IF ((SXYZ(1,I)> BOXX2).OR.(SXYZ(1,I) < -BOXX2).OR.(SXYZ(2,I)> BOXY2) &
          .OR.(SXYZ(2,I)< -BOXY2).OR.(SXYZ(3,I)> BOXZ2).OR.(SXYZ(3,I)<-BOXZ2)) THEN

        RXYZ(1,I) = SXYZ(1,I)/BOXX2
        RXYZ(2,I) = SXYZ(2,I)/BOXY2
        RXYZ(3,I) = SXYZ(3,I)/BOXZ2

        SPXYZ(1) = RXYZ(1,I) - 2.0*INT(RXYZ(1,I)/2.0)
        SPXYZ(2) = RXYZ(2,I) - 2.0*INT(RXYZ(2,I)/2.0)
        SPXYZ(3) = RXYZ(3,I) - 2.0*INT(RXYZ(3,I)/2.0)

        RXYZ(1,I) = (SPXYZ(1) - 2.0*INT(SPXYZ(1)))*BOXX2
        RXYZ(2,I) = (SPXYZ(2) - 2.0*INT(SPXYZ(2)))*BOXY2
        RXYZ(3,I) = (SPXYZ(3) - 2.0*INT(SPXYZ(3)))*BOXZ2
     ELSE

        RXYZ(1,I) = SXYZ(1,I)
        RXYZ(2,I) = SXYZ(2,I)
        RXYZ(3,I) = SXYZ(3,I)
     END IF
  END DO
  !$OMP END PARALLEL DO

  RETURN
END SUBROUTINE SHIFT

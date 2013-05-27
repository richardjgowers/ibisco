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
  REAL*8 ::SPX, SPY, SPZ
  real(kind=rkind) :: t1,t2,tick
  !       *******************************************************************

  ! INT(funzione parte intera) is an intrinsic function that reads a real, cuts the decimal 
  ! part of the number leaving only its integer part (the final format is real)

  ! t1 = omp_get_wtime()
  !$OMP PARALLEL DO SCHEDULE(STATIC,1) NUM_THREADS(7) DEFAULT(NONE) &
  !$OMP& SHARED(SX,SY,SZ,BOXX2,BOXY2,BOXZ2,RX,RY,RZ)&
  !$OMP& private(SPX,SPY,SPZ)
  DO I = 1, NATOMS
     IF ((SX(I)> BOXX2).OR.(SX(I) < -BOXX2).OR.(SY(I)> BOXY2) &
          .OR.(SY(I)< -BOXY2).OR.(SZ(I)> BOXZ2).OR.(SZ(I)<-BOXZ2)) THEN

        RX(I) = SX(I)/BOXX2
        RY(I) = SY(I)/BOXY2
        RZ(I) = SZ(I)/BOXZ2

        SPX = RX(I) - 2.0*INT(RX(I)/2.0)
        SPY = RY(I) - 2.0*INT(RY(I)/2.0)
        SPZ = RZ(I) - 2.0*INT(RZ(I)/2.0)

        RX(I) = (SPX - 2.0*INT(SPX))*BOXX2
        RY(I) = (SPY - 2.0*INT(SPY))*BOXY2
        RZ(I) = (SPZ - 2.0*INT(SPZ))*BOXZ2
     ELSE

        RX(I) = SX(I)
        RY(I) = SY(I)
        RZ(I) = SZ(I)
     END IF
  END DO
  !$OMP END PARALLEL DO

  !  t2=omp_get_wtime()
  !  tick=omp_get_wtick()

  !        write(4000,*) 'Time elapsed ',t2 - t1,' Precision: ',tick
  RETURN
END SUBROUTINE SHIFT

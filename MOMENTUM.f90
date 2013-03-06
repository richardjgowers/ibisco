
      SUBROUTINE MOMENTUM ()

      USE VAR
      IMPLICIT NONE

      INTEGER :: I
      REAL*8 :: SUMX, SUMY, SUMZ
!       ******************************************************************
!	***********	RESET THE NET MOMENTUM TO ZERO    ****************

!	CALCULATE THE TOTAL MOMENTUM FOR EACH DIRECTION

      SUMX = 0.0D0
      SUMY = 0.0D0
      SUMZ = 0.0D0

        DO 205 I = 1, NATOMS
            CM = BEADMASS(I)
            SUMX = SUMX + VX(I) * CM
            SUMY = SUMY + VY(I) * CM
            SUMZ = SUMZ + VZ(I) * CM
205     CONTINUE

      SUMX = SUMX/REAL (NATOMS)
      SUMY = SUMY/REAL (NATOMS)
      SUMZ = SUMZ/REAL (NATOMS)

!	CHANGE THE VELOCITES SUCH THAT TOTAL MOMENTUM IN EACH DIRECTION IS SETED TO ZERO
        DO 300 I = 1, NATOMS

            CM = BEADMASS(I)
            VX(I) = VX(I) - SUMX/CM
            VY(I) = VY(I) - SUMY/CM
            VZ(I) = VZ(I) - SUMZ/CM

300     CONTINUE

      RETURN
      END

!	*********************************************************************************************

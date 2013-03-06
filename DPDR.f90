	SUBROUTINE DPDR ( )

	USE VAR
	IMPLICIT	NONE
        INTEGER       I
!       *******************************************************************

        DO 100 I = 1, NATOMS

           SX(I) = SX(I) + DT * VX(I)
           SY(I) = SY(I) + DT * VY(I)
           SZ(I) = SZ(I) + DT * VZ(I) 

100     CONTINUE

	
        RETURN
        END
!	*********************************************************************************************

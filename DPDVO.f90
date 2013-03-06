	SUBROUTINE DPDVO (VVONST)

	USE VAR
	IMPLICIT	NONE
        INTEGER       I
        REAL*8        VVONST
!       *******************************************************************

        DO 100 I = 1, NATOMS

	   CM = BEADMASS(I)

	   VOX(I) = VX(I) + DT * FX(I) *VVONST / CM 
           VOY(I) = VY(I) + DT * FY(I) *VVONST / CM 
           VOZ(I) = VZ(I) + DT * FZ(I) *VVONST / CM 

100     CONTINUE

        RETURN
        END
!	*********************************************************************************************

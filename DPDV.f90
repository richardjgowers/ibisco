	SUBROUTINE DPDV (VVCONS)

	USE VAR
	IMPLICIT	NONE
        INTEGER       I
        REAL*8        VVCONS
!       *******************************************************************

        DO 100 I = 1, NATOMS

	   CM = BEADMASS(I)

	   VX(I) = VX(I) + DT * FX(I) *VVCONS / CM 
           VY(I) = VY(I) + DT * FY(I) *VVCONS / CM 
           VZ(I) = VZ(I) + DT * FZ(I) *VVCONS / CM 
           
           VTX(I) = VX(I)
           VTY(I) = VY(I)
           VTZ(I) = VZ(I)

100     CONTINUE
        RETURN
        END
!	*********************************************************************************************

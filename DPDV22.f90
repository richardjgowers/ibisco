	SUBROUTINE DPDV22 (VVCOST)

	USE VAR
	IMPLICIT	NONE
        INTEGER       I
        REAL*8        VVCOST,RZI
        real*8       accelx,accely,accelz
!       *******************************************************************

        DO 100 I = 1, NATOMS

	   CM = BEADMASS(I)
           RZI = SZ(I)
           RZI = RZI - BOXZ*DNINT(RZI*BOXZINV)
           
           IF (RZI.GE.0.D0) THEN
           ACCELX = FX(I)/CM + EXTER_FX 
           ELSEIF (RZI.LT.0.D0) THEN
           ACCELX = FX(I)/CM - EXTER_FX
           ENDIF
           ACCELY = FY(I)/CM
           ACCELZ = FZ(I)/CM

	   VX(I) = VX(I) + DT *ACCELX *VVCOST
           VY(I) = VY(I) + DT *ACCELY *VVCOST
           VZ(I) = VZ(I) + DT *ACCELZ *VVCOST
           
           VTX(I) = VX(I)
           VTY(I) = VY(I)
           VTZ(I) = VZ(I)

100     CONTINUE

        RETURN
        END
!	*********************************************************************************************

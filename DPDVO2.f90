	SUBROUTINE DPDVO2 (VVCNST)

	USE VAR
	IMPLICIT	NONE
        INTEGER       I
        REAL*8        VVCNST,RZI
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

	   VOX(I) = VX(I) + DT *ACCELX *VVCNST
           VOY(I) = VY(I) + DT *ACCELY *VVCNST
           VOZ(I) = VZ(I) + DT *ACCELZ *VVCNST

100     CONTINUE

        RETURN
        END
!	*********************************************************************************************

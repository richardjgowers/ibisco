
      SUBROUTINE SCALEBP (TM)

      USE VAR
      IMPLICIT NONE
      INTEGER :: I, TM
      REAL*8 :: P, AVP
      REAL*8 :: MIUX, MIUY, MIUZ
!      REAL*8 :: DBX, DBY, DBZ, NBOXX, NBOXY, NBOXZ, MIUX, MIUY, MIUZ
!      REAL*8 :: PSCALEX, PSCALEY, PSCALEZ, P, AVP
!       ******************************************************************************
!	***** CALCULATE THE NEW VOLUME AND SCALE THE POSITIONS IN NPT ENSEMBLE	******

!	CALCULATE AVERAGE
      IAVP = IAVP + 1
            IF (IAVP .GT. LIMAVP) IAVP = IAVP - LIMAVP

      SP(IAVP) = (PT11 + PT22 + PT33)/3.0d0

      IF (MOD(TM, LIMAVP) == 0) THEN
            AVP = 0.0D0
            DO I = 1, LIMAVP
                  AVP = AVP + SP(I)
            END DO
 
      P = AVP / LIMAVP

      MIUX = (1.0D0 + (P-PRESSURE)*PFAC)**(1.0d0/3.0d0)
      MIUY = (1.0D0 + (P-PRESSURE)*PFAC)**(1.0d0/3.0d0)
      MIUZ = (1.0D0 + (P-PRESSURE)*PFAC)**(1.0d0/3.0d0)
	
      BOXX = BOXX * MIUX
      BOXY = BOXY * MIUY
      BOXZ = BOXZ * MIUZ

      DO I = 1, NATOMS

            SX(I) = SX(I)*MIUX
            SY(I) = SY(I)*MIUY
            SZ(I) = SZ(I)*MIUZ

      END DO

      BOXXINV = 1.0D0 / BOXX
      BOXYINV = 1.0D0 / BOXY
      BOXZINV = 1.0D0 / BOXZ

      BOXX2 = BOXX / 2.0D0
      BOXY2 = BOXY / 2.0D0
      BOXZ2 = BOXZ / 2.0D0

!	P = (PT11 + PT22 + PT33)/3.0d0
	
!	DBX = 2.0D0*(P-PRESSURE)*PFAC
!	DBY = 2.0D0*(P-PRESSURE)*PFAC
!	DBZ = 2.0D0*(P-PRESSURE)*PFAC
	
!	NBOXX = BOXX + DBX
!	NBOXY = BOXY + DBY
!	NBOXZ = BOXZ + DBZ

!	PSCALEX = 1.0D0 + DBX/BOXX
!	PSCALEY = 1.0D0 + DBY/BOXY
!	PSCALEZ = 1.0D0 + DBZ/BOXZ

!	DO I = 1, NATOMS

!		SX(I) = SX(I)*PSCALEX
!		SY(I) = SY(I)*PSCALEY
!		SZ(I) = SZ(I)*PSCALEZ
!	END DO

!	BOXX = NBOXX
!	BOXY = NBOXY
!	BOXZ = NBOXZ
!	BOXXINV = 1.0D0 / BOXX
!	BOXYINV = 1.0D0 / BOXY
!	BOXZINV = 1.0D0 / BOXZ

!	BOXX2 = BOXX / 2.0D0
!	BOXY2 = BOXY / 2.0D0
!	BOXZ2 = BOXZ / 2.0D0

      END IF 

      RETURN
      END
!	*********************************************************************************************

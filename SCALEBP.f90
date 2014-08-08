!> @file
!> @brief Berensden barostat
!> @details Called from NEW_LOOP()
!> @details Rescales positions and box size according to system pressure

SUBROUTINE SCALEBP (TM)

  USE VAR
  IMPLICIT NONE
  INTEGER :: I, TM
  REAL*4 :: P, AVP
  REAL*4, DIMENSION(3) :: MIU

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

     MIU(1) = (1.0 + (P - PRESSURE) * PFAC) ** (1.0 / 3.0)
     MIU(2) = (1.0 + (P - PRESSURE) * PFAC) ** (1.0 / 3.0)
     MIU(3) = (1.0 + (P - PRESSURE) * PFAC) ** (1.0 / 3.0)

     BOXX = BOXX * MIU(1)
     BOXY = BOXY * MIU(2)
     BOXZ = BOXZ * MIU(3)

     DO I = 1, NATOMS
        SXYZ(1,I) = SXYZ(1,I) * MIU(1)
        SXYZ(2,I) = SXYZ(2,I) * MIU(2)
        SXYZ(3,I) = SXYZ(3,I) * MIU(3)
     END DO

     BOXXINV = 1.0 / BOXX
     BOXYINV = 1.0 / BOXY
     BOXZINV = 1.0 / BOXZ

     BOXX2 = BOXX / 2.0
     BOXY2 = BOXY / 2.0
     BOXZ2 = BOXZ / 2.0
  END IF

  RETURN
END SUBROUTINE SCALEBP

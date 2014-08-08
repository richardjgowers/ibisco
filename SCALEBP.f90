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

     MIU(:) = (1.0 + (P - PRESSURE) * PFAC) ** (1.0 / 3.0)

     BOX(:) = BOX(:) * MIU(:)
     BOXINV(:) = 1.0 / BOX(:)
     BOX2(:) = BOX(:) / 2.0

     DO I = 1, NATOMS
        SXYZ(:,I) = SXYZ(:,I) * MIU(:)
     END DO
     ! SXYZ(:,:) = SXYZ(:,:) * MIU(:)
  END IF

  RETURN
END SUBROUTINE SCALEBP

      SUBROUTINE MOVE2 ()

      USE VAR
      IMPLICIT NONE

      INTEGER       I, IZI
      REAL*8        RZI, LCFAC
      REAL*8            VXI, VYI, VZI
      REAL*8            PT11K, PT22K, PT33K
      REAL*8            PT12K, PT23K, PT13K
      REAL*8        ACCELX,ACCELY,ACCELZ


!      INTEGER       I, J, IZI
!      REAL*8        RXI, RYI, RZI, LCFAC
!      REAL*8            VXI, VYI, VZI, EKNEW
!       *******************************************************************
      EK = 0.0D0

      PT11K = 0.0D0
      PT22K = 0.0D0
      PT33K = 0.0D0

      PT12K = 0.0D0
      PT23K = 0.0D0
      PT13K = 0.0D0

      IF ((ENSEMBLE == 1).OR.(ENSEMBLE == 2)) THEN
            LCFAC = SQRT(1.0D0+FAC*((TIN/TEMP)-1.0D0))
            ELSE
            LCFAC = 1.0D0
      END IF

        DO 100 I = 1, NATOMS
           
            RZI = SZ(I)
            RZI = RZI - BOXZ*DNINT(RZI*BOXZINV)
           
            CM = BEADMASS(I)
            
            IF (RZI .GE. 0.d0) THEN
                  ACCELX = FX(I) / CM + EXTER_FX 
            ELSE IF (RZI .LT. 0.d0) THEN
                  ACCELX = FX(I) / CM - EXTER_FX
            END IF

            ACCELY = FY(I) / CM
            ACCELZ = FZ(I) / CM


            VXI = VX(I)
            VYI = VY(I)
            VZI = VZ(I)

            VX(I) = (VXI + DT * ACCELX)*LCFAC
            VY(I) = (VYI + DT * ACCELY)*LCFAC
            VZ(I) = (VZI + DT * ACCELZ)*LCFAC

            SX(I) = SX(I) + DT * VX(I)
            SY(I) = SY(I) + DT * VY(I)
            SZ(I) = SZ(I) + DT * VZ(I)

            VTX(I) = 0.5D0*(VX(I) + VXI)           
            VTY(I) = 0.5D0*(VY(I) + VYI)
            VTZ(I) = 0.5D0*(VZ(I) + VZI)

100     CONTINUE

        CALL PROFILEPPF()

        DO 200 I = 1, NATOMS

            CM = BEADMASS(I)
            VXI = VPX(I)
            PT11K = PT11K + CM *VXI ** 2.0
            PT22K = PT22K + CM *VTY(I) ** 2.0
            PT33K = PT33K + CM *VTZ(I) ** 2.0

            PT12K = PT12K + CM * VXI * VTY(I)
            PT13K = PT13K + CM * VXI * VTZ(I)
            PT23K = PT23K + CM * VTY(I) * VTZ(I)
          
            RZI = SZ(I)
            RZI = RZI - BOXZ*DNINT(RZI*BOXZINV)
            IZI = INT((RZI*BOXZINV+0.5D0)*DBLE(SLIDE))
            IZI = 1 + MOD((IZI+SLIDE),SLIDE)
            PZX(IZI) = PZX(IZI) + CM*VTZ(I) * VXI           

200     CONTINUE

      PT11 = PT11 + PT11K 
      PT22 = PT22 + PT22K
      PT33 = PT33 + PT33K

      PT12 = PT12 + PT12K 
      PT13 = PT13 + PT13K
      PT23 = PT23 + PT23K

      EK = 0.5D0 * (PT11K + PT22K + PT33K)

      TEMP = EK * MKTEMP
      
        RETURN
        END
!      *********************************************************************************************

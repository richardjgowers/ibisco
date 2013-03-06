      SUBROUTINE MOVE ( )

      USE VAR
      IMPLICIT NONE

      INTEGER       I
      REAL*8        LCFAC
      REAL*8            VXI, VYI, VZI
      REAL*8            PT11K, PT22K, PT33K
      REAL*8            PT12K, PT23K, PT13K

!      REAL*8        RXI, RYI, RZI, EKNEW
!      INTEGER       I, J
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
            CM = BEADMASS(I)

            VXI = VX(I)
            VYI = VY(I)
            VZI = VZ(I)

            VX(I) = (VXI + DT * FX(I) / CM)*LCFAC
            VY(I) = (VYI + DT * FY(I) / CM)*LCFAC
            VZ(I) = (VZI + DT * FZ(I) / CM)*LCFAC

            SX(I) = SX(I) + DT * VX(I)
            SY(I) = SY(I) + DT * VY(I)
            SZ(I) = SZ(I) + DT * VZ(I)

            VTX(I) = 0.5D0*(VX(I) + VXI)
            VTY(I) = 0.5D0*(VY(I) + VYI)
            VTZ(I) = 0.5D0*(VZ(I) + VZI)

!100     CONTINUE

!        CALL PROFILERNEMD()

!        DO 200 I = 1, NATOMS
!           CM = BEADMASS(I)
           VXI = VTX(I)

           PT11K = PT11K + CM *VXI ** 2.0
           PT22K = PT22K + CM *VTY(I) ** 2.0
           PT33K = PT33K + CM *VTZ(I) ** 2.0

           PT12K = PT12K + CM * VXI * VTY(I)
           PT13K = PT13K + CM * VXI * VTZ(I)
           PT23K = PT23K + CM * VTY(I) * VTZ(I)

           FX(I) = 0.0D0
           FY(I) = 0.0D0
           FZ(I) = 0.0D0

100     CONTINUE

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

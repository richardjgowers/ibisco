      SUBROUTINE MOVE ( )

      USE VAR
      USE OMP_LIB 
      IMPLICIT NONE

      INTEGER       I
      REAL*8        LCFAC
      REAL*8            VXI, VYI, VZI
      REAL*8            PT11K, PT22K, PT33K
      REAL*8            PT12K, PT23K, PT13K
      real(kind=rkind) :: t1,t2,tick

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
     
t1 = omp_get_wtime()
!$OMP PARALLEL DO  DEFAULT(SHARED) SCHEDULE(GUIDED)&
!$OMP& REDUCTION(+: PT11K,PT22K,PT33K,PT12K,PT13K,PT23K)&
!$OMP& private(CM,VXI,VZI,VYI)


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
           VZI = VTZ(I)
           VYI = VTY(I)  
           PT11K = PT11K + CM *VXI * VXI
           PT22K = PT22K + CM *VYI * VYI  
           PT33K = PT33K + CM *VZI * VZI

           PT12K = PT12K + CM * VXI * VYI
           PT13K = PT13K + CM * VXI * VZI
           PT23K = PT23K + CM * VYI * VZI

          FX(I) = 0.0D0
           FY(I) = 0.0D0
           FZ(I) = 0.0D0

100     CONTINUE
!$OMP END PARALLEL DO
  t2=omp_get_wtime()
  tick=omp_get_wtick()

      write(4001,*) 'Time elapsed ',t2 - t1,' Precision: ',tick

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

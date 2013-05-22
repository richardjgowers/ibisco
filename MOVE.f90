SUBROUTINE MOVE ( )

  USE VAR
  USE OMP_LIB 
  IMPLICIT NONE
  
  INTEGER :: I
  REAL(KIND=RKIND) :: LCFAC, INV_MASS, REG_MASS
  REAL(KIND=RKIND) :: VXI, VYI, VZI
  REAL(KIND=RKIND) :: PT11K, PT22K, PT33K
  REAL(KIND=RKIND) :: PT12K, PT23K, PT13K
!  real(kind=rkind) :: t1,t2,tick

  EK = 0.0D0

  PT11K = 0.0D0
  PT22K = 0.0D0
  PT33K = 0.0D0

  PT12K = 0.0D0
  PT23K = 0.0D0
  PT13K = 0.0D0

  IF ((ENSEMBLE == 1).OR.(ENSEMBLE == 2)) THEN !If NVT or NPT
     LCFAC = SQRT(1.0D0+FAC*((TIN/TEMP)-1.0D0))
  ELSE
     LCFAC = 1.0D0
  END IF

  !t1 = omp_get_wtime()
  !$OMP PARALLEL DO  DEFAULT(SHARED) SCHEDULE(GUIDED) NUM_THREADS(6) &
  !$OMP& REDUCTION(+: PT11K,PT22K,PT33K,PT12K,PT13K,PT23K)&
  !$OMP& private(CM,VXI,VZI,VYI)
  DO I = 1, NATOMS

     REG_MASS = MASS(ITYPE(I))
     INV_MASS = INVMASS(ITYPE(I))

     VXI = VX(I)
     VYI = VY(I)
     VZI = VZ(I)

     VX(I) = (VXI + DT * FX(I) * INV_MASS)*LCFAC
     VY(I) = (VYI + DT * FY(I) * INV_MASS)*LCFAC
     VZ(I) = (VZI + DT * FZ(I) * INV_MASS)*LCFAC

     SX(I) = SX(I) + DT * VX(I)
     SY(I) = SY(I) + DT * VY(I)
     SZ(I) = SZ(I) + DT * VZ(I)

     VTX(I) = 0.5D0*(VX(I) + VXI)
     VTY(I) = 0.5D0*(VY(I) + VYI)
     VTZ(I) = 0.5D0*(VZ(I) + VZI)

     VXI = VTX(I)
     VZI = VTZ(I)
     VYI = VTY(I)
  
     PT11K = PT11K + REG_MASS *VXI * VXI
     PT22K = PT22K + REG_MASS *VYI * VYI  
     PT33K = PT33K + REG_MASS *VZI * VZI

     PT12K = PT12K + REG_MASS * VXI * VYI
     PT13K = PT13K + REG_MASS * VXI * VZI
     PT23K = PT23K + REG_MASS * VYI * VZI
  END DO
  !$OMP END PARALLEL DO
  !  t2=omp_get_wtime()
  !  tick=omp_get_wtick()

  !      write(4001,*) 'Time elapsed ',t2 - t1,' Precision: ',tick

  PT11 = PT11 + PT11K
  PT22 = PT22 + PT22K 
  PT33 = PT33 + PT33K 

  PT12 = PT12 + PT12K
  PT13 = PT13 + PT13K
  PT23 = PT23 + PT23K

  EK = 0.5D0 * (PT11K + PT22K + PT33K)

  TEMP = EK * MKTEMP

  RETURN
END SUBROUTINE MOVE

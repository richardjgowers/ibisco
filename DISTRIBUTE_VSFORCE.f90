SUBROUTINE DISTRIBUTE_VSFORCE()

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER :: I, J, TI, TJ, A, VS_POS

  !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,1) &
  !$OMP& SHARED(NVIRTA,VITYPE,NATOMS,ITYPE)&
  !$OMP& SHARED(VIRT_NUMATOMS,VIRT_ATM_IND,VIRT_MASSCOEFF)&
  !$OMP& SHARED(FX,FY,FZ)&
  !$OMP& PRIVATE(I,TI,VS_POS,A,J,TJ)
  !FX FY and FZ do not need to be reduction variables because only one virtual site will ever address them
  DO I=1,NVIRTA
     TI = VITYPE(I)
     VS_POS = NATOMS + I
     DO A=1,VIRT_NUMATOMS(TI)
        J = VIRT_ATM_IND(I,A)
        TJ = ITYPE(J)
        FX(J) = FX(J) + FX(VS_POS)*VIRT_MASSCOEFF(TI,TJ)
        FY(J) = FY(J) + FY(VS_POS)*VIRT_MASSCOEFF(TI,TJ)
        FZ(J) = FZ(J) + FZ(VS_POS)*VIRT_MASSCOEFF(TI,TJ)
     END DO
     !Reset force on virtual site to 0 once distributed
     FX(VS_POS) = 0.0
     FY(VS_POS) = 0.0
     FZ(VS_POS) = 0.0
  END DO
  !$OMP END PARALLEL DO

  RETURN
END SUBROUTINE DISTRIBUTE_VSFORCE

!> @file
!> @brief Moves forces from VS to atoms underneath
!!
!> @details This is used in hybrid simulations to transfer forces from the virtual sites onto the
!! atoms.  This is done on a mass weighted basis.  This subroutine is called twice, first after 
!! nonbonded forces have been calculated, and then after bonded force has been calculated. This is
!! done so that the nonbonded force can be properly accumulated for pressure calculations.
!!
!> @note Parallelised using OpenMP
!> @author Richard J Gowers

SUBROUTINE DISTRIBUTE_VSFORCE()

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER :: I, J, TI, TJ, A, VS_POS

  !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,1) &
  !$OMP& SHARED(NVIRTA,VITYPE,NATOMS,ITYPE)&
  !$OMP& SHARED(VIRT_NUMATOMS,VIRT_ATM_IND,VIRT_MASSCOEFF)&
  !$OMP& SHARED(FXYZ)&
  !$OMP& PRIVATE(I,TI,VS_POS,A,J,TJ)
  !FX FY and FZ do not need to be reduction variables because only one virtual site will ever address them
  DO I=1,NVIRTA
     TI = VITYPE(I)
     VS_POS = NATOMS + I
     DO A=1,VIRT_NUMATOMS(TI)
        J = VIRT_ATM_IND(I,A)
        TJ = ITYPE(J)
        FXYZ(1,J) = FXYZ(1,J) + FXYZ(1,VS_POS) * VIRT_MASSCOEFF(TI,TJ)
        FXYZ(2,J) = FXYZ(2,J) + FXYZ(2,VS_POS) * VIRT_MASSCOEFF(TI,TJ)
        FXYZ(3,J) = FXYZ(3,J) + FXYZ(3,VS_POS) * VIRT_MASSCOEFF(TI,TJ)
     END DO
     !Reset force on virtual site to 0 once distributed
     FXYZ(1,VS_POS) = 0.0
     FXYZ(2,VS_POS) = 0.0
     FXYZ(3,VS_POS) = 0.0
  END DO
  !$OMP END PARALLEL DO

  RETURN
END SUBROUTINE DISTRIBUTE_VSFORCE

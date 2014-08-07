!> @file
!> @brief The multiple time step (MTS) subroutines
!> @details MTS allows hybrid scale simulations to be ran faster.  Normally, hybrid scale simulations
!!          are ran using an atomistic time step, which is around an order of magnitude smaller than
!!          the native coarse grained counterpart.  Using MTS allows the nonbonded coarse grained 
!!          force to be approximated in some time steps, reducing the computational load.
!!
!!          For further details, see \cite MTS
!!          
!> @author Richard J Gowers

!> @brief The alternate molecular dynamics loop for MTS simulations
SUBROUTINE MTS_LOOP
  ! Could merge this with main loop, and have if statement around MTS_FORCE call to differentiate?
  USE VAR
  USE MTS

  IMPLICIT NONE

  ALLOCATE(MTS_FX(3, NCOARSE), MTS_FY(3, NCOARSE), MTS_FZ(3, NCOARSE))
  ALLOCATE(MTS_PT11(3), MTS_PT22(3), MTS_PT33(3), MTS_PT12(3), MTS_PT13(3), MTS_PT23(3))
  ALLOCATE(MTS_V_NB(3, 3))

  DO STEP = 1, NSTEP

     CALL SHIFT()

     CALL VIRTUAL_DEF()

     IF(MOD(STEP, NUPDATE) .eq. 0) THEN
        CALL UPDATE_NEIGHBOURLIST()
     END IF

     CALL MTS_FORCE(STEP)

     CALL MOVE()

     IF(MOD(STEP, HALT_DRIFT) .eq. 0) THEN
        CALL MOMENTUM()
     END IF

     IF (ENSEMBLE .eq. 2) THEN
        CALL SCALEBP(STEP)
     END IF

     CALL AVERAGE(STEP)

     IF (MOD(STEP, NTRJ) .eq. 0) THEN
        CALL WRITETRJ(STEP)
        CALL OUTPUT(STEP)
     END IF

  END DO

  DEALLOCATE(MTS_FX, MTS_FY, MTS_FZ)
  DEALLOCATE(MTS_PT11, MTS_PT22, MTS_PT33, MTS_PT12, MTS_PT13, MTS_PT23)
  DEALLOCATE(MTS_V_NB)

  RETURN
END SUBROUTINE MTS_LOOP

!> @brief The force loop for MTS simulation
!> @details In this force step, the current step number is used to determine how the coarse grained
!!          nonbonded forces are calculated. The first three are done explicitly, then NMTS (an input 
!!          setting) steps are done using an approximation.
!!
!!          This approximation is a Taylor Expansion
!!
!!          The coarse grained contributions to nonbonded energy (V_NB) and pressure tensors (PTxx)
!!          are also estimated in the same fashion.
SUBROUTINE MTS_FORCE(STEPNO)

  USE VAR
  USE MTS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: STEPNO !< The current step number
  INTEGER :: I, MTS_MOD

  ! Reset force
  FX = 0.0D0
  FY = 0.0D0
  FZ = 0.0D0

  FXNB = 0.0D0
  FYNB = 0.0D0
  FZNB = 0.0D0

  PT11 = 0.0D0
  PT22 = 0.0D0
  PT33 = 0.0D0
  PT12 = 0.0D0
  PT13 = 0.0D0
  PT23 = 0.0D0

  !Reset potential energy measures
  V_NB = 0.0
  V_BOND = 0.0
  V_ANGLE = 0.0
  V_TORSION = 0.0
  V_OOP = 0.0

  MTS_MOD = MOD(STEPNO, (3 + NMTS))
  SELECT CASE(MTS_MOD) 
     CASE(1, 2, 3) ! Explicit steps
        CALL NONBONDED_FORCE(NCOARSE, BEAD, MAXNAB_BEAD, LIST_BEAD, RCUT_BEAD, RCUTSQ_BEAD)

        CALL MTS_SAVEFORCE(MTS_MOD, BEAD, FX, FY, FZ, NITEMS, MTS_FX, MTS_FY, MTS_FZ, NCOARSE)
        CALL MTS_SAVEAUXS(MTS_MOD)

     CASE DEFAULT ! Use approximation
        CALL MTS_APPROX(MTS_MOD, BEAD, FX, FY, FZ, NITEMS, MTS_FX, MTS_FY, MTS_FZ, NCOARSE)
        CALL MTS_LOADAUXS(MTS_MOD)
  END SELECT
  CALL DISTRIBUTE_VSFORCE()

  ! Call atoms second, so that V_NB and PTxx only have CG contributions
  ! Order shouldnt matter anyway...
  CALL NONBONDED_FORCE(NUMATOMS, ATOM, MAXNAB_ATOM, LIST_ATOM, RCUT_ATOM, RCUTSQ_ATOM)

  DO I = 1, NATOMS
     FXNB(I) = FX(I)
     FYNB(I) = FY(I)
     FZNB(I) = FZ(I)
  END DO

  CALL BONDED_FORCE()

  CALL DISTRIBUTE_VSFORCE()
  
  DO I=1,NATOMS
     PT11 = PT11 + (FX(I) - FXNB(I))*SX(I)
     PT22 = PT22 + (FY(I) - FYNB(I))*SY(I)
     PT33 = PT33 + (FZ(I) - FZNB(I))*SZ(I)
     PT12 = PT12 + (FY(I) - FYNB(I))*SX(I)
     PT13 = PT13 + (FZ(I) - FZNB(I))*SX(I)
     PT23 = PT23 + (FZ(I) - FZNB(I))*SY(I)
  END DO

  RETURN
END SUBROUTINE MTS_FORCE

!> @brief Records the coarse grained force from an explicit MTS step.
SUBROUTINE MTS_SAVEFORCE(I, BEAD, FX, FY, FZ, NITEMS, MTS_FX, MTS_FY, MTS_FZ, NCOARSE)

  USE VAR, ONLY : RKIND

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I !< Mod of current MTS step
  INTEGER, INTENT(IN) :: NITEMS, NCOARSE !< Total number of particles (real and virtual)
  REAL*4, DIMENSION(NITEMS), INTENT(IN) :: FX, FY, FZ
  REAL*4, DIMENSION(3, NCOARSE), INTENT(INOUT) :: MTS_FX, MTS_FY, MTS_FZ
  INTEGER, DIMENSION(NCOARSE), INTENT(IN) :: BEAD
  
  INTEGER :: A, ATOM_ID

  DO A = 1, NCOARSE
     ATOM_ID = BEAD(A)
     MTS_FX(I, A) = FX(ATOM_ID)
     MTS_FY(I, A) = FY(ATOM_ID)
     MTS_FZ(I, A) = FZ(ATOM_ID)
  END DO
  
  RETURN
  
END SUBROUTINE MTS_SAVEFORCE

!> @brief Generates approximated coarse grained forces for MTS steps.
SUBROUTINE MTS_APPROX(I, BEAD, FX, FY, FZ, NITEMS, MTS_FX, MTS_FY, MTS_FZ, NCOARSE)

  USE VAR, ONLY : RKIND

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I ! Mod of current MTS step
  INTEGER, INTENT(IN) :: NITEMS, NCOARSE !< Total number of particles (real and virtual)
  REAL*4, DIMENSION(NITEMS), INTENT(INOUT) :: FX, FY, FZ
  REAL*4, DIMENSION(3, NCOARSE), INTENT(IN) :: MTS_FX, MTS_FY, MTS_FZ
  INTEGER, DIMENSION(NCOARSE), INTENT(IN) :: BEAD

  INTEGER :: A, ATOM_ID

  DO A = 1, NCOARSE
     ATOM_ID = BEAD(A)
     ! Forwards approx
     FX(ATOM_ID) = MTS_FX(3, A) + (I-3) * 0.5 * (MTS_FX(3, A) - MTS_FX(1, A)) 
     FY(ATOM_ID) = MTS_FY(3, A) + (I-3) * 0.5 * (MTS_FY(3, A) - MTS_FY(1, A)) 
     FZ(ATOM_ID) = MTS_FZ(3, A) + (I-3) * 0.5 * (MTS_FZ(3, A) - MTS_FZ(1, A)) 
  END DO

  RETURN

END SUBROUTINE MTS_APPROX

!> @brief Records the nonbonded potential and pressure tensors from explicit MTS steps.
SUBROUTINE MTS_SAVEAUXS(I)

  USE MTS
  USE VAR

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I !< The current mod of the MTS step (1, 2, or 3 here)
  
  MTS_V_NB(I, 1) = V_NB(1)
  MTS_V_NB(I, 2) = V_NB(2)
  MTS_V_NB(I, 3) = V_NB(3)

  MTS_PT11(I) = PT11
  MTS_PT22(I) = PT22
  MTS_PT33(I) = PT33
  MTS_PT12(I) = PT12
  MTS_PT13(I) = PT13
  
  RETURN
END SUBROUTINE MTS_SAVEAUXS

!> @brief Approximates the nonbonded potential and pressure tensors in approximate MTS steps.
SUBROUTINE MTS_LOADAUXS(I)

  USE MTS
  USE VAR

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I

  V_NB(1) = MTS_V_NB(3, 1) + (I-3) * 0.5 * (MTS_V_NB(3, 1) - MTS_V_NB(1, 1))
  V_NB(2) = MTS_V_NB(3, 2) + (I-3) * 0.5 * (MTS_V_NB(3, 2) - MTS_V_NB(1, 2))
  V_NB(3) = MTS_V_NB(3, 3) + (I-3) * 0.5 * (MTS_V_NB(3, 3) - MTS_V_NB(1, 3))

  PT11 = MTS_PT11(3) + (I-3) * 0.5 * (MTS_PT11(3) - MTS_PT11(1))
  PT22 = MTS_PT22(3) + (I-3) * 0.5 * (MTS_PT22(3) - MTS_PT22(1))
  PT33 = MTS_PT33(3) + (I-3) * 0.5 * (MTS_PT33(3) - MTS_PT33(1))
  PT12 = MTS_PT12(3) + (I-3) * 0.5 * (MTS_PT12(3) - MTS_PT12(1))
  PT13 = MTS_PT13(3) + (I-3) * 0.5 * (MTS_PT13(3) - MTS_PT13(1))
  PT23 = MTS_PT23(3) + (I-3) * 0.5 * (MTS_PT23(3) - MTS_PT23(1))

  RETURN
END SUBROUTINE MTS_LOADAUXS

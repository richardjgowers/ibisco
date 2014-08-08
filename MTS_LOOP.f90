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

  ALLOCATE(MTS_FXYZ(3, 3, NCOARSE))
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

  DEALLOCATE(MTS_FXYZ)
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
  FXYZ = 0.0

  FXYZNB = 0.0

  PT11 = 0.0
  PT22 = 0.0
  PT33 = 0.0
  PT12 = 0.0
  PT13 = 0.0
  PT23 = 0.0

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

        CALL MTS_SAVEFORCE(MTS_MOD, BEAD, FXYZ, NITEMS, MTS_FXYZ, NCOARSE)
        CALL MTS_SAVEAUXS(MTS_MOD)

     CASE DEFAULT ! Use approximation
        CALL MTS_APPROX(MTS_MOD, BEAD, FXYZ, NITEMS, MTS_FXYZ, NCOARSE)
        CALL MTS_LOADAUXS(MTS_MOD)
  END SELECT
  CALL DISTRIBUTE_VSFORCE()

  ! Call atoms second, so that V_NB and PTxx only have CG contributions
  ! Order shouldnt matter anyway...
  CALL NONBONDED_FORCE(NUMATOMS, ATOM, MAXNAB_ATOM, LIST_ATOM, RCUT_ATOM, RCUTSQ_ATOM)

  DO I = 1, NATOMS
     FXYZNB(1,I) = FXYZ(1,I)
     FXYZNB(2,I) = FXYZ(2,I)
     FXYZNB(3,I) = FXYZ(3,I)
  END DO

  CALL BONDED_FORCE()

  CALL DISTRIBUTE_VSFORCE()
  
  DO I=1,NATOMS
     PT11 = PT11 + (FXYZ(1,I) - FXYZNB(1,I))*SXYZ(1,I)
     PT22 = PT22 + (FXYZ(2,I) - FXYZNB(2,I))*SXYZ(2,I)
     PT33 = PT33 + (FXYZ(3,I) - FXYZNB(3,I))*SXYZ(3,I)
     PT12 = PT12 + (FXYZ(2,I) - FXYZNB(2,I))*SXYZ(1,I)
     PT13 = PT13 + (FXYZ(3,I) - FXYZNB(3,I))*SXYZ(1,I)
     PT23 = PT23 + (FXYZ(3,I) - FXYZNB(3,I))*SXYZ(2,I)
  END DO

  RETURN
END SUBROUTINE MTS_FORCE

!> @brief Records the coarse grained force from an explicit MTS step.
SUBROUTINE MTS_SAVEFORCE(I, BEAD, FXYZ, NITEMS, MTS_FXYZ, NCOARSE)

  USE VAR, ONLY : RKIND

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I !< Mod of current MTS step
  INTEGER, INTENT(IN) :: NITEMS, NCOARSE !< Total number of particles (real and virtual)
  REAL*4, DIMENSION(3, NITEMS), INTENT(IN) :: FXYZ
  REAL*4, DIMENSION(3, 3, NCOARSE), INTENT(INOUT) :: MTS_FXYZ
  INTEGER, DIMENSION(NCOARSE), INTENT(IN) :: BEAD
  
  INTEGER :: A, ATOM_ID

  DO A = 1, NCOARSE
     ATOM_ID = BEAD(A)
     MTS_FXYZ(1,I, A) = FXYZ(1,ATOM_ID)
     MTS_FXYZ(2,I, A) = FXYZ(2,ATOM_ID)
     MTS_FXYZ(3,I, A) = FXYZ(3,ATOM_ID)
  END DO
  
  RETURN
  
END SUBROUTINE MTS_SAVEFORCE

!> @brief Generates approximated coarse grained forces for MTS steps.
SUBROUTINE MTS_APPROX(I, BEAD, FXYZ, NITEMS, MTS_FXYZ, NCOARSE)

  USE VAR, ONLY : RKIND

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I ! Mod of current MTS step
  INTEGER, INTENT(IN) :: NITEMS, NCOARSE !< Total number of particles (real and virtual)
  REAL*4, DIMENSION(3, NITEMS), INTENT(INOUT) :: FXYZ
  REAL*4, DIMENSION(3, 3, NCOARSE), INTENT(IN) :: MTS_FXYZ
  INTEGER, DIMENSION(NCOARSE), INTENT(IN) :: BEAD

  INTEGER :: A, ATOM_ID

  DO A = 1, NCOARSE
     ATOM_ID = BEAD(A)
     ! Forwards approx
     FXYZ(1,ATOM_ID) = MTS_FXYZ(1,3,A) + (I-3) * 0.5 * (MTS_FXYZ(1,3,A) - MTS_FXYZ(1,1,A)) 
     FXYZ(2,ATOM_ID) = MTS_FXYZ(2,3,A) + (I-3) * 0.5 * (MTS_FXYZ(2,3,A) - MTS_FXYZ(2,1,A)) 
     FXYZ(3,ATOM_ID) = MTS_FXYZ(3,3,A) + (I-3) * 0.5 * (MTS_FXYZ(3,3,A) - MTS_FXYZ(3,1,A)) 
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

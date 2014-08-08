#ifndef DOXY_SKIP
#include "ibi-preprocess.h"
#endif
!> @file
!> @brief The loop which calculates all force information on all particles

SUBROUTINE NEW_FORCE()

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER :: I

  !Reset all forces and pressure
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

#ifdef TIMING_ON
  t_NONBONDED_atom(1) = OMP_GET_WTIME()
#endif
  !Nonbonded forces
  CALL NONBONDED_FORCE(NUMATOMS,ATOM,MAXNAB_ATOM,LIST_ATOM,RCUT_ATOM,RCUTSQ_ATOM)
#ifdef TIMING_ON
  t_NONBONDED_atom(2) = OMP_GET_WTIME()
#endif

  IF(IBRDESCR .eq. 0) THEN
#ifdef TIMING_ON
     t_NONBONDED_bead(1) = OMP_GET_WTIME()
#endif
     CALL NONBONDED_FORCE(NCOARSE,BEAD,MAXNAB_BEAD,LIST_BEAD,RCUT_BEAD,RCUTSQ_BEAD)
#ifdef TIMING_ON
     t_NONBONDED_bead(2) = OMP_GET_WTIME()
     t_DISTRIBUTE_VSFORCE(1) = OMP_GET_WTIME()
#endif
     CALL DISTRIBUTE_VSFORCE() !Distributes forces from VS onto the atoms underneath
#ifdef TIMING_ON
     t_DISTRIBUTE_VSFORCE(2) = OMP_GET_WTIME()
#endif
  END IF

  DO I=1,NATOMS
     FXYZNB(1,I) = FXYZ(1,I)
     FXYZNB(2,I) = FXYZ(2,I)
     FXYZNB(3,I) = FXYZ(3,I)
  END DO

#ifdef TIMING_ON
  t_BONDED_FORCE(1) = OMP_GET_WTIME()
#endif
  !Bonded forces
  CALL BONDED_FORCE()
#ifdef TIMING_ON
  t_BONDED_FORCE(2) = OMP_GET_WTIME()
#endif
  CALL DISTRIBUTE_VSFORCE() !Distribute forces on virtual sites arising from bonded interactions

  DO I=1,NATOMS
     PT11 = PT11 + (FXYZ(1,I) - FXYZNB(1,I)) * SXYZ(1,I)
     PT22 = PT22 + (FXYZ(2,I) - FXYZNB(2,I)) * SXYZ(2,I)
     PT33 = PT33 + (FXYZ(3,I) - FXYZNB(3,I)) * SXYZ(3,I)
     PT12 = PT12 + (FXYZ(2,I) - FXYZNB(2,I)) * SXYZ(1,I)
     PT13 = PT13 + (FXYZ(3,I) - FXYZNB(3,I)) * SXYZ(1,I)
     PT23 = PT23 + (FXYZ(3,I) - FXYZNB(3,I)) * SXYZ(2,I)
  END DO

  RETURN
END SUBROUTINE NEW_FORCE

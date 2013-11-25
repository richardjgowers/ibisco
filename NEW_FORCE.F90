#include "ibi-preprocess.h"

SUBROUTINE NEW_FORCE()

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER :: I

  !Reset all forces and pressure
  FX = 0.0D0
  FY = 0.0D0
  FZ = 0.0D0

  FXNB      = 0.0D0
  FYNB      = 0.0D0
  FZNB      = 0.0D0

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

#ifdef TIMING_ON
  t_NONBONDED_atom(1) = OMP_GET_WTIME()
#endif
  !Nonbonded forces
  CALL NONBONDED_FORCE(NUMATOMS,ATOM,MAXNAB_ATOM,LIST_ATOM,RCUT_ATOM,RCUTSQ_ATOM)

#ifdef TIMING_ON
  t_NONBONDED_atom(2) = OMP_GET_WTIME()
  t_NONBONDED_bead(1) = OMP_GET_WTIME()
#endif

  IF(IBRDESCR .eq. 0) THEN
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
     FXNB(I) = FX(I)
     FYNB(I) = FY(I)
     FZNB(I) = FZ(I)
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
     PT11 = PT11 + (FX(I) - FXNB(I))*SX(I)
     PT22 = PT22 + (FY(I) - FYNB(I))*SY(I)
     PT33 = PT33 + (FZ(I) - FZNB(I))*SZ(I)
     PT12 = PT12 + (FY(I) - FYNB(I))*SX(I)
     PT13 = PT13 + (FZ(I) - FZNB(I))*SX(I)
     PT23 = PT23 + (FZ(I) - FZNB(I))*SY(I)
  END DO

  RETURN
END SUBROUTINE NEW_FORCE

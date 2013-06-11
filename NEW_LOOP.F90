#include "ibi-preprocess.h"

SUBROUTINE NEW_LOOP()

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  DO STEP=1,NSTEP
#ifdef TIMING_ON
!Reset timer functions to zero at start of each loop
     t_LOOP = 0.0D0

     t_FORCECALC  = 0.0D0
     t_VIRTUAL_DEF = 0.0D0
     t_UPDATE_NEIGHBOURLIST  = 0.0D0
     t_MAPS_atom = 0.0D0
     t_MAPS_bead = 0.0D0
     t_LINKS_atom = 0.0D0
     t_LINKS_bead = 0.0D0
     t_NEW_NEIGHBOUR_atom = 0.0D0
     t_NEW_NEIGHBOUR_bead = 0.0D0
     t_NEW_FORCE = 0.0D0
     t_NONBONDED_atom = 0.0D0
     t_NONBONDED_bead = 0.0D0
     t_DISTRIBUTE_VSFORCE = 0.0D0

     t_BONDED_FORCE = 0.0D0
     t_BONDS = 0.0D0
     t_ANGLES = 0.0D0
     t_TORSIONS = 0.0D0
     t_OUTOFPLANES = 0.0D0

     t_SHIFT = 0.0D0
     t_MOVE = 0.0D0
     t_MOMENTUM = 0.0D0
     t_SCALEBP = 0.0D0

     t_REPORTRESULTS = 0.0D0
     t_AVERAGE = 0.0D0
     t_WRITETRJ = 0.0D0
     t_OUTPUT = 0.0D0

     t_LOOP(1) = OMP_GET_WTIME()
     t_SHIFT(1) = OMP_GET_WTIME()
#endif

     !Update shifted positions of atoms in box
     CALL SHIFT() !Applies PBC

#ifdef TIMING_ON
     t_SHIFT(2) = OMP_GET_WTIME()
     t_FORCECALC(1) = OMP_GET_WTIME()
#endif

     IF(IBRDESCR .eq. 0 .and. MOD(STEP,VUPDATE) .eq. 0) THEN

#ifdef TIMING_ON
        t_VIRTUAL_DEF(1) = OMP_GET_WTIME()
#endif
        CALL VIRTUAL_DEF() !Defines the position of virtual sites

#ifdef TIMING_ON
        t_VIRTUAL_DEF(2)= OMP_GET_WTIME()    
#endif
     END IF

     !If needed, update neighbour list
     IF(MOD(STEP,NUPDATE) .eq. 0) THEN

#ifdef TIMING_ON
        t_UPDATE_NEIGHBOURLIST(1) = OMP_GET_WTIME()
#endif

        CALL UPDATE_NEIGHBOURLIST()

#ifdef TIMING_ON
        t_UPDATE_NEIGHBOURLIST(2) = OMP_GET_WTIME()
#endif

     END IF

     !Calculate all forces on atoms
#ifdef TIMING_ON
        t_NEW_FORCE(1) = OMP_GET_WTIME()
#endif
     CALL NEW_FORCE()

#ifdef TIMING_ON
        t_NEW_FORCE(2) = OMP_GET_WTIME()
        t_FORCECALC(2) = OMP_GET_WTIME()
        t_MOVE(1) = OMP_GET_WTIME()
#endif

     !Move atoms within box
     CALL MOVE()

     VOLUME = BOXX * BOXY * BOXZ
     INV_VOLUME = 1.0D0 / VOLUME

     PT11 = PT11 * INV_VOLUME
     PT22 = PT22 * INV_VOLUME
     PT33 = PT33 * INV_VOLUME

     PT12 = PT12 * INV_VOLUME
     PT13 = PT13 * INV_VOLUME
     PT23 = PT23 * INV_VOLUME

#ifdef TIMING_ON
        t_MOVE(2) = OMP_GET_WTIME()
#endif

     !If needed, halt the net drift of box
     IF(MOD(STEP,HALT_DRIFT) .eq. 0) THEN
#ifdef TIMING_ON
        t_MOMENTUM(1) = OMP_GET_WTIME()
#endif
        CALL MOMENTUM()
#ifdef TIMING_ON
        t_MOMENTUM(2) = OMP_GET_WTIME()
#endif
     END IF

     IF(ENSEMBLE .eq. 2) THEN !If NPT
#ifdef TIMING_ON
        t_SCALEBP(1) = OMP_GET_WTIME()
#endif
        CALL SCALEBP(STEP) !Change box size and scale positions
#ifdef TIMING_ON
        t_SCALEBP(2) = OMP_GET_WTIME()
#endif
     END IF

#ifdef TIMING_ON
     t_REPORTRESULTS(1) = OMP_GET_WTIME()
#endif

     !STORING AVERAGE DATA AND RESTART FILE
     IF ((STEP .EQ. 1).OR.(MOD(STEP, NAVERAGE) == 0)) THEN
#ifdef TIMING_ON
        t_AVERAGE(1) = OMP_GET_WTIME()
#endif
        CALL AVERAGE (STEP)
#ifdef TIMING_ON
        t_AVERAGE(2) = OMP_GET_WTIME()
#endif
     END IF

     !STORING THE TRAJECTORY FILE
     IF (MOD(STEP, NTRJ) == 0) THEN
#ifdef TIMING_ON
        t_WRITETRJ(1) = OMP_GET_WTIME()
#endif
        CALL WRITETRJ (STEP)
#ifdef TIMING_ON
        t_WRITETRJ(2) = OMP_GET_WTIME()
#endif
     END IF

     IF (MOD(STEP, NSAMPLING) == 0) THEN
#ifdef TIMING_ON
        t_OUTPUT(1) = OMP_GET_WTIME()
#endif
        CALL OUTPUT(STEP)
#ifdef TIMING_ON
        t_OUTPUT(2) = OMP_GET_WTIME()
#endif
     END IF

#ifdef TIMING_ON
     t_REPORTRESULTS(2) = OMP_GET_WTIME()
     t_LOOP(2) = OMP_GET_WTIME()

     CALL TIMING()
#endif

  END DO

  RETURN
END SUBROUTINE NEW_LOOP

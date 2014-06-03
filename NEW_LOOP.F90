#ifndef DOXY_SKIP
#include "ibi-preprocess.h"
#endif
!> @file
!> @brief The molecular dynamics loop
!> @details Called from IBIsCO()
!> @details Calls a series of subroutines to perform the MD loop \n
!!          SHIFT() \n
!!          VIRTUAL_DEF() If using hybrid system \n
!!          UPDATE_NEIGHBOURLIST() If required (updates every n steps) \n
!!          NEW_FORCE() Main calculation step \n
!!          MOMENTUM() If required \n
!!          SCALEBP() If NPT \n
!!          AVERAGE() Keeps track of average physical values for reporting \n
!!          WRITETRJ() Writes the trajectory \n
!!          OUTPUT() Writes the s-md.out file \n
!!          TIMING() If compiled using TIMING flag, reports on time in each subroutine \n

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

     t_LOOP(1) = OMP_GET_WTIME()
     t_SHIFT(1) = OMP_GET_WTIME()
#endif

     !Update shifted positions of atoms in box
     CALL SHIFT() !Applies PBC

#ifdef TIMING_ON
     t_SHIFT(2) = OMP_GET_WTIME()
     t_FORCECALC(1) = OMP_GET_WTIME()
     t_VIRTUAL_DEF(1) = OMP_GET_WTIME()
#endif

     IF(IBRDESCR .eq. 0) THEN
        CALL VIRTUAL_DEF() !Defines the position of virtual sites
     END IF

#ifdef TIMING_ON
     t_VIRTUAL_DEF(2)= OMP_GET_WTIME() 
     t_UPDATE_NEIGHBOURLIST(1) = OMP_GET_WTIME()
#endif

     !If needed, update neighbour list
     IF(MOD(STEP,NUPDATE) .eq. 0) THEN
        CALL UPDATE_NEIGHBOURLIST()
     END IF

#ifdef TIMING_ON
     t_UPDATE_NEIGHBOURLIST(2) = OMP_GET_WTIME()
     t_NEW_FORCE(1) = OMP_GET_WTIME()
#endif

     !Calculate all forces on atoms
     CALL NEW_FORCE()

#ifdef TIMING_ON
     t_NEW_FORCE(2) = OMP_GET_WTIME()
     t_FORCECALC(2) = OMP_GET_WTIME()
     t_MOVE(1) = OMP_GET_WTIME()
#endif

     !Move atoms within box
     CALL MOVE()

#ifdef TIMING_ON
     t_MOVE(2) = OMP_GET_WTIME()
     t_MOMENTUM(1) = OMP_GET_WTIME()
#endif

     !If needed, halt the net drift of box
     IF(MOD(STEP,HALT_DRIFT) .eq. 0) THEN
        CALL MOMENTUM()
     END IF

#ifdef TIMING_ON
     t_MOMENTUM(2) = OMP_GET_WTIME()
     t_SCALEBP(1) = OMP_GET_WTIME()
#endif

     IF(ENSEMBLE .eq. 2) THEN !If NPT
        CALL SCALEBP(STEP) !Change box size and scale positions
     END IF

#ifdef TIMING_ON
     t_SCALEBP(2) = OMP_GET_WTIME()
     t_REPORTRESULTS(1) = OMP_GET_WTIME()
     t_AVERAGE(1) = OMP_GET_WTIME()
#endif

     !STORING AVERAGE DATA AND RESTART FILE
     CALL AVERAGE (STEP)


#ifdef TIMING_ON
     t_AVERAGE(2) = OMP_GET_WTIME()
     t_WRITETRJ(1) = OMP_GET_WTIME()
#endif

     !STORING THE TRAJECTORY FILE
     IF (MOD(STEP, NTRJ) == 0) THEN
        CALL WRITETRJ (STEP)
        CALL OUTPUT(STEP)
     END IF

#ifdef TIMING_ON
     t_WRITETRJ(2) = OMP_GET_WTIME()
     t_REPORTRESULTS(2) = OMP_GET_WTIME()
     t_LOOP(2) = OMP_GET_WTIME()
     !Write all timing results to file
     CALL TIMING()
#endif

  END DO

  RETURN
END SUBROUTINE NEW_LOOP

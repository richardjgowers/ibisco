SUBROUTINE NEW_LOOP()

  USE VAR

  IMPLICIT NONE

  INTEGER :: I

  DO STEP=1,NSTEP

     !Update shifted positions of atoms in box
     CALL SHIFT() !Applies PBC
     IF(IBRDESCR .eq. 0) CALL VIRTUAL_DEF() !Defines the position of virtual sites

     !If needed, update neighbour list
     IF(MOD(STEP,NUPDATE) .eq. 0) THEN
        CALL UPDATE_NEIGHBOURLIST()
     END IF

     !Calculate all forces on atoms
     CALL NEW_FORCE()

     !Move atoms within box
     CALL MOVE()

     !If needed, halt the net drift of box
     IF(MOD(STEP,HALT_DRIFT) .eq. 0) THEN
        CALL MOMENTUM()
     END IF

     IF(ENSEMBLE .eq. 2) THEN !If NPT
        CALL SCALEBP(STEP) !Change box size and scale positions
     END IF

     !STORING AVERAGE DATA AND RESTART FILE
     IF ((STEP .EQ. 1).OR.(MOD(STEP, NAVERAGE) == 0)) THEN
        CALL AVERAGE (STEP)
     END IF

     !STORING THE TRAJECTORY FILE
     IF (MOD(STEP, NTRJ) == 0) THEN
        CALL WRITETRJ (STEP)
     END IF

     IF (MOD(STEP, NSAMPLING) == 0) THEN
        CALL OUTPUT(STEP)
     END IF

  END DO

  RETURN
END SUBROUTINE NEW_LOOP

        SUBROUTINE LOOP (I)
        USE VAR
        USE RNEMD
        IMPLICIT        NONE
        INTEGER I

!       DO 25 I = 1, NSTEP
		PT11 = 0.0D0
		PT22 = 0.0D0
		PT33 = 0.0D0

		PT12 = 0.0D0
		PT13 = 0.0D0
		PT23 = 0.0D0

!		SHIFT THE ATOMS INSIDE THE BOX 
		CALL SHIFT ()

!		MAKEING THE NEIGHBOUR LIST
    IF (MOD(I, NUPDATE) == 0) THEN
        IF (NEIGHBORLIST == 'NEIGHBOR_NOLIST') THEN
            CALL NEIGHBOR_NOLIST ()
        ELSE
            CALL NEIGHBOR_WITHLIST ()
        END IF
    END IF
	
!		CALCULATE THE BONDED AND NON-BONDED INTERACTIONS
    IF (PPF_INPUT.EQ.1) THEN
        CALL FORCE2 ()
        call MOVE2()
    ELSE
        CALL FORCE()
        if (VISCINPUT == 1) then
            CALL MOVERNEMD()
        ELSE
            if(.not. shakeOK)then
                CALL MOVE()
            else
                CALL SHAKE()
            end if
        ENDIF
    ENDIF
                      
		
!		UPDATE THE POSITIONS AND VELOCITIES OF THE ATOMS AND CALCULATE THE KINETIC TERM OF PRESSUR
             
                IF (VISCINPUT == 1) then
                 IF (I.EQ.TSET) THEN
                    NUMPRO = 0
                    NSWP = 0
                    TTRANSF = 0.D0
                    VXMEAN_SLAB = 0.D0
                    SLAB_MDENS = 0.D0
                    SLAB_MTEMP = 0.D0
                 ENDIF
                 call VISC_NEMD (I)                                    ! RNEMD shear viscosity
                ENDIF
              

	!		CALCULATE THE DIAGONAL COMPONENTS OF PRESSURE TENSOR
			VOLUME = BOXX * BOXY * BOXZ
			PT11 = PT11 / VOLUME
			PT22 = PT22 / VOLUME
			PT33 = PT33 / VOLUME

			PT12 = PT12 / VOLUME
			PT13 = PT13 / VOLUME
			PT23 = PT23 / VOLUME

	!		RESET THE NET MOMENTUM TO ZERO
                        IF (DPDINPUT.NE.1.AND.LAINPUT.NE.1) THEN
                        IF (MOD(I, HALT_DRIFT) == 0) THEN
			CALL MOMENTUM()
			END IF
                        END IF

	!		CALCULATE THE NEW VOLUME AND SCALE THE POSITIONS IN NPT ENSEMBLE
	 		IF (ENSEMBLE == 2) THEN
			CALL SCALEBP (I)
	 		END IF

	!		STORING AVERAGE DATA AND RESTART FILE
			IF ((I .EQ. 1).OR.(MOD(I, NAVERAGE) == 0)) THEN
			CALL AVERAGE (I)
			END IF

	!		STORING THE TRAJECTORY FILE
			IF (MOD(I, NTRJ) == 0) THEN		
			CALL WRITETRJ (I)
			END IF

			IF (MOD(I, NSAMPLING) == 0) THEN
			CALL OUTPUT (I)
                        END IF
                   
                        IF (PPF_INPUT.EQ.1) THEN
                        
                        CALL PBCPOISEUILLE (I)
        
                        IF (MOD(I, NSAMPLING) == 0) THEN
	                write(173,371) I,I*DT*TIMESCALE,VISCOSITY,RVISCOSITY
                        call flush (173)
                        END IF
                        END IF
!	25	CONTINUE
      
 371 FORMAT (I10,3(1X,1E12.5))
        RETURN
        END


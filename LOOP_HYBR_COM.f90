!  
!	Modified by Nicodemo
!		Dec 2010
!
!
!	UPDATE:
!	line 34: Neighbour list updating: added the choose 
!	beetwen the neighbour list with or without atoms depending on the description
!	(hybrid or not)
!
!     NOUP = if the position of the virtual site is updated because the neighbour 
!           list is useless to update again the position of virtual site.
!     NOUP = 1, the virtual site has just been updated when is updated the neighbour list
!     NOUP = 0, the virtual site is updated 

      SUBROUTINE LOOP_HYBR_COM ()
      USE VAR
      USE RNEMD
      IMPLICIT NONE
      INTEGER :: I,IL, JL, HL
      INTEGER :: NOUP = 0
      REAL(KIND=RKIND) :: SUMTOTX, SUMTOTY, SUMTOTZ
      REAL(KIND=RKIND) :: SPX, SPY, SPZ
      REAL(KIND=RKIND) :: RPX,RPY,RPZ


do i = 1, NSTEP

    timestepcheck = i
    PT11 = 0.0D0
    PT22 = 0.0D0
    PT33 = 0.0D0
    PT12 = 0.0D0
    PT13 = 0.0D0
    PT23 = 0.0D0

!		SHIFT THE ATOMS INSIDE THE BOX 
    CALL SHIFT ()

!		MAKING THE NEIGHBOUR LIST
    IF (MOD(I, NUPDATE) == 0) THEN
        IF (NEIGHBORLIST == 'NEIGHBOR_NOLIST') THEN
            write(*,*)' DEVI ANCORA SCRIVERE LA FUNZIONE VIRT_NEIGHBOUR_NOLIST_COM!!!'
            stop
!                             CALL VIRT_NEIGHBOR_NOLIST_COM ()
        ELSE
! Before the calling to the function: neighbour_withlist, the program update the positions of the center of mass
            DO IL=1,NVIRTA
                SUMTOTX = 0.0D0
                SUMTOTY = 0.0D0
                SUMTOTZ = 0.0D0
                DO HL = 1,init_numbcomp(il)
                    SUMTOTX = SUMTOTX+MASS(ITYPE(COMPCOM(IL,HL)))*SX(COMPCOM(IL,HL))
                    SUMTOTY = SUMTOTY+MASS(ITYPE(COMPCOM(IL,HL)))*SY(COMPCOM(IL,HL))
                    SUMTOTZ = SUMTOTZ+MASS(ITYPE(COMPCOM(IL,HL)))*SZ(COMPCOM(IL,HL))
                END DO
                VIRTRX(IL) = SUMTOTX*INVTOTBMASS(IL)
                VIRTRY(IL) = SUMTOTY*INVTOTBMASS(IL)
                VIRTRZ(IL) = SUMTOTZ*INVTOTBMASS(IL)
                IF ((VIRTRX(IL)> BOXX2).OR.(VIRTRX(IL) < -BOXX2).OR.(VIRTRY(IL)> BOXY2) &
                      .OR.(VIRTRY(IL)< -BOXY2).OR.(VIRTRZ(IL)> BOXZ2).OR.(VIRTRZ(IL)<-BOXZ2)) THEN
                   RPX = VIRTRX(IL)/BOXX2
                   RPY = VIRTRY(IL)/BOXY2
                   RPZ = VIRTRZ(IL)/BOXZ2
                   SPX = RPX - 2.0D0*DBLE(REAL(INT(RPX/2.0D0)))
                   SPY = RPY - 2.0D0*DBLE(REAL(INT(RPY/2.0D0)))
                   SPZ = RPZ - 2.0D0*DBLE(REAL(INT(RPZ/2.0D0)))
                   RPX = (SPX - 2.0D0*DBLE(INT(SPX)))*BOXX2
                   RPY = (SPY - 2.0D0*DBLE(INT(SPY)))*BOXY2
                   RPZ = (SPZ - 2.0D0*DBLE(INT(SPZ)))*BOXZ2
                   VIRTRX(IL) = RPX
                   VIRTRY(IL) = RPY
                   VIRTRZ(IL) = RPZ 
                END IF 
            END DO
            CALL VIRT_NEIGHBOR_WITHLIST_COM ()
            NOUP = 1
        END IF
    END IF


IF (MOD(I, VUPDATE) .EQ. 0) THEN
    IF (NOUP .EQ. 0) THEN
        DO IL=1,NVIRTA
            SUMTOTX = 0.0D0
            SUMTOTY = 0.0D0
            SUMTOTZ = 0.0D0
            DO HL = 1,init_numbcomp(il)
                SUMTOTX = SUMTOTX+MASS(ITYPE(COMPCOM(IL,HL)))*SX(COMPCOM(IL,HL))
                SUMTOTY = SUMTOTY+MASS(ITYPE(COMPCOM(IL,HL)))*SY(COMPCOM(IL,HL))
                SUMTOTZ = SUMTOTZ+MASS(ITYPE(COMPCOM(IL,HL)))*SZ(COMPCOM(IL,HL))
            END DO
            VIRTRX(IL) = SUMTOTX*INVTOTBMASS(IL)
            VIRTRY(IL) = SUMTOTY*INVTOTBMASS(IL)
            VIRTRZ(IL) = SUMTOTZ*INVTOTBMASS(IL)
            IF ((VIRTRX(IL)> BOXX2).OR.(VIRTRX(IL) < -BOXX2).OR.(VIRTRY(IL)> BOXY2) &
                  .OR.(VIRTRY(IL)< -BOXY2).OR.(VIRTRZ(IL)> BOXZ2).OR.(VIRTRZ(IL)<-BOXZ2)) THEN
                RPX = VIRTRX(IL)/BOXX2
                RPY = VIRTRY(IL)/BOXY2
                RPZ = VIRTRZ(IL)/BOXZ2
                SPX = RPX - 2.0D0*DBLE(REAL(INT(RPX/2.0D0)))
                SPY = RPY - 2.0D0*DBLE(REAL(INT(RPY/2.0D0)))
                SPZ = RPZ - 2.0D0*DBLE(REAL(INT(RPZ/2.0D0)))
                RPX = (SPX - 2.0D0*DBLE(INT(SPX)))*BOXX2
                RPY = (SPY - 2.0D0*DBLE(INT(SPY)))*BOXY2
                RPZ = (SPZ - 2.0D0*DBLE(INT(SPZ)))*BOXZ2
                VIRTRX(IL) = RPX
                VIRTRY(IL) = RPY
                VIRTRZ(IL) = RPZ 	
            END IF 
        END DO
    END IF
    NOUP = 0
END IF


!		CALCULATE THE BONDED AND NON-BONDED INTERACTIONS

    IF (PPF_INPUT.EQ.1) THEN
        CALL FORCE2 ()
        CALL MOVE2()
    ELSE
        CALL FORCE_HYBR_COM()
        IF (VISCINPUT == 1) THEN
            CALL MOVERNEMD()
        ELSE
            CALL MOVE()
        END IF
    END IF  

!		UPDATE THE POSITIONS AND VELOCITIES OF THE ATOMS AND CALCULATE THE KINETIC TERM OF PRESSUR
             
            IF (VISCINPUT == 1) THEN
                  IF (I.EQ.TSET) THEN
                        NUMPRO = 0
                        NSWP = 0
                        TTRANSF = 0.D0
                        VXMEAN_SLAB = 0.D0
                        SLAB_MDENS = 0.D0
                        SLAB_MTEMP = 0.D0
                  ENDIF
                  CALL VISC_NEMD (I)                                    ! RNEMD shear viscosity
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
                        CALL HAVERAGE (I)
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
                              WRITE(173,371) I,DBLE(I)*DT*TIMESCALE,VISCOSITY,RVISCOSITY
                              CALL FLUSH (173)
                        END IF
                  END IF

end do

      371 FORMAT (I10,3(1X,1E12.5))
      RETURN
      END


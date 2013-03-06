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


      SUBROUTINE LOOP_HYBR_MTS4_COM ()
      USE VAR
      USE RNEMD
      IMPLICIT NONE
      INTEGER :: I,IL, JL, HL,IP, HH
      INTEGER :: NOUP = 0
      REAL(KIND=RKIND) :: SUMTOTX, SUMTOTY, SUMTOTZ, K1,K2,K3,K4
      REAL(KIND=RKIND) :: SPX, SPY, SPZ
      REAL(KIND=RKIND) :: RPX,RPY,RPZ

      REAL(KIND=RKIND) :: FXfirst,FXsec,FXthird,FXfourth,FYfirst,FYsec,FYthird, &
                              FYfourth,FZfirst,FZsec,FZthird,FZfourth,fxfifth,fyfifth,fzfifth

      !real(KIND=RKIND) :: ex,ey,ez,moderr,ggg,Fmod
      !real,dimension(10000) :: errFx,errFy,errFz,merr

do I = 1, NSTEP

    timestepcheck = i
    Kmts = Kmts+1
    K1   = real(Kmts)
    K2   = (real(Kmts)**2)
    K3   = (real(Kmts)**3)
    K4   = (real(Kmts)**4)

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
!write(*,*)i,kmts
    IF (I .LT. MTSPARAM2) THEN
        CALL FORCE_HYBR_PREMTS_COM ()
        IF(timestepcheck .EQ. I0) THEN
            I0 = I0+MTSPARAM+TSMTS
! write(*,*)i0,'I0',i
        ELSE IF (timestepcheck .EQ. I1) THEN
            I1 = I1+MTSPARAM+TSMTS
!write(*,*)i1,'I1',i
        ELSE IF (timestepcheck .EQ. I2) THEN
            I2 = I2+MTSPARAM+TSMTS
!write(*,*)i2,'I2',i
        ELSE IF (timestepcheck .EQ. I3) THEN 
            I3 = I3+MTSPARAM+TSMTS
!write(*,*)i3,'I3',i
        ELSE IF (timestepcheck .EQ. I4) THEN
            I4 = I4+MTSPARAM+TSMTS
!write(*,*)i4,'I4',i
        ELSE IF (timestepcheck .EQ. I5) THEN
            I5 = I5+MTSPARAM+TSMTS
!write(*,*)i5,'I5',i
        END IF
    ELSE
        DO HH =1,NUM_BEAD
           IP = INDEX_AB(HH)
           CALL VIRT_FORCE_COM_MTS(IP,FCUTB)
! TAYLOR 4: 5 points
            FXfirst  = (-0.25*FXiiii(IP)+FXiii(IP)*1.333333333333333-3*FXii(IP)+ &
                            4*FXi1(IP)-FXi0(IP)*2.0833333333333)*idt
            FYfirst  = (-0.25*FYiiii(IP)+FYiii(IP)*1.333333333333333-3*FYii(IP)+ &
                            4*FYi1(IP)-FYi0(IP)*2.0833333333333)*idt
            FYfirst  = (-0.25*FZiiii(IP)+FZiii(IP)*1.333333333333333-3*FZii(IP)+ & 
                            4*FZi1(IP)-FZi0(IP)*2.0833333333333)*idt
            FXsec    = (FXiiii(IP)*0.916666667-FXiii(IP)*4.666666667+FXii(IP)*9.5-&
                            FXi1(IP)*8.666666667+FXi0(IP)*2.916666667)*idt2
            FYsec    = (FYiiii(IP)*0.916666667-FYiii(IP)*4.666666667+FYii(IP)*9.5-&
                            FYi1(IP)*8.666666667+FYi0(IP)*2.916666667)*idt2
            FZsec    = (FZiiii(IP)*0.916666667-FZiii(IP)*4.666666667+FZii(IP)*9.5-&
                            FZi1(IP)*8.666666667+FZi0(IP)*2.916666667)*idt2
            FXthird  = (-1.5*FXiiii(IP)+7*FXiii(IP)-12*FXii(IP)+9*FXi1(IP)-2.5*FXi0(IP))*idt3
            FYthird  = (-1.5*FYiiii(IP)+7*FYiii(IP)-12*FYii(IP)+9*FYi1(IP)-2.5*FYi0(IP))*idt3
            FZthird  = (-1.5*FZiiii(IP)+7*FZiii(IP)-12*FZii(IP)+9*FZi1(IP)-2.5*FZi0(IP))*idt3

            FXfourth = (FXiiii(IP) - 4*FXiii(IP) + 6*FXii(IP) - 4*FXi1(IP) + FXi0(IP))*idt4
            FYfourth = (FYiiii(IP) - 4*FYiii(IP) + 6*FYii(IP) - 4*FYi1(IP) + FYi0(IP))*idt4
            FZfourth = (FZiiii(IP) - 4*FZiii(IP) + 6*FZii(IP) - 4*FZi1(IP) + FZi0(IP))*idt4

! Force calculation with Taylor 4

            FX(IP) = FXiiii(IP) + FXfirst*DT*K1 + FXsec*mtsDT2*K2 + &
                        FXthird*DT3*K3 + FXfourth*DT4*K4
            FY(IP) = FYiiii(IP) + FYfirst*DT*K1 + FYsec*mtsDT2*K2 + &
                        FYthird*DT3*K3 + FYfourth*DT4*K4
            FZ(IP) = FZiiii(IP) + FZfirst*DT*K1 + FZsec*mtsDT2*K2 + &
                        FZthird*DT3*K3 + FZfourth*DT4*K4
!        write(1000,*)ip,FX(IP),FY(IP),FZ(IP),ip

        END DO
        CALL FORCE_HYBR_MTS_COM()

    END IF

    if (Kmts .ge. mtsparam) then
        mtsparam2 = mtsparam+mtsparam2+TSMTS
      	Kmts = -TSMTS
    end if

    CALL MOVE()

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


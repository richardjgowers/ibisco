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


      SUBROUTINE LOOP_HYBR_MTS3_COM ()
      USE VAR
      USE RNEMD
      IMPLICIT NONE
      INTEGER :: I,IL, JL, HL,IP, HH, TI, KL
      INTEGER :: NOUP = 0
      REAL(KIND=RKIND) :: SUMTOTX, SUMTOTY, SUMTOTZ, K1,K2,K3
      REAL(KIND=RKIND) :: SPX, SPY, SPZ
      REAL(KIND=RKIND) :: RPX,RPY,RPZ

      REAL(KIND=RKIND) :: FXfirst,FXsec,FXthird,FXfourth,FYfirst,FYsec,FYthird, &
                              FYfourth,FZfirst,FZsec,FZthird,FZfourth,fxfifth,fyfifth,fzfifth

!      real(KIND=RKIND) :: ex,ey,ez,moderr,ggg,Fmod
!      real,dimension(10000) :: errFx,errFy,errFz,merr

do I = 1, NSTEP

!write(*,*)'**********',i

    timestepcheck = i
    Kmts = Kmts+1
    K1   = real(Kmts)
    K2   = (real(Kmts)**2)
    K3   = (real(Kmts)**3)

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
              DO JL = 1,VIRT_NUMATOMS(TI)
                 KL = VIRT_ATM_IND(IL,JL)
                 SUMTOTX = SUMTOTX+MASS(ITYPE(KL))*SX(KL)
                 SUMTOTY = SUMTOTY+MASS(ITYPE(KL))*SY(KL)
                 SUMTOTZ = SUMTOTZ+MASS(ITYPE(KL))*SZ(KL)
              END DO
              VIRTRX(IL) = SUMTOTX*VIRT_INVMASS(TI)
              VIRTRY(IL) = SUMTOTY*VIRT_INVMASS(TI)
              VIRTRZ(IL) = SUMTOTZ*VIRT_INVMASS(TI)

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

              DO JL = 1,VIRT_NUMATOMS(TI)
                 KL = VIRT_ATM_IND(IL,JL)
                 SUMTOTX = SUMTOTX+MASS(ITYPE(KL))*SX(KL)
                 SUMTOTY = SUMTOTY+MASS(ITYPE(KL))*SY(KL)
                 SUMTOTZ = SUMTOTZ+MASS(ITYPE(KL))*SZ(KL)
              END DO

              VIRTRX(IL) = SUMTOTX*VIRT_INVMASS(TI)
              VIRTRY(IL) = SUMTOTY*VIRT_INVMASS(TI)
              VIRTRZ(IL) = SUMTOTZ*VIRT_INVMASS(TI)

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
    IF (I .LT. MTSPARAM2) THEN
        CALL FORCE_HYBR_PREMTS_COM ()
        IF(timestepcheck .EQ. I0) THEN
            I0 = I0+MTSPARAM+TSMTS
! write(*,*)i0,'I0'
        ELSE IF (timestepcheck .EQ. I1) THEN
            I1 = I1+MTSPARAM+TSMTS
!write(*,*)i1,'I1'
        ELSE IF (timestepcheck .EQ. I2) THEN
            I2 = I2+MTSPARAM+TSMTS
!write(*,*)i2,'I2'
        ELSE IF (timestepcheck .EQ. I3) THEN 
            I3 = I3+MTSPARAM+TSMTS
!write(*,*)i3,'I3'
        ELSE IF (timestepcheck .EQ. I4) THEN
            I4 = I4+MTSPARAM+TSMTS
!write(*,*)i4,'I4'
        ELSE IF (timestepcheck .EQ. I5) THEN
            I5 = I5+MTSPARAM+TSMTS
!write(*,*)i5,'I5'
        END IF

    ELSE

        DO HH =1,NUM_BEAD
           IP = INDEX_AB(HH)
           CALL VIRT_FORCE_COM_MTS(IP,FCUTB)
! TAYLOR 3: 4 points

            FXfirst = (FXiii(IP)*0.333333333 - 1.5*FXii(IP) + 3*FXi1(IP) - FXi0(IP)*1.833333333)*idt
            FYfirst = (FYiii(IP)*0.333333333 - 1.5*FYii(IP) + 3*FYi1(IP) - FYi0(IP)*1.833333333)*idt
            FYfirst = (FZiii(IP)*0.333333333 - 1.5*FZii(IP) + 3*FZi1(IP) - FZi0(IP)*1.833333333)*idt
            FXsec   = (-FXiii(IP) + 4*FXii(IP) - 5*FXi1(IP) +  2*FXi0(IP))*idt2
            FYsec   = (-FYiii(IP) + 4*FYii(IP) - 5*FYi1(IP) +  2*FYi0(IP))*idt2
            FZsec   = (-FZiii(IP) + 4*FZii(IP) - 5*FZi1(IP) +  2*FZi0(IP))*idt2
            FXthird = (FXiii(IP) - 3*FXii(IP)+ 3*FXi1(IP) - FXi0(IP))*idt3
            FYthird = (FYiii(IP) - 3*FYii(IP)+ 3*FYi1(IP) - FYi0(IP))*idt3
            FZthird = (FZiii(IP) - 3*FZii(IP)+ 3*FZi1(IP) - FZi0(IP))*idt3

! Force calculation with Taylor 3

            FX(IP) = FXiii(IP)+FXfirst*DT*K1 + FXsec*mtsDT2*K2 + FXthird*DT3*K3
            FY(IP) = FYiii(IP)+FYfirst*DT*K1 + FYsec*mtsDT2*K2 + FYthird*DT3*K3
            FZ(IP) = FZiii(IP)+FZfirst*DT*K1 + FZsec*mtsDT2*K2 + FZthird*DT3*K3

        END DO
        CALL FORCE_HYBR_MTS_COM()
    END IF

    if (kmts .ge. mtsparam) then
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



      SUBROUTINE FORCE_HYBR_MTS ()

      USE VAR
      IMPLICIT NONE

      INTEGER :: I, J, K, M, L, hh, h
	  INTEGER :: JBEG, JEND, JNAB, TI, TJ, TIJ, NI
      REAL(KIND=RKIND) :: RCUTSQ, RCUTIBRSQ, ALPHA, FCUT
      REAL(KIND=RKIND) :: RXI, RYI, RZI, FXI, FYI, FZI, FXIJa, FYIJa, FZIJa
      REAL(KIND=RKIND) :: RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
      REAL(KIND=RKIND) :: VIJ, FIJ, RIJ
!      *******************************************************************

VNBOND = 0.0D0
VBOND  = 0.0D0
VANGLE = 0.0D0
VTOR   = 0.0D0
VOOP   = 0.0D0
VFXNB  = 0.0D0
VFYNB  = 0.0D0
VFZNB  = 0.0D0

!	If we choose the ibrid description IBRDESCR = 0. Then the program chek 
!	if the I-th atom is a Bead or not. If yes program use the VIRTSITE 
!	variable to decide to calculate the force on the virtual site using the 
!	real atoms or the Center of Mass


!       ** USE THE LIST TO FIND THE NEIGHBOURS **
	DO 200 hh = NUM_BEAD+1, NUM_BA 
        i = INDEX_AB(hh)

!***************************************************************
        JBEG = POINT(i)
        JEND = POINT(i + 1) - 1

!       ** CHECK THAT BEAD I HAS NEIGHBOURS **
!write(2000,*)hh,INDEX_AB(hh),JBEG,JEND
IF( JBEG .LE. JEND ) THEN

    RXI = RX(I)
    RYI = RY(I)
	RZI = RZ(I)
	FXI = FX(I)
    FYI = FY(I)
    FZI = FZ(I)
	TI = ITYPE(I)

    DO 199 JNAB = JBEG, JEND

!		TAKE THE INDEX OF NEIGHBOUR ATOMS
		J = LIST(JNAB)

	    if(type_label(j) .eq. 2)then
!		TAKE THE TYPE OF NEIGHBOUR ATOM AND NON-BONDED INTERACTIONS	
		    TJ = ITYPE(J)
	   	    TIJ = INBONDT(TI, TJ)
		
 	        IF( TIJ .NE. 0) THEN	
                RXIJ = RXI - RX(J)
           	    RYIJ = RYI - RY(J)
           	    RZIJ = RZI - RZ(J)
           	    RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
           	    RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
           	    RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ
           	    RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0
 		        IF ( RIJSQ < FCUTB ) THEN
		            RIJ = DSQRT(RIJSQ)
                    NI = INT (RIJ / BINNB(TIJ))
                        IF(NI .GT. NDATNB(TIJ)) THEN
                            WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                            WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ
                            WRITE(1,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                            WRITE(*,*)'Simulation stopped at time step: ', timestepcheck
                            STOP
                        END IF
               
!		LINEAR INTEPOLATION

                    ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)
                    FIJ = NBOND_FORCE(TIJ,NI)*(1.0-ALPHA) &
                              + ALPHA*NBOND_FORCE(TIJ,NI+1) 
		            VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
		                      + ALPHA*NBOND_POT(TIJ,NI+1)
    		        VNBOND = VNBOND + VIJ
    		        FXIJ  = FIJ * RXIJ
                    FYIJ  = FIJ * RYIJ
                    FZIJ  = FIJ * RZIJ
   
!     QUI HO TOLTO LA PARTE DOVE C'ERA IL DPD, TANTO PER ORA NON SERVE

                    FXI   = FXI + FXIJ
                    FYI   = FYI + FYIJ
                    FZI   = FZI + FZIJ
                    FX(J) = FX(J) - FXIJ
                    FY(J) = FY(J) - FYIJ
                    FZ(J) = FZ(J) - FZIJ
		
!		ADD THE NON-BONDED PART OF PRESSURE
                    PT11 = PT11 + FXIJ * RXIJ
                    PT22 = PT22 + FYIJ * RYIJ
                    PT33 = PT33 + FZIJ * RZIJ
     	            PT12 = PT12 + FYIJ * RXIJ
                    PT13 = PT13 + FZIJ * RXIJ
                    PT23 = PT23 + FZIJ * RYIJ


                END IF   ! endif  IF ( RIJSQ < RCUTSQ )
            END IF     ! endif   IF( TIJ .NE. 0) THEN
    else
	   	    TIJ = INBONDT(TI, vitype(virtual_center(J)))
!write(3000,*)'BEAD', I,J
 	        IF( TIJ .NE. 0) THEN	
                RXIJ = RXI - RX(J)
           	    RYIJ = RYI - RY(J)
           	    RZIJ = RZI - RZ(J)
           	    RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
           	    RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
           	    RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ
           	    RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0
 		        IF ( RIJSQ < FCUTB ) THEN
		            RIJ = DSQRT(RIJSQ)
                    NI = INT (RIJ / BINNB(TIJ))
                        IF(NI .GT. NDATNB(TIJ)) THEN
                            WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                            WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ
                            WRITE(1,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                            WRITE(*,*)'Simulation stopped at time step: ', timestepcheck
                            STOP
                        END IF
               
!		LINEAR INTEPOLATION
                    ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)
                    FIJ = NBOND_FORCE(TIJ,NI)*(1.0-ALPHA) &
                              + ALPHA*NBOND_FORCE(TIJ,NI+1) 
		            VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
		                      + ALPHA*NBOND_POT(TIJ,NI+1)
		            VNBOND = VNBOND + VIJ
		            FXIJ  = FIJ * RXIJ
                    FYIJ  = FIJ * RYIJ
                    FZIJ  = FIJ * RZIJ

                    FXI   = FXI + FXIJ
                    FYI   = FYI + FYIJ
                    FZI   = FZI + FZIJ

                    DO H = 1,init_numbcomp(virtual_center(J))
!write(3000,*)INDX_ATM(virtual_center(j),H),I,J
			            FXIJa = FXIJ*MASS(ITYPE(INDX_ATM(virtual_center(J),H)))*INVTOTBMASS(virtual_center(J))
			            FYIJa = FYIJ*MASS(ITYPE(INDX_ATM(virtual_center(J),H)))*INVTOTBMASS(virtual_center(J))
			            FZIJa = FZIJ*MASS(ITYPE(INDX_ATM(virtual_center(J),H)))*INVTOTBMASS(virtual_center(J))
                        FX(INDX_ATM(virtual_center(J),H)) = FX(INDX_ATM(virtual_center(J),H)) - FXIJa
                        FY(INDX_ATM(virtual_center(J),H)) = FY(INDX_ATM(virtual_center(J),H)) - FYIJa
                        FZ(INDX_ATM(virtual_center(J),H)) = FZ(INDX_ATM(virtual_center(J),H)) - FZIJa
                    END DO
		
!		ADD THE NON-BONDED PART OF PRESSURE
                    PT11 = PT11 + FXIJ * RXIJ
                    PT22 = PT22 + FYIJ * RYIJ
                    PT33 = PT33 + FZIJ * RZIJ
 	                PT12 = PT12 + FYIJ * RXIJ
                    PT13 = PT13 + FZIJ * RXIJ
                    PT23 = PT23 + FZIJ * RYIJ

                END IF   ! endif  IF ( RIJSQ < RCUTB )
            END IF     ! endif   IF( TIJ .NE. 0) THEN
        end if ! if(type_label(j) .eq. 1)
199     CONTINUE

        FX(I) = FXI
        FY(I) = FYI
        FZ(I) = FZI

    ENDIF
200     CONTINUE

!###############################################################################
!###############################################################################
! Loop through the VS

do 220 hh=num_BA+1,num_vs
    i = INDEX_AB(hh)
    JBEG = POINT(I)
    JEND = POINT(I+1) - 1

!       ** CHECK THAT BEAD I HAS NEIGHBOURS **

IF( JBEG .LE. JEND ) THEN

    RXI = RX(I)
    RYI = RY(I)
	RZI = RZ(I)
	FXI = FX(I)
    FYI = FY(I)
    FZI = FZ(I)
	TI = ITYPE(I)

    DO 219 JNAB = JBEG, JEND

!		TAKE THE INDEX OF NEIGHBOUR ATOMS or BEADS
        J = LIST(JNAB)
        if(type_label(j) .eq. 1)then
	
!		TAKE THE TYPE OF NEIGHBOUR ATOM AND NON-BONDED INTERACTIONS	
		    TJ = ITYPE(J)
	   	    TIJ = INBONDT(TI, TJ)
 	        IF( TIJ .NE. 0) THEN	
                RXIJ = RXI - RX(J)
           	    RYIJ = RYI - RY(J)
           	    RZIJ = RZI - RZ(J)
           	    RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
           	    RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
           	    RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ
           	    RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0
 		        IF ( RIJSQ < FCUTA ) THEN
		            RIJ = DSQRT(RIJSQ)
                    NI = INT (RIJ / BINNB(TIJ))
                        IF(NI .GT. NDATNB(TIJ)) THEN
                            WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                            WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ
                            WRITE(1,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                            WRITE(*,*)'Simulation stopped at time step: ', timestepcheck
                            STOP
                        END IF
               
!		LINEAR INTEPOLATION
                    ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)
                    FIJ = NBOND_FORCE(TIJ,NI)*(1.0-ALPHA) &
                              + ALPHA*NBOND_FORCE(TIJ,NI+1) 
		            VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
		                      + ALPHA*NBOND_POT(TIJ,NI+1)
		            VNBOND = VNBOND + VIJ
		            FXIJ  = FIJ * RXIJ
                    FYIJ  = FIJ * RYIJ
                    FZIJ  = FIJ * RZIJ
   
!     QUI HO TOLTO LA PARTE DOVE C'ERA IL DPD, TANTO PER ORA NON SERVE
                    FXI   = FXI + FXIJ
                    FYI   = FYI + FYIJ
                    FZI   = FZI + FZIJ
                    FX(J) = FX(J) - FXIJ
                    FY(J) = FY(J) - FYIJ
                    FZ(J) = FZ(J) - FZIJ
		
!		ADD THE NON-BONDED PART OF PRESSURE
                    PT11 = PT11 + FXIJ * RXIJ
                    PT22 = PT22 + FYIJ * RYIJ
                    PT33 = PT33 + FZIJ * RZIJ
 	                PT12 = PT12 + FYIJ * RXIJ
                    PT13 = PT13 + FZIJ * RXIJ
                    PT23 = PT23 + FZIJ * RYIJ

            END IF   ! endif  IF ( RIJSQ < RCUTSQ )
        END IF     ! endif   IF( TIJ .NE. 0) THEN

        else ! Here we have the interaction with beads

		    TJ = ITYPE(J)
	   	    TIJ = INBONDT(vitype(virtual_center(I)), TJ)

 	        IF( TIJ .NE. 0) THEN	
                RXIJ = RXI - RX(J)
           	    RYIJ = RYI - RY(J)
           	    RZIJ = RZI - RZ(J)
           	    RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
           	    RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
           	    RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ
           	    RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0
 		        IF ( RIJSQ < FCUTB ) THEN
		            RIJ = DSQRT(RIJSQ)
                    NI = INT (RIJ / BINNB(TIJ))
                        IF(NI .GT. NDATNB(TIJ)) THEN
                            WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                            WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ
                            WRITE(1,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                            WRITE(*,*)'Simulation stopped at time step: ', timestepcheck
                            STOP
                        END IF
               
!		LINEAR INTEPOLATION
                    ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)
                    FIJ = NBOND_FORCE(TIJ,NI)*(1.0-ALPHA) &
                              + ALPHA*NBOND_FORCE(TIJ,NI+1) 
		            VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
		                      + ALPHA*NBOND_POT(TIJ,NI+1)
		            VNBOND = VNBOND + VIJ
		            FXIJ  = FIJ * RXIJ
                    FYIJ  = FIJ * RYIJ
                    FZIJ  = FIJ * RZIJ

                    DO H = 1,init_numbcomp(virtual_center(I))
!write(1000,*)INDX_ATM(virtual_center(I),H),I,J
			            FXIJa = FXIJ*MASS(ITYPE(INDX_ATM(virtual_center(I),H)))*INVTOTBMASS(virtual_center(I))
			            FYIJa = FYIJ*MASS(ITYPE(INDX_ATM(virtual_center(I),H)))*INVTOTBMASS(virtual_center(I))
			            FZIJa = FZIJ*MASS(ITYPE(INDX_ATM(virtual_center(I),H)))*INVTOTBMASS(virtual_center(I))
                        FX(INDX_ATM(virtual_center(I),H)) = FX(INDX_ATM(virtual_center(I),H)) + FXIJa
                        FY(INDX_ATM(virtual_center(I),H)) = FY(INDX_ATM(virtual_center(I),H)) + FYIJa
                        FZ(INDX_ATM(virtual_center(I),H)) = FZ(INDX_ATM(virtual_center(I),H)) + FZIJa
                    END DO
		
!		ADD THE NON-BONDED PART OF PRESSURE
                    PT11 = PT11 + FXIJ * RXIJ
                    PT22 = PT22 + FYIJ * RYIJ
                    PT33 = PT33 + FZIJ * RZIJ
 	                PT12 = PT12 + FYIJ * RXIJ
                    PT13 = PT13 + FZIJ * RXIJ
                    PT23 = PT23 + FZIJ * RYIJ

                END IF   ! endif  IF ( RIJSQ < RCUTSQ )
            END IF     ! endif   IF( TIJ .NE. 0) THEN
        end if ! if(type_label(j) .eq. 1)
219     CONTINUE

            FX(I) = FXI
            FY(I) = FYI
            FZ(I) = FZI

    ENDIF
220     CONTINUE

!###############################################################################
!###############################################################################

! Loop through atoms

do 210 hh=num_vs+1,natoms
    i = INDEX_AB(hh)
    JBEG = POINT(I)
    JEND = POINT(I+1) - 1

!       ** CHECK THAT BEAD I HAS NEIGHBOURS **

IF( JBEG .LE. JEND ) THEN

    RXI = RX(I)
    RYI = RY(I)
	RZI = RZ(I)
	FXI = FX(I)
    FYI = FY(I)
    FZI = FZ(I)
	TI = ITYPE(I)

    DO 209 JNAB = JBEG, JEND

!		TAKE THE INDEX OF NEIGHBOUR ATOMS
		J = LIST(JNAB)
	
!		TAKE THE TYPE OF NEIGHBOUR ATOM AND NON-BONDED INTERACTIONS	
		TJ = ITYPE(J)
	   	TIJ = INBONDT(TI, TJ)
		
 	    IF( TIJ .NE. 0) THEN	
            RXIJ = RXI - RX(J)
           	RYIJ = RYI - RY(J)
           	RZIJ = RZI - RZ(J)
           	RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
           	RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
           	RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ
           	RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0
 		    IF ( RIJSQ < FCUTA ) THEN
		        RIJ = DSQRT(RIJSQ)
                NI = INT (RIJ / BINNB(TIJ))
                    IF(NI .GT. NDATNB(TIJ)) THEN
                        WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                        WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ
                        WRITE(1,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                        WRITE(*,*)'Simulation stopped at time step: ', timestepcheck
                        STOP
                    END IF
               
!		LINEAR INTEPOLATION

                ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)
                FIJ = NBOND_FORCE(TIJ,NI)*(1.0-ALPHA) &
                          + ALPHA*NBOND_FORCE(TIJ,NI+1) 
		        VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
		                  + ALPHA*NBOND_POT(TIJ,NI+1)
		        VNBOND = VNBOND + VIJ
		        FXIJ  = FIJ * RXIJ
                FYIJ  = FIJ * RYIJ
                FZIJ  = FIJ * RZIJ
   
!     QUI HO TOLTO LA PARTE DOVE C'ERA IL DPD, TANTO PER ORA NON SERVE

                FXI   = FXI + FXIJ
                FYI   = FYI + FYIJ
                FZI   = FZI + FZIJ
                FX(J) = FX(J) - FXIJ
                FY(J) = FY(J) - FYIJ
                FZ(J) = FZ(J) - FZIJ
		
!		ADD THE NON-BONDED PART OF PRESSURE
                PT11 = PT11 + FXIJ * RXIJ
                PT22 = PT22 + FYIJ * RYIJ
                PT33 = PT33 + FZIJ * RZIJ
 	            PT12 = PT12 + FYIJ * RXIJ
                PT13 = PT13 + FZIJ * RXIJ
                PT23 = PT23 + FZIJ * RYIJ


            END IF   ! endif  IF ( RIJSQ < RCUTSQ )
        END IF     ! endif   IF( TIJ .NE. 0) THEN
209     CONTINUE

        FX(I) = FXI
        FY(I) = FYI
        FZ(I) = FZI

    ENDIF

210     CONTINUE
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################



!	SAVE THE NON-BONDED PART OF FORCE
	DO hh=NUM_BEAD+1, NATOMS
	    i = INDEX_AB(hh)
	    FXNB(i) = FX(i) + VFXNB(i)
	    FYNB(i) = FY(i) + VFYNB(i)
	    FZNB(i) = FZ(i) + VFZNB(i)
	END DO	
	
    DO 300 hh = NUM_BEAD+1, NATOMS
        i = INDEX_AB(hh)
        RXI = SX(i)
        RYI = SY(i)
        RZI = SZ(i)
        FXI = FX(i)
        FYI = FY(i)
        FZI = FZ(i)
        TI = ITYPE(i)

!	*******************************************************************************************
!	*********************CALCULATE THE BONDED FORCE AND POTENTIAL******************************
      CALL HFPBOND (i, RXI, RYI, RZI, FXI, FYI, FZI, TI)
	
!	*******************************************************************************************
!	*************************CALCULATE THE ANGLE FORCE AND POTENTIAL***************************
      CALL HFPANGLE (i, RXI, RYI, RZI, FXI, FYI, FZI, TI)

!	*******************************************************************************************
!	**********************CALCULATE THE TORSION FORCE AND POTENTIAL****************************
      CALL HFPTOR (i, RXI, RYI, RZI, FXI, FYI, FZI, TI)

!	*******************************************************************************************
!	**********************CALCULATE THE IMPROPER TORSION FORCE AND POTENTIAL********************
      CALL FPOUTPLANE (i, RXI, RYI, RZI, FXI, FYI, FZI, TI)

      FX(i) = FXI
      FY(i) = FYI
      FZ(i) = FZI

300	CONTINUE

      DO I = 1, NATOMS

            PT11 = PT11 + (FX(i) - FXNB(i))*SX(i)
            PT22 = PT22 + (FY(i) - FYNB(i))*SY(i)
            PT33 = PT33 + (FZ(i) - FZNB(i))*SZ(i)

            PT12 = PT12 + (FY(i) - FYNB(i))*SX(i)
            PT13 = PT13 + (FZ(i) - FZNB(i))*SX(i)
            PT23 = PT23 + (FZ(i) - FZNB(i))*SY(i)

      END DO


      RETURN
      END

!	*********************************************************************************************

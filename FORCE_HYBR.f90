
SUBROUTINE FORCE_HYBR ( )

USE VAR
IMPLICIT NONE

INTEGER :: I, J, hh, h
INTEGER :: JBEG, JEND, JNAB, TI, TJ, TIJ, NI
REAL(kind=rkind) :: RCUTSQ, RCUTIBRSQ, ALPHA, FCUT
REAL(kind=rkind) :: RXI, RYI, RZI, FXI, FYI, FZI, FXIJa, FYIJa, FZIJa
REAL(kind=rkind) :: RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
REAL(kind=rkind) :: VIJ, FIJ, RIJ

!      REAL,PARAMETER :: hrij = 0.4

!DO 100 I = 1, NATOMS

!           FX(I) = 0.0D0
!           FY(I) = 0.0D0
!           FZ(I) = 0.0D0

!100     CONTINUE

    VNBOND     = 0.0D0
	VBOND      = 0.0D0
	VANGLE     = 0.0D0
	VTOR       = 0.0D0
	VOOP       = 0.0D0
	VFXNB      = 0.0D0
	VFYNB      = 0.0D0
	VFZNB      = 0.0D0
	VBOND_MIX  = 0.0D0
	VANGLE_MIX = 0.0D0
	VTOR_MIX   = 0.0D0
	VOOP_MIX   = 0.0D0
    VNBOND_MIX = 0.0D0
	VBOND_CG   = 0.0D0
	VANGLE_CG  = 0.0D0
	VTOR_CG    = 0.0D0
	VOOP_CG    = 0.0D0
    VNBOND_CG  = 0.0D0

!       ** USE THE LIST TO FIND THE NEIGHBOURS **

!###############################################################################
! Loop through BEADS
!write(*,*)num_bead,num_vs
do 200 hh=1,num_bead
    i = INDEX_AB(hh)

    JBEG = POINT(I)
    JEND = POINT(I+1) - 1

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
!write(2000,*)'BEAD', I,j,type_label(i),type_label(J)
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

    		        VNBOND_CG = VNBOND_CG + VIJ

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

		            VNBOND_MIX = VNBOND_MIX + VIJ

		            FXIJ  = FIJ * RXIJ
                    FYIJ  = FIJ * RYIJ
                    FZIJ  = FIJ * RZIJ

                    FXI   = FXI + FXIJ
                    FYI   = FYI + FYIJ
                    FZI   = FZI + FZIJ

                    DO H = 1,init_numbcomp(virtual_center(J))
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

do 220 hh=num_bead+1,num_vs
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
!if(j .eq. 1 .and. i .eq. 3727)write(*,*)j,i,RIJSQ,FCUTB
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

		            VNBOND_MIX = VNBOND_MIX + VIJ

		            FXIJ  = FIJ * RXIJ
                    FYIJ  = FIJ * RYIJ
                    FZIJ  = FIJ * RZIJ

                    FX(J) = FX(J) - FXIJ
                    FY(J) = FY(J) - FYIJ
                    FZ(J) = FZ(J) - FZIJ

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
	DO I=1, NATOMS
	
	FXNB(I) = FX(I) + VFXNB(I)
	FYNB(I) = FY(I) + VFYNB(I)
	FZNB(I) = FZ(I) + VFZNB(I)

	END DO	

      DO 300 I = 1, NATOMS

      RXI = SX(I)
      RYI = SY(I)
      RZI = SZ(I)
      FXI = FX(I)
      FYI = FY(I)
      FZI = FZ(I)
	
      TI = ITYPE(I)

!	*******************************************************************************************
!	*********************CALCULATE THE BONDED FORCE AND POTENTIAL******************************
      CALL HFPBOND (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

!	*******************************************************************************************
!	*************************CALCULATE THE ANGLE FORCE AND POTENTIAL***************************
      CALL HFPANGLE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

!	*******************************************************************************************
!	**********************CALCULATE THE TORSION FORCE AND POTENTIAL****************************
      CALL HFPTOR (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

!	*******************************************************************************************
!	**********************CALCULATE THE IMPROPER TORSION FORCE AND POTENTIAL********************
     CALL FPOUTPLANE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

      FX(I) = FXI
      FY(I) = FYI
      FZ(I) = FZI

300	CONTINUE

      DO I = 1, NATOMS

            PT11 = PT11 + (FX(I) - FXNB(I))*SX(I)
            PT22 = PT22 + (FY(I) - FYNB(I))*SY(I)
            PT33 = PT33 + (FZ(I) - FZNB(I))*SZ(I)

            PT12 = PT12 + (FY(I) - FYNB(I))*SX(I)
            PT13 = PT13 + (FZ(I) - FZNB(I))*SX(I)
            PT23 = PT23 + (FZ(I) - FZNB(I))*SY(I)

      END DO


      RETURN
      END

!	*********************************************************************************************


SUBROUTINE FORCE_HYBR_PREMTS ()

USE VAR
IMPLICIT NONE

      INTEGER :: I, J, K, M, L , hh, h, atm, tatm
      INTEGER :: JBEG, JEND, JNAB, TI, TJ, TIJ, NI
      REAL(KIND=RKIND) :: RCUTSQ, RCUTIBRSQ, ALPHA, FCUT
      REAL(KIND=RKIND) :: RXI, RYI, RZI, FXI, FYI, FZI, FXIJa, FYIJa, FZIJa
      REAL(KIND=RKIND) :: RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
      REAL(KIND=RKIND) :: VIJ, FIJ, RIJ

!      *******************************************************************

!	RCUTSQ = RCUT*RCUT
 !	RCUTIBRSQ = RCUTIBR*RCUTIBR

!do I = 1, NATOMS
!   FX(INDEX_AB(I)) = 0.0D0
!   FY(INDEX_AB(I)) = 0.0D0
!   FZ(INDEX_AB(I)) = 0.0D0
!end do

VNBOND = 0.0D0
VBOND  = 0.0D0
VANGLE = 0.0D0
VTOR   = 0.0D0
VOOP   = 0.0D0
VFXNB  = 0.0D0
VFYNB  = 0.0D0
VFZNB  = 0.0D0

!       ** USE THE LIST TO FIND THE NEIGHBOURS **

DO 400 hh = 1,NUM_BA 
    i = INDEX_AB(hh)
JBEG = POINT(i)
JEND = POINT(i + 1) - 1

!** CHECK THAT ATOM I HAS NEIGHBOURS **

IF( JBEG .LE. JEND ) THEN

    RXI = RX(i)
    RYI = RY(i)
    RZI = RZ(i)
    FXI = FX(i)
    FYI = FY(i)
    FZI = FZ(i)
    TI = ITYPE(i)

    DO 399 JNAB = JBEG, JEND

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
! IN this part we have only beads so we can use directly FCUTB
                IF ( RIJSQ < FCUTB ) THEN
            		RIJ = DSQRT(RIJSQ)
		            NI = INT (RIJ / BINNB(TIJ))
            		IF(NI .GT. NDATNB(TIJ)) THEN
			            WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist1'
			            WRITE(1,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'		
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
           K = VIRT_CENTER(J) !K is the VS
           TJ = VITYPE(K)
           TIJ = INBONDT(TI,TJ)
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

                    DO H = 1,VIRT_NUMATOMS(TJ)
                       atm = VIRT_ATM_IND(K,H)
                       tatm = ITYPE(atm)
                       FX(atm) = FX(atm) - FXIJ*VIRT_MASSCOEFF(TJ,tatm)
                       FY(atm) = FY(atm) - FYIJ*VIRT_MASSCOEFF(TJ,tatm)
                       FZ(atm) = FZ(atm) - FZIJ*VIRT_MASSCOEFF(TJ,tatm)
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

        end if  !if(type_label(j) .eq. 2)      
399     CONTINUE

        FX(i) = FXI
        FY(i) = FYI
        FZ(i) = FZI

        ENDIF
400     CONTINUE
	
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
    FXI = 0.0D0
    FYI = 0.0D0
    FZI = 0.0D0
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
           K = VIRT_CENTER(I)
           TI = VITYPE(K)
           TJ = ITYPE(J)
           TIJ = INBONDT(TI,TJ)

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

                    FX(J) = FX(J) - FXIJ
                    FY(J) = FY(J) - FYIJ
                    FZ(J) = FZ(J) - FZIJ

                    DO H = 1,VIRT_NUMATOMS(TI)
                       atm = VIRT_ATM_IND(K,H)
                       tatm = ITYPE(atm)
                       FX(atm) = FX(atm) + FXIJ*VIRT_MASSCOEFF(TI,tatm)
                       FY(atm) = FY(atm) + FYIJ*VIRT_MASSCOEFF(TI,tatm)
                       FZ(atm) = FZ(atm) + FZIJ*VIRT_MASSCOEFF(TI,tatm)
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

        FX(I) = FX(I) + FXI
        FY(I) = FY(I) + FYI
        FZ(I) = FZ(I) + FZI

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
	DO hh=1, NATOMS
	    i = INDEX_AB(hh)
	    FXNB(i) = FX(i) + VFXNB(i)
	    FYNB(i) = FY(i) + VFYNB(i)
	    FZNB(i) = FZ(i) + VFZNB(i)
	END DO	
	
    DO 300 hh = 1, NATOMS
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

    DO hh = 1, num_bead
        PT11 = PT11 + (FX(hh) - FXNB(hh))*SX(hh)
        PT22 = PT22 + (FY(hh) - FYNB(hh))*SY(hh)
        PT33 = PT33 + (FZ(hh) - FZNB(hh))*SZ(hh)
        PT12 = PT12 + (FY(hh) - FYNB(hh))*SX(hh)
        PT13 = PT13 + (FZ(hh) - FZNB(hh))*SX(hh)
        PT23 = PT23 + (FZ(hh) - FZNB(hh))*SY(hh)
        i = INDEX_AB(hh)
            IF(timestepcheck .EQ. I0) THEN
            FXi0(I) = FX(I)
            FYi0(I) = FY(I) 
            FZi0(I) = FZ(I) 
! write(*,*)i0,'I0'
         ELSE IF (timestepcheck .EQ. I1) THEN
            FXi1(I) = FX(I)
            FYi1(I) = FY(I) 
            FZi1(I) = FZ(I) 
!write(*,*)i1,'I1'
        ELSE IF (timestepcheck .EQ. I2) THEN
            FXii(I) = FX(I)
            FYii(I) = FY(I) 
            FZii(I) = FZ(I) 
!write(*,*)i2,'I2'
        ELSE IF (timestepcheck .EQ. I3) THEN
            FXiii(I) = FX(I)
            FYiii(I) = FY(I) 
            FZiii(I) = FZ(I) 
!write(*,*)i3,'I3'
        ELSE IF (timestepcheck .EQ. I4) THEN
            FXiiii(I) = FX(I)
            FYiiii(I) = FY(I) 
            FZiiii(I) = FZ(I) 
!write(*,*)i4,'I4'
        ELSE IF (timestepcheck .EQ. I5) THEN
            FXv(I) = FX(I)
            FYv(I) = FY(I) 
            FZv(I) = FZ(I) 
!write(*,*)i5,'I5'
        END IF

    END DO

do hh = num_bead+1,natoms
        PT11 = PT11 + (FX(hh) - FXNB(hh))*SX(hh)
        PT22 = PT22 + (FY(hh) - FYNB(hh))*SY(hh)
        PT33 = PT33 + (FZ(hh) - FZNB(hh))*SZ(hh)
        PT12 = PT12 + (FY(hh) - FYNB(hh))*SX(hh)
        PT13 = PT13 + (FZ(hh) - FZNB(hh))*SX(hh)
        PT23 = PT23 + (FZ(hh) - FZNB(hh))*SY(hh)
end do


      RETURN
      END

!	*********************************************************************************************

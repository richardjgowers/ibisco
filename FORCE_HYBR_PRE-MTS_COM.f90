
SUBROUTINE FORCE_HYBR_PREMTS_COM ()

USE VAR
IMPLICIT NONE

INTEGER       I, J, K, M, L, hh
REAL*8        RCUTSQ, RCUTIBRSQ, ALPHA, FCUT
REAL*8        RXI, RYI, RZI, FXI, FYI, FZI
REAL*8        RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
REAL*8        VIJ, FIJ, RIJ
INTEGER       JBEG, JEND, JNAB, TI, TJ, TIJ, NI
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

!      FCUT = fcutB
!###############################################################################
! Loop through BEADS

DO 400 hh = 1,NUM_BEAD 
    i = INDEX_AB(hh)
    CALL VIRT_FORCE_COM(I,FCUTB)
    JBEG = POINT(I)
    JEND = POINT(I + 1) - 1

!** CHECK THAT ATOM I HAS NEIGHBOURS **

    IF( JBEG .LE. JEND ) THEN
        RXI = RX(I)
        RYI = RY(I)
        RZI = RZ(I)
        FXI = FX(I)
        FYI = FY(I)
        FZI = FZ(I)
        TI = ITYPE(I)
      DO 399 JNAB = JBEG, JEND
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
399     CONTINUE

        FX(I) = FXI
        FY(I) = FYI
        FZ(I) = FZI

        ENDIF
400     CONTINUE

!###############################################################################
!###############################################################################
! Loop through Beads connected to atoms

!       ** USE THE LIST TO FIND THE NEIGHBOURS **
    DO 600 hh = NUM_BEAD+1, NUM_BA

        I = INDEX_AB(hh)
        CALL VIRT_FORCE_COM(I,FCUTB)

        JBEG = POINT(I)
        JEND = POINT(I + 1) - 1

!       ** CHECK THAT ATOM I HAS NEIGHBOURS **

    IF( JBEG .LE. JEND ) THEN
        RXI = RX(I)
        RYI = RY(I)
        RZI = RZ(I)
        FXI = FX(I)
        FYI = FY(I)
        FZI = FZ(I)
        TI = ITYPE(I)
        DO 599 JNAB = JBEG, JEND

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
            IF ( RIJSQ < FCUTB ) THEN
    		    RIJ = DSQRT(RIJSQ)
    		    NI = INT (RIJ / BINNB(TIJ))
    		    IF(NI .GT. NDATNB(TIJ)) THEN
                   WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                   WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ,NI
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
599     CONTINUE

        FX(I) = FXI
        FY(I) = FYI
        FZ(I) = FZI

        ENDIF
600     CONTINUE

!###############################################################################

        CALL VIRT_FORCE_COM_SEC (FCUTB)

!###############################################################################
! Loop through ATOMS

!       ** USE THE LIST TO FIND THE NEIGHBOURS **
    DO 200 hh = NUM_BA+1, NATOMS
        I = INDEX_AB(hh)

        JBEG = POINT(I)
        JEND = POINT(I + 1) - 1

!       ** CHECK THAT ATOM I HAS NEIGHBOURS **

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
                   WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ,NI
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
199     CONTINUE

        FX(I) = FXI
        FY(I) = FYI
        FZ(I) = FZI

        ENDIF
200     CONTINUE

!###############################################################################
!###############################################################################
! Common part, for beads, beads connected to atoms and atoms

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


      SUBROUTINE FORCE_HYBR_MTS_COM ()

      USE VAR
      IMPLICIT NONE

      INTEGER       I, J, K, M, L, HH
      REAL*8        RCUTSQ, RCUTIBRSQ, ALPHA, FCUT
      REAL*8        RXI, RYI, RZI, FXI, FYI, FZI
      REAL*8        RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
      REAL*8        VIJ, FIJ, RIJ
	INTEGER       JBEG, JEND, JNAB, TI, TJ, TIJ, NI
!      *******************************************************************

VNBOND = 0.0D0
VBOND  = 0.0D0
VANGLE = 0.0D0
VTOR   = 0.0D0
VOOP   = 0.0D0
VFXNB  = 0.0D0
VFYNB  = 0.0D0
VFZNB  = 0.0D0

!       ** USE THE LIST TO FIND THE NEIGHBOURS **
!###############################################################################
! Loop through BEADS connected to atoms

DO 200 HH = NUM_BEAD+1, NUM_BA
    I = INDEX_AB(HH)
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
		        IF ( RIJSQ < FCUTB ) THEN
            		RIJ = DSQRT(RIJSQ)
            		NI = INT (RIJ / BINNB(TIJ))
            		IF(NI .GT. NDATNB(TIJ)) THEN
			            WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
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
199     CONTINUE

        FX(I) = FXI
        FY(I) = FYI
        FZ(I) = FZI

    END IF
200     CONTINUE

!###############################################################################

    CALL VIRT_SEC_FORCE_COM_MTS(FCUTB)

!###############################################################################
! Loop through atoms

DO 300 HH = NUM_BA+1, NATOMS
    I = INDEX_AB(HH)

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
        DO 299 JNAB = JBEG, JEND

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
299     CONTINUE

        FX(I) = FXI
        FY(I) = FYI
        FZ(I) = FZI

    END IF
300     CONTINUE

!###############################################################################
! Common part to atoms and beads connected to atoms


!	SAVE THE NON-BONDED PART OF FORCE
DO HH=NUM_BEAD+1, NATOMS
    I = INDEX_AB(HH)
	FXNB(I) = FX(I) + VFXNB(I)
	FYNB(I) = FY(I) + VFYNB(I)
	FZNB(I) = FZ(I) + VFZNB(I)
END DO	
	
DO 400 HH = NUM_BEAD+1, NATOMS
    I = INDEX_AB(HH)
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

400	CONTINUE

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

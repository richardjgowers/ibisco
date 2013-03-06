
        SUBROUTINE FDR_STEP_ADDITION (I,RXI,RYI,RZI,FXI,FYI,FZI )
	USE VAR
	IMPLICIT	NONE
        INTEGER       I, J, K, M, L, H 
        REAL*8        RCUTSQ, ALPHA
        REAL*8        RXI, RYI, RZI, FXI, FYI, FZI
        REAL*8        RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
        REAL*8        VIJ, FIJ, RIJ
	REAL(KIND=RKIND) :: FXIJa, FYIJa, FZIJa
	INTEGER       JBEG, JEND, JNAB, TI, TJ, TIJ, NI
!      *******************************************************************

	
	RCUTSQ = RCUTDPD**2.d0

!       ** USE THE ADDITIONAL LIST TO CALCULATE THE DPD/LA THERMOSTATS **

        JBEG = DPDPOINT(I)
        JEND = DPDPOINT(I+1) - 1

!       ** CHECK THAT ATOM I HAS NEIGHBOURS **

        IF( JBEG .LE. JEND ) THEN

        DO 199 JNAB = JBEG, JEND

!		TAKE THE INDEX OF NEIGHBOUR ATOMS
		J = DPDLIST(JNAB)
	
!		TAKE THE TYPE OF NEIGHBOUR ATOM AND NON-BONDED INTERACTIONS	
              	RXIJ = RXI - RX(J)
              	RYIJ = RYI - RY(J)
              	RZIJ = RZI - RZ(J)

              	RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
              	RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
              	RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ

              	RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0

             	IF ( RIJSQ < RCUTSQ ) THEN                             !!!!!bbbbbbbbbbbbbbbbbbbbbbbbbbbb!!!!!!!!!!!

                FXIJ = 0.D0
                FYIJ = 0.D0
                FZIJ = 0.D0
                RIJ = SQRT(RIJSQ)
            
                CALL FDR_STEP(I,J,FXIJ,FYIJ,FZIJ,RXIJ,RYIJ,RZIJ,RIJ)

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
		END IF                                             !!!!!!!!bbbbbbbbbbbbbbbbbbbbbbbbb!!!!!!!!!
199     CONTINUE

        ENDIF
	
        RETURN
        END
!	*********************************************************************************************

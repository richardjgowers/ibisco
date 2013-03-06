!
!
!	Modified by Nicodemo
!	Dec 2010
!
!
!
!
        SUBROUTINE LA_ADDITION (I,RXI,RYI,RZI)
	USE VAR
	IMPLICIT	NONE
        INTEGER       I, J, K, M, L 
        REAL*8        RCUTSQ, ALPHA
        REAL*8        RXI, RYI, RZI, FXI, FYI, FZI
        REAL*8        RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
        REAL*8        VIJ, FIJ, RIJ
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

                CALL LOWEAND (RIJSQ,RXIJ,RYIJ,RZIJ,I,J)

		END IF                                             !!!!!!!!bbbbbbbbbbbbbbbbbbbbbbbbb!!!!!!!!!
199     CONTINUE

	END IF
	
        RETURN
        END
!	*********************************************************************************************

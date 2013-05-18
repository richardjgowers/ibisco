   SUBROUTINE VIRT_SEC_FORCE_COM_MTS (FCUT)
	
	USE VAR
	IMPLICIT NONE

        INTEGER :: I,J,K,H, atm, tatm
	INTEGER :: VJBEG, VJEND, VJNAB
	INTEGER       JBEG, JEND, JNAB, TI, TJ, TIJ, NI
        REAL(KIND=RKIND) :: RCUTSQ, ALPHA, FXIJa, FYIJa, FZIJa,FCUT
        REAL(KIND=RKIND) :: RXI, RYI, RZI, FXI, FYI, FZI
        REAL(KIND=RKIND) :: RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
        REAL(KIND=RKIND) :: VIJ, FIJ, RIJ

!***********************************************************************************

DO I = 1, NVIRTA

        VJBEG = VIRT_POINT_SEC(I)
        VJEND = VIRT_POINT_SEC(I+1) - 1

!       ** CHECK THAT ATOM I HAS NEIGHBOURS **

        IF( VJBEG .LE. VJEND ) THEN

      RXI = VIRTRX(I)
      RYI = VIRTRY(I)
      RZI = VIRTRZ(I)
      TI = VITYPE(I)


        DO VJNAB = VJBEG, VJEND

!		TAKE THE INDEX OF NEIGHBOUR ATOMS
           J = VLIST_SEC(VJNAB)
	
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

		IF ( RIJSQ .LE. FCUT ) THEN

		RIJ = DSQRT(RIJSQ)
		NI = INT (RIJ / BINNB(TIJ))

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

	!	VFXNB(J) = VFXNB(J) + FXIJ
	!	VFYNB(J) = VFYNB(J) + FYIJ
	!	VFZNB(J) = VFZNB(J) + FZIJ
	
!		ADD THE NON-BONDED PART OF PRESSURE
            PT11 = PT11 + FXIJ * RXIJ
            PT22 = PT22 + FYIJ * RYIJ
            PT33 = PT33 + FZIJ * RZIJ
 	
            PT12 = PT12 + FYIJ * RXIJ
            PT13 = PT13 + FZIJ * RXIJ
            PT23 = PT23 + FZIJ * RYIJ

    		DO H = 1,VIRT_NUMATOMS(TI)
         atm = VIRT_ATM_IND(I,H)
         tatm = ITYPE(atm)
         FX(atm) = FX(atm) + FXIJ*VIRT_MASSCOEFF(TI,tatm)
         FY(atm) = FY(atm) + FYIJ*VIRT_MASSCOEFF(TI,tatm)
         FZ(atm) = FZ(atm) + FZIJ*VIRT_MASSCOEFF(TI,tatm)
		END DO
	

	END IF
	END IF

END DO

END IF
END DO

RETURN

END

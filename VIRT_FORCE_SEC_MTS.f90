   SUBROUTINE VIRT_FORCE_SEC_MTS (FCUT)
	
	USE VAR
	IMPLICIT NONE

        INTEGER :: I,J,K,H
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

      RXI = RX(INDEX_VSITE(I))
      RYI = RY(INDEX_VSITE(I))
      RZI = RZ(INDEX_VSITE(I))
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

!                FXI   = FXI + FXIJ
!                FYI   = FYI + FYIJ
!                FZI   = FZI + FZIJ    

!            FX(J) = FX(J) - FXIJ
 !           FY(J) = FY(J) - FYIJ
  !          FZ(J) = FZ(J) - FZIJ

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


!*****************************************************************************************************************
! In the case the virtual site is a real atom the force is weighted on the weight of each single atoms in the bead

        DO H = 1,init_numbcomp(j)
			FXIJa = FXIJ*MASS(ITYPE(INDX_ATM(J,H)))*INVTOTBMASS(J)
			FYIJa = FYIJ*MASS(ITYPE(INDX_ATM(J,H)))*INVTOTBMASS(J)
			FZIJa = FZIJ*MASS(ITYPE(INDX_ATM(J,H)))*INVTOTBMASS(J)

	              FX(INDX_ATM(J,H)) = FX(INDX_ATM(J,H)) - FXIJa
 		        FY(INDX_ATM(J,H)) = FY(INDX_ATM(J,H)) - FYIJa
		        FZ(INDX_ATM(J,H)) = FZ(INDX_ATM(J,H)) - FZIJa

!                  RXIJ = RX(J) - RX(INDX_ATM(J,H))
!              	RYIJ = RY(J) - RY(INDX_ATM(J,H))
!              	RZIJ = RZ(J) - RZ(INDX_ATM(J,H))

!              	RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
!              	RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
!              	RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ

!		ADD THE NON-BONDED PART OF PRESSURE
!        	PT11 = PT11 + FXIJa * RXIJ
!            PT22 = PT22 + FYIJa * RYIJ
!        	PT33 = PT33 + FZIJa * RZIJ
 	
!        	PT12 = PT12 + FYIJa * RXIJ
!        	PT13 = PT13 + FZIJa * RXIJ
!          	PT23 = PT23 + FZIJa * RYIJ




!			VFXNB(INDX_ATM(J,H)) = FX(INDX_ATM(J,H))
!			VFYNB(INDX_ATM(J,H)) = FY(INDX_ATM(J,H))
!			VFZNB(INDX_ATM(J,H)) = FZ(INDX_ATM(J,H))

		END DO
	

	END IF
	END IF

END DO

END IF
END DO

RETURN

END

   SUBROUTINE VIRT_FORCE_COM_MTS (I,FCUT)
	
	USE VAR
	IMPLICIT NONE

      INTEGER,INTENT(IN) :: I
      INTEGER :: J,K,H, atm, tatm
      INTEGER :: VJBEG, VJEND, VJNAB
      INTEGER       JBEG, JEND, JNAB, TI, TJ, TIJ, NI
      REAL(KIND=RKIND) :: RCUTSQ, ALPHA, FXIJa, FYIJa, FZIJa, FCUT
      REAL(KIND=RKIND) :: RXI, RYI, RZI, FXI, FYI, FZI
      REAL(KIND=RKIND) :: RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
      REAL(KIND=RKIND) :: VIJ, FIJ, RIJ
      REAL(KIND=RKIND) :: somma = 0, sommay = 0, sommaz = 0

!***********************************************************************************

        VJBEG = VIRT_POINT(I)
        VJEND = VIRT_POINT(I+1) - 1

!       ** CHECK THAT ATOM I HAS NEIGHBOURS **

        IF( VJBEG .LE. VJEND ) THEN

      RXI = RX(I)
      RYI = RY(I)
      RZI = RZ(I)
      FXI = FX(I)
      FYI = FY(I)
      FZI = FZ(I)
      TI = ITYPE(I)

        DO VJNAB = VJBEG, VJEND

!		TAKE THE INDEX OF NEIGHBOUR ATOMS
		J = VLIST(VJNAB)
	
!		TAKE THE TYPE OF NEIGHBOUR ATOM AND NON-BONDED INTERACTIONS	
		TJ = VITYPE(J)
	      	TIJ = INBONDT(TI, TJ)
		
 	        IF( TIJ .NE. 0) THEN	
              	RXIJ = RXI - VIRTRX(J)
              	RYIJ = RYI - VIRTRY(J)
              	RZIJ = RZI - VIRTRZ(J)

              	RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
              	RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
              	RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ

              	RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0

     !        	IF ( RIJSQ < RCUTSQ ) THEN

			IF (RIJSQ.LT.FCUT) THEN
	!	IF ( RIJSQ .LE. RCUTSQ ) THEN

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
 !               FYI   = FYI + FYIJ
  !              FZI   = FZI + FZIJ

	!        FX(I) = FXI
	 !       FY(I) = FYI
	  !      FZ(I) = FZI

		VFXNB(I) = VFXNB(I) + FXIJ
		VFYNB(I) = VFYNB(I) + FYIJ
		VFZNB(I) = VFZNB(I) + FZIJ


!*****************************************************************************************************************
! The force is weighted on the weight of each single atoms in the bead
        DO H = 1,VIRT_NUMATOMS(TJ)
           atm = VIRT_ATM_IND(J,H)
           tatm = ITYPE(atm)
           FX(atm) = FX(atm) - FXIJ*VIRT_MASSCOEFF(TJ,tatm)
           FY(atm) = FY(atm) - FYIJ*VIRT_MASSCOEFF(TJ,tatm)
           FZ(atm) = FZ(atm) - FZIJ*VIRT_MASSCOEFF(TJ,tatm)
        END DO


	END IF
	END IF

END DO

END IF

     RETURN
        END

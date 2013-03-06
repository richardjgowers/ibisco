
        SUBROUTINE FORCE2 ( )
	USE VAR
	IMPLICIT	NONE
        INTEGER       I, J, K, M, L 
        REAL*8        RCUTSQ, ALPHA
        REAL*8        RXI, RYI, RZI, FXI, FYI, FZI
        REAL*8        RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
        REAL*8        VIJ, FIJ, RIJ, stressxz
	INTEGER       JBEG, JEND, JNAB, TI, TJ, TIJ, NI,IZI,IZJ
!      *******************************************************************

        do IZI = 1,SLIDE
           PZX(IZI) = 0.D0
        ENDDO
	
	RCUTSQ = RCUT*RCUT

        DO 100 I = 1, NATOMS

           FX(I) = 0.0D0
           FY(I) = 0.0D0
           FZ(I) = 0.0D0

100     CONTINUE

        VNBOND = 0.0D0
	VBOND  = 0.0D0
	VANGLE = 0.0D0
	VTOR   = 0.0D0
	VOOP   = 0.0D0

!       ** USE THE LIST TO FIND THE NEIGHBOURS **
	DO 200 I = 1, NATOMS 

        JBEG = POINT(I)
        JEND = POINT(I+1) - 1

!       ** CHECK THAT ATOM I HAS NEIGHBOURS **

        IF( JBEG .LE. JEND ) THEN

        RXI = RX(I)
        RYI = RY(I)
	RZI = RZ(I)
	FXI = FX(I)
        FYI = FY(I)
        FZI = FZ(I)
	TI = ITYPE(I)

        IZI = INT((RZI*BOXZINV+0.5D0)*DBLE(SLIDE)) 
        IZI = 1 + MOD((IZI+SLIDE),SLIDE)
     
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

             	IF ( RIJSQ < RCUTSQ ) THEN

		RIJ = SQRT(RIJSQ)
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

                IF (DPDINPUT.EQ.1.and.RIJ.LT.RCUTDPD) THEN
                   IF (WDPDT.EQ.2) THEN
                    call FDR(I,J,FXIJ,FYIJ,FZIJ,RXIJ,RYIJ,RZIJ,RIJ)
                   ELSEIF (WDPDT.EQ.1) THEN
                    call FDR_STEP(I,J,FXIJ,FYIJ,FZIJ,RXIJ,RYIJ,RZIJ,RIJ)
                   ENDIF
                ENDIF


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
 	
                stressxz = FZIJ*RXIJ            

        	PT12 = PT12 + FYIJ * RXIJ
!       	PT13 = PT13 + FZIJ * RXIJ
                PT13 = PT13 + stressxz
           	PT23 = PT23 + FZIJ * RYIJ

!               PZX(IZI) = PZX(IZI) + FZIJ*RXIJ       
                PZX(IZI) = PZX(IZI) + stressxz     
  
                IF (LAINPUT.EQ.1.AND.RIJ.LT.RCUTDPD) THEN
                call LOWEAND (RIJSQ,RXIJ,RYIJ,RZIJ,I,J)
                ENDIF

                END IF
		END IF   
199     CONTINUE
        
        IF (DPD_BONDED.EQ.1) THEN
           IF (LAINPUT.EQ.1) THEN
           call  LA_ADDITION (I,RXI,RYI,RZI)
           ELSEIF (DPDINPUT.EQ.1) THEN
            IF(WDPDT.EQ.2) THEN
             call FDR_ADDITION (I,RXI,RYI,RZI,FXI,FYI,FZI )
            ELSEIF(WDPDT.EQ.1) THEN
             call FDR_STEP_ADDITION (I,RXI,RYI,RZI,FXI,FYI,FZI )
            ENDIF
           ENDIF
        ENDIF

        FX(I) = FXI
        FY(I) = FYI
        FZ(I) = FZI

        ENDIF
200     CONTINUE
	
!	SAVE THE NON-BONDED PART OF FORCE
	DO I=1, NATOMS
	
	FXNB(I) = FX(I)
	FYNB(I) = FY(I)
	FZNB(I) = FZ(I)
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
	CALL FPBOND (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)
	
!	*******************************************************************************************
!	*************************CALCULATE THE ANGLE FORCE AND POTENTIAL***************************
	CALL FPANGLE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

!	*******************************************************************************************
!	**********************CALCULATE THE TORSION FORCE AND POTENTIAL****************************
	CALL FPTOR (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

!	*******************************************************************************************
!	**********************CALCULATE THE IMPROPERTORSION FORCE AND POTENTIAL********************
	CALL FPOUTPLANE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

        FX(I) = FXI
        FY(I) = FYI
        FZ(I) = FZI

300	CONTINUE

        DO I = 1, NATOMS                      !BONDED PART

                RZI = RZ(I)
                IZI = INT((RZI*BOXZINV+0.5D0)*DBLE(SLIDE))
                IZI = 1 + MOD((IZI+SLIDE),SLIDE)                       

                PT11 = PT11 + (FX(I) - FXNB(I))*SX(I)
                PT22 = PT22 + (FY(I) - FYNB(I))*SY(I)
                PT33 = PT33 + (FZ(I) - FZNB(I))*SZ(I)

                stressxz = (FZ(I) - FZNB(I))*SX(I)        
    
                PT12 = PT12 + (FY(I) - FYNB(I))*SX(I)
!               PT13 = PT13 + (FZ(I) - FZNB(I))*SX(I)
                PT13 = PT13 +  stressxz
                PT23 = PT23 + (FZ(I) - FZNB(I))*SY(I)
               
!               PZX(IZI) = PZX(IZI) + (FZ(I) - FZNB(I))*RXI
                PZX(IZI) = PZX(IZI) + stressxz
        END DO

        RETURN
        END

!	*********************************************************************************************

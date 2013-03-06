
SUBROUTINE FORCE_HYBR_COM ( )

USE VAR
IMPLICIT NONE

INTEGER :: I, J, hh
INTEGER :: JBEG, JEND, JNAB, TI, TJ, TIJ, NI
REAL(kind=rkind) :: RCUTSQ, RCUTIBRSQ, ALPHA, FCUT
REAL(kind=rkind) :: RXI, RYI, RZI, FXI, FYI, FZI
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

do 200 hh=1,num_bead
    i = INDEX_AB(hh)

	CALL VIRT_FORCE_COM(I,FCUTB)

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
 
!***************LOWE ANDERSON***************************************

                IF (LAINPUT.EQ.1.AND.RIJ.LT.RCUTDPD) THEN
			        CALL LOWEAND (RIJSQ,RXIJ,RYIJ,RZIJ,I,J)
                ENDIF

!***************LOWE ANDERSON***************************************


            END IF   ! endif  IF ( RIJSQ < RCUTSQ )
        END IF     ! endif   IF( TIJ .NE. 0) THEN
199     CONTINUE

! DPD_INPUT is the keyword to activate the DPD (or LA) thermostat on the particles (atoms or beads). The thermostat acts only (ricontrollare bene questa affermazione) within the particles that are connected with bond, bend or torsion
        
	IF (DPD_BONDED.EQ.1) THEN
		IF (LAINPUT.EQ.1) THEN
			CALL  LA_ADDITION (I,RXI,RYI,RZI)
		ELSEIF (DPDINPUT.EQ.1) THEN
           		IF(WDPDT.EQ.2) THEN
				CALL FDR_ADDITION (I,RXI,RYI,RZI,FXI,FYI,FZI )
			ELSE IF(WDPDT.EQ.1) THEN
				CALL FDR_STEP_ADDITION (I,RXI,RYI,RZI,FXI,FYI,FZI )
			END IF
        ENDIF
    ENDIF

        FX(I) = FXI
        FY(I) = FYI
        FZ(I) = FZI

    ENDIF
200     CONTINUE
!###############################################################################
!###############################################################################
! Loop through ATOMS

do 210 hh=num_bead+1,natoms
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
 
!***************LOWE ANDERSON***************************************

                IF (LAINPUT.EQ.1.AND.RIJ.LT.RCUTDPD) THEN
			        CALL LOWEAND (RIJSQ,RXIJ,RYIJ,RZIJ,I,J)
                ENDIF

!***************LOWE ANDERSON***************************************


            END IF   ! endif  IF ( RIJSQ < RCUTSQ )
        END IF     ! endif   IF( TIJ .NE. 0) THEN
209     CONTINUE

! DPD_INPUT is the keyword to activate the DPD (or LA) thermostat on the particles (atoms or beads). The thermostat acts only (ricontrollare bene questa affermazione) within the particles that are connected with bond, bend or torsion
        
	IF (DPD_BONDED.EQ.1) THEN
		IF (LAINPUT.EQ.1) THEN
			CALL  LA_ADDITION (I,RXI,RYI,RZI)
		ELSEIF (DPDINPUT.EQ.1) THEN
           		IF(WDPDT.EQ.2) THEN
				CALL FDR_ADDITION (I,RXI,RYI,RZI,FXI,FYI,FZI )
			ELSE IF(WDPDT.EQ.1) THEN
				CALL FDR_STEP_ADDITION (I,RXI,RYI,RZI,FXI,FYI,FZI )
			END IF
        ENDIF
    ENDIF

        FX(I) = FXI
        FY(I) = FYI
        FZ(I) = FZI

    ENDIF
210     CONTINUE

!###############################################################################

    CALL VIRT_FORCE_COM_SEC (FCUTB)

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

  !Inline method of calculating bonded forces
  
  DO A=1,NBONDS(I)

     J = JBOND(I,A)

     !            IF ANY OF NEIGHBOUR OF ITH SITE HAS INDEX J < I THEN WE 
     !            SKIP THE BONDED INTERATION OF THIS GROUP AS IT HAS BEEN CALCULATED BEFORE 
     IF (J > I) THEN
        TJ = ITYPE(J)
        TIJ = IBONDT(TI, TJ)
        RXIJ = RXI - SX(J)
        RYIJ = RYI - SY(J)
        RZIJ = RZI - SZ(J)

        RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0

        RIJ = SQRT(RIJSQ)
        NI = INT((RIJ - RBOND(TIJ,0)) / BINB(TIJ))

        !If doesn't fall within table, use last/first entry
        IF(NI .GT. NDATB(TIJ)) THEN
!           WRITE(*,*)'FATAL ERROR: Entry in bond table', TIJ,' does not exist'
!	   WRITE(*,*) I, J, RIJ
!           WRITE(1,*)'FATAL ERROR: Entry in bond table', TIJ,' does not exist'                   
!           WRITE(*,*)'Simulation stopped at Time Step: ',STEP
!           STOP 'Failure in BONDS.inc'
           WRITE(*,*) 'Overextended bond between atoms ', I, ' and ', J
           WRITE(*,*) 'Distance ', RIJ, ' Timestep ', STEP
           WRITE(*,*) 'Continuing on, but this might be a problem!'
           WRITE(1,*) 'Overextended bond between atoms ', I, ' and ', J
           WRITE(1,*) 'Distance ', RIJ, ' Timestep ', STEP
           WRITE(1,*) 'Continuing on, but this might be a problem!'
           NI = NDATB(TIJ)
        ELSE IF(NI .LT. 0) THEN
           NI = 0
        END IF

        !            LINEAR INTEPOLATION

        ALPHA=(RIJ-RBOND(TIJ,NI))/BINB(TIJ)

        FIJ = BOND_FORCE(TIJ,NI)*(1.0D0-ALPHA) &
             + ALPHA*BOND_FORCE(TIJ,NI+1) 

        VIJ = BOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
             + ALPHA*BOND_POT(TIJ,NI+1) 

        VBOND = VBOND + VIJ

        FXIJ  = FIJ * RXIJ
        FYIJ  = FIJ * RYIJ
        FZIJ  = FIJ * RZIJ

        FXI   = FXI + FXIJ
        FYI   = FYI + FYIJ
        FZI   = FZI + FZIJ

        FXL(J) = FXL(J) - FXIJ
        FYL(J) = FYL(J) - FYIJ
        FZL(J) = FZL(J) - FZIJ

     END IF
  END DO

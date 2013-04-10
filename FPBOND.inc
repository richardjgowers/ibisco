  !inline version of FPBOND.f90
  
  DO L = 1, NBONDS(I) !299

     !     TAKE THE INDEX OF BONDED ATOM

     J = JBOND(I, L)

     !            IF ANY OF NEIGHBOUR OF ITH SITE HAS INDEX J < I THEN WE 
     !            SKIP THE BONDED INTERATION OF THIS GROUP AS IT HAS BEEN CALCULATED BEFORE
     IF (J > I) THEN
        TJ = ITYPE(J)
        TIJ = IBONDT(TI, TJ)
        if(.not. typeBond(tij))then
           RXIJ = RXI - SX(J)
           RYIJ = RYI - SY(J)
           RZIJ = RZI - SZ(J)

           RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0

           RIJ = SQRT(RIJSQ)
           NI = INT (RIJ / BINB(TIJ))

           IF(NI .GT. NDATB(TIJ)) THEN
              WRITE(*,*)'FATAL ERROR: Entry in bond table', TIJ,' does not exist'
              WRITE(1,*)'FATAL ERROR: Entry in bond table', TIJ,' does not exist'     
              write(*,*)i,j,rij                 
              WRITE(*,*)'Simulation stopped at Time Step: ',timestepcheck
              !                        call config()    
              STOP
           END IF

           !            LINEAR INTEPOLATION
           IF (RIJ == RCUT) THEN
              FIJ = BOND_FORCE(TIJ,NI)
              VIJ = BOND_POT(TIJ,NI)
           ELSE
              ALPHA=(RIJ-RBOND(TIJ,NI))/BINB(TIJ)
              FIJ = BOND_FORCE(TIJ,NI)*(1.0D0-ALPHA) &
                   + ALPHA*BOND_FORCE(TIJ,NI+1) 

              VIJ = BOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
                   + ALPHA*BOND_POT(TIJ,NI+1) 
           END IF

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
     end if
  END DO !299      CONTINUE
!      *********************************************************************************************

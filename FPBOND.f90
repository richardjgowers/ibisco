!      
!      Modified by Nicodemo
!      Jan. 2011
!
!      Added: lines 42-46

      SUBROUTINE FPBOND (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

      USE VAR

      IMPLICIT NONE

      INTEGER       I, J, L 
      REAL*8        ALPHA
      REAL*8        RXI, RYI, RZI, FXI, FYI, FZI
      REAL*8        RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
      REAL*8        VIJ, FIJ, RIJ
      INTEGER       TI, TJ, TIJ, NI

!     INTEGER       K, M
!     REAL*8        WIJ
!      *******************************************************************

DO 299 L = 1, NBONDS(I)

!            TAKE THE INDEX OF BONDED ATOM

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
                        call config()    
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

                  FX(J) = FX(J) - FXIJ
                  FY(J) = FY(J) - FYIJ
                  FZ(J) = FZ(J) - FZIJ

            END IF
        end if
299      CONTINUE
      
        RETURN
        END
!      *********************************************************************************************

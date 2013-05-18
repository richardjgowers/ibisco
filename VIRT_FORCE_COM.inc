!inline version of VIRT_FORCE_COM

VJBEG = VIRT_POINT(I)
VJEND = VIRT_POINT(I+1) - 1
!       ** CHECK THAT ATOM I HAS NEIGHBOURS **

      IF( VJBEG .LE. VJEND ) THEN

         RXI = RX(I)
         RYI = RY(I)
         RZI = RZ(I)

         FXI=0.0D0
         FYI=0.0D0
         FZI=0.0D0

         TI = ITYPE(I)

         DO VJNAB = VJBEG, VJEND

!     TAKE THE INDEX OF NEIGHBOUR VIRTUAL SITE
            J = VLIST(VJNAB)	
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

               IF (RIJSQ.LT.FCUT) THEN

                  RIJ = DSQRT(RIJSQ)
                  NI = INT (RIJ / BINNB(TIJ))
                  ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)

                  FIJ = NBOND_FORCE(TIJ,NI)*(1.0D0-ALPHA) &
                  + ALPHA*NBOND_FORCE(TIJ,NI+1) 
                  VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
                  + ALPHA*NBOND_POT(TIJ,NI+1)

                  VNBOND_MIX = VNBOND_MIX + VIJ

                  FXIJ  = FIJ * RXIJ
                  FYIJ  = FIJ * RYIJ
                  FZIJ  = FIJ * RZIJ

                  FXI   = FXI + FXIJ
                  FYI   = FYI + FYIJ
                  FZI   = FZI + FZIJ

                  VFXNB(I) = VFXNB(I) + FXIJ
                  VFYNB(I) = VFYNB(I) + FYIJ
                  VFZNB(I) = VFZNB(I) + FZIJ

                  FXVL(J) = FXVL(J) - FXIJ
                  FYVL(J) = FYVL(J) - FYIJ
                  FZVL(J) = FZVL(J) - FZIJ

               END IF !If RIJ lt RCUT
	    END IF !If TIJ ne 0


         END DO !Do over JBEG JEND

         FXL(I) = FXL(I) + FXI
         FYL(I) = FYL(I) + FYI
         FZL(I) = FZL(I) + FZI

      END IF

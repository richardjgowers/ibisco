!Nonbonded pair include
DO I=1,NITEMS

   JBEG = POINT(I)
   JEND = POINT(I+1) - 1
   IF(JBEG .LE. JEND) THEN
      RXI = RX(I)
      RYI = RY(I)
      RZI = RZ(I)

      FXI = 0.0
      FYI = 0.0
      FZI = 0.0

      TI = ITYPE(I)
      DO JNAB = JBEG, JEND

         J = LIST(JNAB)

         TJ = ITYPE(J)
         
         TIJ = INBONDT(TI,TJ)
         IF(TIJ .NE. 0) THEN
            RXIJ = RXI - RX(J)
            RYIJ = RYI - RY(J)
            RZIJ = RZI - RZ(J)

            RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
            RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
            RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ

            RIJSQ = RXIJ**2.0 + RYIJ**2.0 + RZIJ**2.0
            IF(RIJSQ .LT. RCUTSQ) THEN
               RIJ = DSQRT(RIJSQ)

               NI = INT(RIJ/BINNB(TIJ))

               ALPHA = (RIJ - RNBOND(TIJ,NI))/BINNB(TIJ)

               FIJ = NBOND_FORCE(TIJ,NI)*(1.0D0 - ALPHA) &
                    + NBOND_FORCE(TIJ,NI+1)*ALPHA

               VIJ = NBOND_POT(TIJ,NI)*(1.0D0 - ALPHA) &
                    + NBOND_POT(TIJ,NI+1)*ALPHA

               VNBOND = VNBOND + VIJ

               FXIJ = FIJ * RXIJ
               FYIJ = FIJ * RYIJ
               FZIJ = FIJ * RZIJ

               FXI = FXI + FXIJ
               FYI = FYI + FYIJ
               FZI = FZI + FZIJ

               FXL(J) = FXL(J) - FXIJ
               FYL(J) = FYL(J) - FYIJ
               FZL(J) = FZL(J) - FZIJ

               PT11 = PT11 + FXIJ * RXIJ
               PT22 = PT22 + FYIJ * RYIJ
               PT33 = PT33 + FZIJ * RZIJ

               PT12 = PT12 + FYIJ * RXIJ
               PT13 = PT13 + FZIJ * RXIJ
               PT23 = PT23 + FZIJ * RYIJ

            END IF !End if RIJSQ lt RCUTSQ
         END IF !End if TIJ ne 0
      END DO !End loop over neighbours

      FXL(I) = FXL(I) + FXI
      FYL(I) = FYL(I) + FYI
      FZL(I) = FZL(I) + FZI

   END IF !End if has neighbours
END DO !End loop over all items

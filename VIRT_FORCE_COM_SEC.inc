!inline version of VIRT_SEC_FORCE_COM.f90
      
      VJBEG = VIRT_POINT_SEC(I)
      VJEND = VIRT_POINT_SEC(I+1) - 1

!     ** CHECK THAT ATOM I HAS NEIGHBOURS **

      IF( VJBEG .LE. VJEND ) THEN
         RXI = VIRTRX(I)
         RYI = VIRTRY(I)
         RZI = VIRTRZ(I)
         TI  = VITYPE(I)

         DO VJNAB = VJBEG, VJEND
!     TAKE THE INDEX OF NEIGHBOUR ATOMS
            J = VLIST_SEC(VJNAB)
!     TAKE THE TYPE OF NEIGHBOUR ATOM AND NON-BONDED INTERACTIONS	
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
!     LINEAR INTEPOLATION
                  ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)

                  FIJ = NBOND_FORCE(TIJ,NI)*(1.0D0-ALPHA) &
                  + ALPHA*NBOND_FORCE(TIJ,NI+1) 

                  VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
                  + ALPHA*NBOND_POT(TIJ,NI+1)

                  VNBOND_MIX = VNBOND_MIX + VIJ

                  FXIJ  = FIJ * RXIJ
                  FYIJ  = FIJ * RYIJ
                  FZIJ  = FIJ * RZIJ

                  FXL(J) = FXL(J) - FXIJ
                  FYL(J) = FYL(J) - FYIJ
                  FZL(J) = FZL(J) - FZIJ

                  VFXNB(J) = VFXNB(J) + FXIJ !FIX THIS
                  VFYNB(J) = VFYNB(J) + FYIJ !Shouldn't this be -FXIJ?
                  VFZNB(J) = VFZNB(J) + FZIJ
                  
!     ADD THE NON-BONDED PART OF PRESSURE
                  PT11 = PT11 + FXIJ * RXIJ
                  PT22 = PT22 + FYIJ * RYIJ
                  PT33 = PT33 + FZIJ * RZIJ
                  
                  PT12 = PT12 + FYIJ * RXIJ
                  PT13 = PT13 + FZIJ * RXIJ
                  PT23 = PT23 + FZIJ * RYIJ

                  FXVL(I) = FXVL(I) + FXIJ
                  FYVL(I) = FYVL(I) + FYIJ
                  FZVL(I) = FZVL(I) + FZIJ   

               END IF
            END IF

         END DO

      END IF


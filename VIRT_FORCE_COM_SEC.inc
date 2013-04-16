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
                  FIJ = NBOND_FORCE(TIJ,NI)*(1.0-ALPHA) &
                  + ALPHA*NBOND_FORCE(TIJ,NI+1) 
                  VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
                  + ALPHA*NBOND_POT(TIJ,NI+1)
                  VNBOND_MIX = VNBOND_MIX + VIJ

                  FXIJ  = FIJ * RXIJ
                  FYIJ  = FIJ * RYIJ
                  FZIJ  = FIJ * RZIJ
                  FX(J) = FX(J) - FXIJ
                  FY(J) = FY(J) - FYIJ
                  FZ(J) = FZ(J) - FZIJ
                  VFXNB(J) = VFXNB(J) + FXIJ
                  VFYNB(J) = VFYNB(J) + FYIJ
                  VFZNB(J) = VFZNB(J) + FZIJ
                  
!     ADD THE NON-BONDED PART OF PRESSURE
                  PT11 = PT11 + FXIJ * RXIJ
                  PT22 = PT22 + FYIJ * RYIJ
                  PT33 = PT33 + FZIJ * RZIJ
                  
                  PT12 = PT12 + FYIJ * RXIJ
                  PT13 = PT13 + FZIJ * RXIJ
                  PT23 = PT23 + FZIJ * RYIJ

                  DO H = 1,init_numbcomp(i)
                     K = COMPCOM(I,H)
!     FXIJa = FXIJ*MASS(ITYPE(COMPCOM(I,H)))*INVTOTBMASS(I)
!     FYIJa = FYIJ*MASS(ITYPE(COMPCOM(I,H)))*INVTOTBMASS(I)
!     FZIJa = FZIJ*MASS(ITYPE(COMPCOM(I,H)))*INVTOTBMASS(I)
!     FX(COMPCOM(I,H)) = FX(COMPCOM(I,H)) + FXIJa
!     FY(COMPCOM(I,H)) = FY(COMPCOM(I,H)) + FYIJa
!     FZ(COMPCOM(I,H)) = FZ(COMPCOM(I,H)) + FZIJa
                     FXL(K) = FXL(K) + FXIJ*masscoeff(I,H)
                     FYL(K) = FYL(K) + FYIJ*masscoeff(I,H)
                     FZL(K) = FZL(K) + FZIJ*masscoeff(I,H)
!     VFXNB(COMPCOM(J,H)) = FX(COMPCOM(J,H))
!     VFYNB(COMPCOM(J,H)) = FY(COMPCOM(J,H))
!     VFZNB(COMPCOM(J,H)) = FZ(COMPCOM(J,H))


!     RXIJ = RX(J) - RX(COMPCOM(I,H))
!     RYIJ = RY(J) - RY(COMPCOM(I,H))
!     RZIJ = RZ(J) - RZ(COMPCOM(I,H))

!     RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
!     RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
!     RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ

!     ADD THE NON-BONDED PART OF PRESSURE
!     PT11 = PT11 + FXIJa * RXIJ
!     PT22 = PT22 + FYIJa * RYIJ
!     PT33 = PT33 + FZIJa * RZIJ
                     
!     PT12 = PT12 + FYIJa * RXIJ
!     PT13 = PT13 + FZIJa * RXIJ
!     PT23 = PT23 + FZIJa * RYIJ


                  END DO

                  

               END IF
            END IF

         END DO

      END IF


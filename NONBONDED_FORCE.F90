#include "ibi-preprocess.h"

SUBROUTINE NONBONDED_FORCE(N,INDEX_LIST,POINT,MAXNAB,LIST,VNBOND,RCUT,RCUTSQ)

  USE VAR

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N, INDEX_LIST(N), POINT(N+1), MAXNAB, LIST(MAXNAB)
  REAL(KIND=RKIND), INTENT(INOUT) :: VNBOND
  REAL(KIND=RKIND), INTENT(IN) :: RCUT, RCUTSQ
  INTEGER :: A, I, J, TI, TJ, TIJ
  INTEGER :: JNAB, JBEG, JEND
  REAL(KIND=RKIND) :: FXI,FYI,FZI, FXIJ, FYIJ, FZIJ
  REAL(KIND=RKIND) :: RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,RIJSQ,RIJ
  REAL(KIND=RKIND) :: ALPHA, FIJ, VIJ
  INTEGER :: NI

  DO A=1,N
     I = INDEX_LIST(A) !I is the index of atom being considered
     JBEG = POINT(A)
     JEND = POINT(A+1) - 1

     IF(JBEG .LE. JEND) THEN
        RXI = RX(I)
        RYI = RY(I)
        RZI = RZ(I)

        FXI=0.0D0
        FYI=0.0D0
        FZI=0.0D0
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
#ifdef DEBUG_DETAILEDNB
                 IF(TYPE_LABEL(I) .eq. 1) THEN !ATOM
                    WRITE(7100,*) I, J, VIJ, RIJ, FIJ, ALPHA, NI
                 ELSE IF(TYPE_LABEL(I) .eq. 2) THEN !BEAD/VS
                    IF(I .le. NATOMS) THEN !I BEAD
                       IF(J .le. NATOMS) THEN !J BEAD
                          WRITE(7200,*) I, J, VIJ, RIJ, FIJ
                       ELSE ! J VS
                          WRITE(7300,*) I, J, VIJ, RIJ, FIJ
                       END IF
                    ELSE !I VS
                       WRITE(7300,*) I, J, VIJ, RIJ, FIJ
                    END IF
                 END IF
#endif


                 VNBOND = VNBOND + VIJ

                 FXIJ = FIJ * RXIJ
                 FYIJ = FIJ * RYIJ
                 FZIJ = FIJ * RZIJ

                 FXI = FXI + FXIJ
                 FYI = FYI + FYIJ
                 FZI = FZI + FZIJ

                 FX(J) = FX(J) - FXIJ
                 FY(J) = FY(J) - FYIJ
                 FZ(J) = FZ(J) - FZIJ

                 PT11 = PT11 + FXIJ * RXIJ
                 PT22 = PT22 + FYIJ * RYIJ
                 PT33 = PT33 + FZIJ * RZIJ

                 PT12 = PT12 + FYIJ * RXIJ
                 PT13 = PT13 + FZIJ * RXIJ
                 PT23 = PT23 + FZIJ * RYIJ

              END IF !End if RIJSQ lt RCUTSQ
           END IF !End if TIJ ne 0
        END DO !End loop over neighbours
     END IF
     FX(I) = FX(I) + FXI
     FY(I) = FY(I) + FYI
     FZ(I) = FZ(I) + FZI
  END DO !End loop over all items

  RETURN

END SUBROUTINE NONBONDED_FORCE

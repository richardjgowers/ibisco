#include "ibi-preprocess.h"

SUBROUTINE NONBONDED_FORCE(N,INDEX_LIST,POINT,MAXNAB,LIST,RCUT,RCUTSQ)

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N, INDEX_LIST(N), POINT(N+1), MAXNAB, LIST(MAXNAB)
  REAL(KIND=RKIND), INTENT(IN) :: RCUT, RCUTSQ
  INTEGER :: A, I, J, TI, TJ, TIJ
  INTEGER :: JNAB, JBEG, JEND
  INTEGER :: NI
  REAL(KIND=RKIND) :: FXI,FYI,FZI, FXIJ, FYIJ, FZIJ
  REAL(KIND=RKIND) :: RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,RIJSQ,RIJ
  REAL(KIND=RKIND) :: ALPHA, FIJ, VIJ
  REAL(KIND=RKIND), DIMENSION(NITEMS) :: FXL, FYL, FZL !Local arrays for parallelisation

!$OMP PARALLEL DEFAULT(NONE)&
!$OMP& SHARED(N,NITEMS,INDEX_LIST,POINT,LIST,ITYPE,INBONDT)&
!$OMP& SHARED(RX,RY,RZ,BOXXINV,BOXYINV,BOXZINV,BOXX,BOXY,BOXZ,RCUTSQ)&
!$OMP& SHARED(BINNB,RNBOND,NBOND_FORCE,NBOND_POT,TYPE_LABEL)&
!$OMP& PRIVATE(A,I,TI,J,TJ,TIJ,JBEG,JEND,JNAB)&
!$OMP& PRIVATE(FXI,FYI,FZI,RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,RIJSQ,RIJ)&
!$OMP& PRIVATE(NI,ALPHA,FIJ,VIJ,FXIJ,FYIJ,FZIJ)&
!$OMP& PRIVATE(FXL,FYL,FZL)&
!$OMP& REDUCTION(+:VNBOND_TOTAL,VNBOND_ATOM,VNBOND_BEAD,VNBOND_MIX)&
!$OMP& REDUCTION(+:PT11,PT22,PT33,PT12,PT13,PT23) &
!$OMP& REDUCTION(+:FX,FY,FZ)

DO I=1,NITEMS
   FXL(I) = 0.0D0
   FYL(I) = 0.0D0
   FZL(I) = 0.0D0
END DO

!$OMP DO SCHEDULE(STATIC,1)
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
                    VNBOND_ATOM = VNBOND_ATOM + VIJ
!                    WRITE(5100,*) I, J, VIJ*conv
                 ELSE IF(TYPE_LABEL(I) .eq. 2) THEN !BEAD
                    IF(TYPE_LABEL(J) .eq.2) THEN !I BEAD
                       VNBOND_BEAD = VNBOND_BEAD + VIJ
!                       WRITE(5200,*) I, J, VIJ*conv
                    ELSE IF(TYPE_LABEL(J) .eq. 3) THEN
                       VNBOND_MIX = VNBOND_MIX + VIJ
 !                      WRITE(5300,*) I, J, VIJ*conv
                    ELSE
                       WRITE(*,*) 'EXCEPTION!', I, J
                    END IF
                 ELSE IF(TYPE_LABEL(I) .eq. 3) THEN !VS
                    IF(TYPE_LABEL(J) .eq. 2) THEN
                       VNBOND_MIX = VNBOND_MIX + VIJ
  !                     WRITE(5300,*) I, J, VIJ*conv
                    ELSE
                       WRITE(*,*) 'EXCEPTION!', I, J
                    END IF
                 ELSE
                    WRITE(*,*) 'EXCEPTION', I, J
                 END IF

   !              WRITE(5000,*) I, J, VIJ*conv
#endif

                 VNBOND_TOTAL = VNBOND_TOTAL + VIJ

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
     END IF
     FXL(I) = FXL(I) + FXI
     FYL(I) = FYL(I) + FYI
     FZL(I) = FZL(I) + FZI
  END DO !End loop over all items
  !$OMP END DO

  !Collate forces from all threads
  DO I=1,NITEMS
     FX(I) = FX(I) + FXL(I)
     FY(I) = FY(I) + FYL(I)
     FZ(I) = FZ(I) + FZL(I)
  END DO

  !$OMP END PARALLEL

  RETURN

END SUBROUTINE NONBONDED_FORCE

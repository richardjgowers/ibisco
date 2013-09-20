#include "ibi-preprocess.h"

SUBROUTINE NONBONDED_FORCE(N,INDEX_LIST,MAXNAB,LIST,RCUT,RCUTSQ)

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N, INDEX_LIST(N), MAXNAB, LIST(MAXNAB,N)
  REAL(KIND=RKIND), INTENT(IN) :: RCUT, RCUTSQ
  INTEGER :: A, B, I, J, TI, TJ, TIJ
  INTEGER :: JNAB
  INTEGER :: NI
  REAL(KIND=RKIND) :: FXI,FYI,FZI, FXIJ, FYIJ, FZIJ
  REAL(KIND=RKIND) :: RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,RIJSQ,RIJ
  REAL(KIND=RKIND) :: ALPHA, FIJ, VIJ
  REAL(KIND=RKIND), DIMENSION(NITEMS) :: FXL, FYL, FZL !Local arrays for parallelisation
#ifdef DETAILED_ATOM
  INTEGER :: ME, GROUPSIZE, SPEC_ATOM = 1837, OUTFILE(0:7) !Will print out all information on a certain atom to sdout
#endif

#ifdef DETAILED_ATOM
  WRITE(*,*) FX(SPEC_ATOM), FY(SPEC_ATOM), FZ(SPEC_ATOM), 'Start of NONBONDED'

  OUTFILE(0) = 1001
  OUTFILE(1) = 1002
  OUTFILE(2) = 1003
  OUTFILE(3) = 1004
  OUTFILE(4) = 1005
  OUTFILE(5) = 1006
  OUTFILE(6) = 1007
  OUTFILE(7) = 1008
#endif

  !$OMP PARALLEL DEFAULT(NONE)&
  !$OMP& SHARED(N,B,NITEMS,INDEX_LIST,LIST,ITYPE,INBONDT,MAXNAB)&
  !$OMP& SHARED(RX,RY,RZ,BOXXINV,BOXYINV,BOXZINV,BOXX,BOXY,BOXZ,RCUTSQ)&
  !$OMP& SHARED(BINNB,RNBOND,NBOND_FORCE,NBOND_POT,TYPE_LABEL)&
#ifdef DETAILED_ATOM
  !$OMP& SHARED(SPEC_ATOM,GROUPSIZE,OUTFILE)&
  !$OMP& PRIVATE(ME)&
#endif
  !$OMP& PRIVATE(A,I,TI,J,TJ,TIJ,JNAB)&
  !$OMP& PRIVATE(FXI,FYI,FZI,RXI,RYI,RZI,RXIJ,RYIJ,RZIJ,RIJSQ,RIJ)&
  !$OMP& PRIVATE(NI,ALPHA,FIJ,VIJ,FXIJ,FYIJ,FZIJ)&
  !$OMP& PRIVATE(FXL,FYL,FZL)&
  !$OMP& REDUCTION(+:VNBOND_TOTAL,VNBOND_ATOM,VNBOND_BEAD,VNBOND_MIX)&
  !$OMP& REDUCTION(+:PT11,PT22,PT33,PT12,PT13,PT23) &
  !$OMP& REDUCTION(+:FX,FY,FZ)

#ifdef DETAILED_ATOM
ME = OMP_GET_THREAD_NUM()
GROUPSIZE = OMP_GET_NUM_THREADS()
#endif
  

FXL = 0.0D0
FYL = 0.0D0
FZL = 0.0D0

#ifdef DETAILED_ATOM
  WRITE(OUTFILE(ME),*) FXL(SPEC_ATOM), FYL(SPEC_ATOM), FZL(SPEC_ATOM), ME, GROUPSIZE, 'Start of NB, in parallel'
#endif

!$OMP DO SCHEDULE(STATIC,1)
  DO A=1,N
     I = INDEX_LIST(A) !I is the index of atom being considered
     RXI = RX(I)
     RYI = RY(I)
     RZI = RZ(I)

     FXI=0.0D0
     FYI=0.0D0
     FZI=0.0D0
     TI = ITYPE(I)

        DO JNAB = 1, MAXNAB

           J = LIST(JNAB,A)
           IF(J .EQ. 0) EXIT !Detect end of list


           TJ = ITYPE(J)
           TIJ = INBONDT(TI,TJ)

           IF(TIJ .NE. 0) THEN
              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J)
              RZIJ = RZI - RZ(J)

              RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
              RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
              RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ

              RIJSQ = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ

              IF(RIJSQ .LT. RCUTSQ) THEN
                 RIJ = SQRT(RIJSQ) !!double precision

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

                 FX(J) = FX(J) - FXIJ
                 FY(J) = FY(J) - FYIJ
                 FZ(J) = FZ(J) - FZIJ

#ifdef DETAILED_ATOM
    IF(J .eq. SPEC_ATOM) THEN
       WRITE(OUTFILE(ME),*) FX(J), -FXIJ, ME, 'Nonbonded small'
    END IF
#endif

                 PT11 = PT11 + FXIJ * RXIJ
                 PT22 = PT22 + FYIJ * RYIJ
                 PT33 = PT33 + FZIJ * RZIJ

                 PT12 = PT12 + FYIJ * RXIJ
                 PT13 = PT13 + FZIJ * RXIJ
                 PT23 = PT23 + FZIJ * RYIJ

              END IF !End if RIJSQ lt RCUTSQ	
           END IF !End if TIJ ne 0
        END DO !End loop over neighbours

     FX(I) = FX(I) + FXI
     FY(I) = FY(I) + FYI
     FZ(I) = FZ(I) + FZI

#ifdef DETAILED_ATOM
    IF(I .eq. SPEC_ATOM) THEN
       WRITE(OUTFILE(ME),*) FX(I), FXI, ME, 'Nonbonded large'
    END IF
#endif
 END DO !End loop over all items
  !$OMP END DO


#ifdef DETAILED_ATOM
       WRITE(OUTFILE(ME),*) FX(SPEC_ATOM), 'End of loop'
#endif

  !$OMP END PARALLEL

#ifdef DETAILED_ATOM
  WRITE(*,*) FX(SPEC_ATOM), FY(SPEC_ATOM), FZ(SPEC_ATOM), 'End of NONBONDED'
#endif

  RETURN

END SUBROUTINE NONBONDED_FORCE

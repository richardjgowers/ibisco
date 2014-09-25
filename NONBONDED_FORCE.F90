#ifndef DOXY_SKIP
#include "ibi-preprocess.h"
#endif
  !> @file
  !> @brief Calculates the nonbonded force in a passed particle list
  !> @details This subroutine is designed to work with any particle list that is passed to it.
  !!         In this way, atoms and beads can use the same subroutine.
  !> @note This subroutine is parallelised using OpenMP
  !> @author Richard J Gowers
  
  SUBROUTINE NONBONDED_FORCE(N,INDEX_LIST,MAXNAB,LIST,RCUT,RCUTSQ)

    USE VAR
    USE OMP_LIB

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N !< Size of the particle list you have passed to it
    INTEGER, INTENT(IN) :: INDEX_LIST(N) !< Index containing the address in master array of 
    !!each particle in this group
    INTEGER, INTENT(IN) :: MAXNAB !< The maximum number of neighbours that a particle could ever have
    INTEGER, INTENT(IN) :: LIST(MAXNAB,N) !< 2d array, first dimension goes over all particles, 
    !! second dimension contains neighbours for 
    !!             this particle
    REAL*4, INTENT(IN) :: RCUT !< Cutoff radius for nonbonded interactions
    REAL*4, INTENT(IN) :: RCUTSQ !< Cutoff radius squared
    INTEGER :: A, B, I, J, TI, TJ, TIJ
    INTEGER :: JNAB
    INTEGER :: NI
    REAL*4, DIMENSION(3) :: FXYZ_I
    REAL*4, DIMENSION(3) :: RXYZ_I, RXYZ_IJ
    REAL*4 :: RIJSQ, RIJ
    REAL*4 :: ALPHA, FIJ, VIJ

    !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC,1)&
    !$OMP& SHARED(N,B,NITEMS,INDEX_LIST,LIST,ITYPE,INBONDT,NNEBS)&
    !$OMP& SHARED(RXYZ,BOXINV,BOX,RCUTSQ)&
    !$OMP& SHARED(BINNB,RNBOND,NBOND_FORCE,NBOND_POT,TYPE_LABEL)&
    !$OMP& PRIVATE(A,I,TI,J,TJ,TIJ,JNAB)&
    !$OMP& PRIVATE(FXYZ_I,RXYZ_I,RXYZ_IJ,RIJSQ,RIJ)&
    !$OMP& PRIVATE(NI,ALPHA,FIJ,VIJ)&
    !$OMP& REDUCTION(+:V_NB)&
    !$OMP& REDUCTION(+:PT11,PT22,PT33,PT12,PT13,PT23) &
    !$OMP& REDUCTION(+:FXYZ)
    DO A=1,N
       I = INDEX_LIST(A) !I is the index of atom being considered
       RXYZ_I(:) = RXYZ(:,I)
       FXYZ_I(:) = 0.0

       TI = ITYPE(I)
       
       !DIR$ IVDEP
       DO JNAB = 1, NNEBS(I)
          J = LIST(JNAB,A)
          TJ = ITYPE(J)
          TIJ = INBONDT(TI,TJ)

          RXYZ_IJ(:) = RXYZ_I(:) - RXYZ(:,J)
          RXYZ_IJ(:) = RXYZ_IJ(:) - ANINT(RXYZ_IJ(:) * BOXINV(:)) * BOX(:)

          RIJSQ = SUM(RXYZ_IJ(:) * RXYZ_IJ(:)) ! R(1)*R(1) + R(2)*R(2) + ..etc

          IF(RIJSQ .LT. RCUTSQ) THEN
             RIJ = SQRT(RIJSQ) 

             NI = INT(RIJ / BINNB(TIJ))

             ALPHA = (RIJ - RNBOND(NI,TIJ)) / BINNB(TIJ)

             FIJ = NBOND_FORCE(NI,TIJ)*(1.0 - ALPHA) &
                  + NBOND_FORCE(NI+1,TIJ)*ALPHA

             VIJ = NBOND_POT(NI,TIJ)*(1.0 - ALPHA) &
                  + NBOND_POT(NI+1,TIJ)*ALPHA

             IF(TYPE_LABEL(I) .eq. 1) THEN !ATOM
                V_NB(1) = V_NB(1) + VIJ
             ELSE IF(TYPE_LABEL(I) .eq. 2) THEN !BEAD
                IF(TYPE_LABEL(J) .eq.2) THEN !BOTH BEADS
                   V_NB(2) = V_NB(2) + VIJ
                ELSE IF(TYPE_LABEL(J) .eq. 3) THEN!MIXED
                   V_NB(3) = V_NB(3) + VIJ
                END IF
             ELSE IF(TYPE_LABEL(I) .eq. 3) THEN !VS
                V_NB(3) = V_NB(3) + VIJ
             END IF

             FXYZ_I(:) = FXYZ_I(:) + FIJ * RXYZ_IJ(:)

             FXYZ(:,J) = FXYZ(:,J) - FIJ * RXYZ_IJ(:)

             PT11 = PT11 + FIJ * RXYZ_IJ(1) * RXYZ_IJ(1)
             PT22 = PT22 + FIJ * RXYZ_IJ(2) * RXYZ_IJ(2)
             PT33 = PT33 + FIJ * RXYZ_IJ(3) * RXYZ_IJ(3)

             PT12 = PT12 + FIJ * RXYZ_IJ(2) * RXYZ_IJ(1)
             PT13 = PT13 + FIJ * RXYZ_IJ(3) * RXYZ_IJ(1)
             PT23 = PT23 + FIJ * RXYZ_IJ(3) * RXYZ_IJ(2)

          END IF !End if RIJSQ lt RCUTSQ	
       END DO !End loop over neighbours

       FXYZ(:,I) = FXYZ(:,I) + FXYZ_I(:)
    END DO !End loop over all items
    !$OMP END PARALLEL DO

    RETURN

  END SUBROUTINE NONBONDED_FORCE

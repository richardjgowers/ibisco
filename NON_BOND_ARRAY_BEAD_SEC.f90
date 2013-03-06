
        SUBROUTINE NON_BOND_ARRAY_BEAD_SEC(BI,BK)

      USE VAR
        IMPLICIT NONE

      INTEGER            BI,BK,I, J,K,J1,K1,L1
!       *******************************************************************
                NONBOND = 1
                I = MIN(BI,BK)
                K = MAX(BI,BK)
!            IF (NONBEXC == 4) THEN
!!!!!!!!!!!!!!!!BOND PART: EXLUDE THE NONBONDED INTERACTION BETWEEN THE ATOMS BELONG TO THE SAME BOND!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO J = 1, NBONDS(I)
                   J1 = JBOND(I,J)
                   IF (K.EQ.J1) NONBOND = 0
            END DO

      RETURN
        END

!      ************************************************************************************

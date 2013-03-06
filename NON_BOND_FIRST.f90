        SUBROUTINE NON_BOND_FIRST(BI,BK)

      USE VAR
        IMPLICIT NONE

      INTEGER            BI,BK,I, J,K,J1,K1,L1
!       *******************************************************************
                NONBONDVIRT = 1
                I = MIN(BI,BK)
                K = MAX(BI,BK)
!            IF (NONBEXC == 4) THEN
!!!!!!!!!!!!!!!!BOND PART: EXLUDE THE NONBONDED INTERACTION BETWEEN THE ATOMS BELONG TO THE SAME BOND!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO J = 1, NBONDS(I)
                   J1 = JBOND(I,J)
                   IF (K.EQ.J1) NONBONDVIRT = 0
            END DO

!!!!!!!!!!!!!!!!!ANGLE PART: EXLUDE THE NONBONDED INTERACTION BETWEEN THE ATOMS BELONG TO THE SAME ANGLE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO J = 1, NIJK(I)
                   J1 = JANGLEIJK(I,J)
                   K1 = KANGLEIJK(I,J)
                   IF (K.EQ.J1.OR.K.EQ.K1) NONBONDVIRT = 0
            END DO

            DO J = 1, NOANGLEIJK(I)
                   J1 = NOJANGLEIJK(I,J)
                   K1 = NOKANGLEIJK(I,J)
                   IF (K.EQ.J1.OR.K.EQ.K1) NONBONDVIRT = 0
            END DO



      RETURN
        END

!      ************************************************************************************

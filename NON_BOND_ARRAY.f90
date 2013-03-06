
        SUBROUTINE NON_BOND_ARRAY(BI,BK)

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

!!!!!!!!!!!!!!!!!ANGLE PART: EXLUDE THE NONBONDED INTERACTION BETWEEN THE ATOMS BELONG TO THE SAME ANGLE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DO J = 1, NIJK(I)
                   J1 = JANGLEIJK(I,J)
                   K1 = KANGLEIJK(I,J)
                   IF (K.EQ.J1.OR.K.EQ.K1) NONBOND = 0
            END DO

            DO J = 1, NOANGLEIJK(I)
                   J1 = NOJANGLEIJK(I,J)
                   K1 = NOKANGLEIJK(I,J)
                   IF (K.EQ.J1.OR.K.EQ.K1) NONBOND = 0
            END DO

            IF (NONBEXC == 5) THEN
                   DO J = 1, NIJKL(I)
                      J1 = JTORIJKL(I,J)
                      K1 = KTORIJKL(I,J)
                      L1 = LTORIJKL(I,J)
                      IF(K.EQ.J1.OR.K.EQ.K1.OR.K.EQ.L1) NONBOND = 0
                   ENDDO
              DO J = 1, NOOPIJKL(I)
                      J1 = JOOPIJKL(I,J)
                      K1 = KOOPIJKL(I,J)
                      L1 = LOOPIJKL(I,J)
                      IF(K.EQ.J1.OR.K.EQ.K1.OR.K.EQ.L1) NONBOND = 0
                   ENDDO
                 DO J = 1, FNIJKL(I)
                      J1 = FJTORIJKL(I,J)
                      K1 = FKTORIJKL(I,J)
                      L1 = FLTORIJKL(I,J)
                      IF(K.EQ.J1.OR.K.EQ.K1.OR.K.EQ.L1) NONBOND = 0
                   ENDDO
                ENDIF

      RETURN
        END

!      ************************************************************************************

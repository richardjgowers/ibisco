      SUBROUTINE WRITETP ()

      USE VAR
      IMPLICIT NONE
      INTEGER :: I, N, NUMANG, NUMTOR, M , J, NUMMDNB
!       *******************************************************************

      NUMANG = 0
      NUMTOR = 0
      NUMMDNB = 0
      N = NATOMS

      DO I = 1, N 
            NUMANG = NUMANG + NIJK(I)
            NUMTOR = NUMTOR + NIJKL(I)
      END DO

      DO I = 1, N - 1
            DO J = 1 + I, N
!                  CALL NON_BOND_ARRAY(I,J) !FIX THIS
                  IF (NONBOND == 0) NUMMDNB = NUMMDNB + 1
            END DO
      END DO

9046 FORMAT (1('atoms:'))
9047 FORMAT (1('title:'))
9048 FORMAT (1('basta:'))
9049 FORMAT (1('angles:'))
9050 FORMAT (1('modified_nonbonded:'))
9051 FORMAT (1('torsions:'))

      WRITE(116,9047)
      WRITE(116,*)!TITLE
      WRITE(116,9046)
      WRITE(116,*)N

9045 FORMAT (I8,1X,A3,1X,E14.7,1X,3 (I4,1X))
      DO I = 1, N
            WRITE(116,9045)I,LABEL(ITYPE(I)),MASS0(ITYPE(I)), 0, 0, 0
      END DO

9052 FORMAT (6 (I8,1X))
      IF (NUMANG > 0) THEN
            WRITE(116,9049)
            WRITE(116,*)NUMANG
            M = 0

            DO I = 1, N
                  DO J = 1, NIJK(I)
                        M = M + 1
                        WRITE(116,9052)M, I, JANGLEIJK(I,J), KANGLEIJK(I,J), 0, 0
                  END DO
            END DO
      END IF 
	
9053 FORMAT (8 (I8,1X))

      IF (NUMTOR > 0) THEN 
            WRITE(116,9051)
            WRITE(116,*)NUMTOR
            M = 0
            DO I = 1, N

                  DO J = 1, NIJKL(I)
                        M = M + 1
                        WRITE(116,9053)M, I, JTORIJKL(I,J), KTORIJKL(I,J), LTORIJKL(I, J), 0, 0, 0
                  END DO
            END DO
      END IF

9054 FORMAT (7 (I8,1X))
      IF (NUMMDNB > 0) THEN
            WRITE(116,9050)
            WRITE(116,*)NUMMDNB
            M = 0
            
            DO I = 1, N - 1
                  DO J = I + 1, N

!                        CALL NON_BOND_ARRAY(I,J) !FIX THIS
                        IF (NONBOND == 0) THEN
                              M = M + 1
                              WRITE(116,9054)M, I, J, 0, 0, 0, 0
                        END IF
                  END DO
            END DO
      END IF 

      WRITE(116,9048)
      CLOSE(116)

      RETURN
      END
!	*********************************************************************************************

SUBROUTINE BUILD_CONNECTIVITY(N,INDEX_LIST,EXCLUSION)

USE VAR

IMPLICIT NONE

INTEGER, INTENT(IN) :: N, INDEX_LIST(N), EXCLUSION
INTEGER :: A, B, COUNTER
INTEGER :: I, J, K, L

CONNECTIONS = 0
CONNECTED_TO = 0

DO A=1,N !Cycle over 1 to NATOMS, NATOMS+1 to NITEMS is 0 because they are VS
   I = INDEX_LIST(A)
   COUNTER = 0

   !Loop over bonds
   DO B=1,NBONDS(I)
      J = JBOND(I,B)
      COUNTER = COUNTER +1
      IF(COUNTER .gt. MAXCONNECTIONS) THEN
         WRITE(1,*) 'ERROR, Number of connections greater than max connections'
         WRITE(*,*) 'ERROR, Number of connections greater than max connections'
         ISTOP = 1
         RETURN
      END IF
      CONNECTED_TO(I,COUNTER) = J
   END DO

   !Loop over angles
   DO B=1,NIJK(I)
      K = KANGLEIJK(I,B)
      COUNTER = COUNTER +1
      IF(COUNTER .gt. MAXCONNECTIONS) THEN
         WRITE(1,*) 'ERROR, Number of connections greater than max connections'
         WRITE(*,*) 'ERROR, Number of connections greater than max connections'
         ISTOP = 1
         RETURN
      END IF
      CONNECTED_TO(I,COUNTER) = K
   END DO

   DO B=1,NOANGLEIJK(I)
      K = NOKANGLEIJK(I,B)
      COUNTER = COUNTER +1
      IF(COUNTER .gt. MAXCONNECTIONS) THEN
         WRITE(1,*) 'ERROR, Number of connections greater than max connections'
         WRITE(*,*) 'ERROR, Number of connections greater than max connections'
         ISTOP = 1
         RETURN
      END IF
      CONNECTED_TO(I,COUNTER) = K
   END DO

   IF(EXCLUSION .eq. 5) THEN !If 1-4 ARE included
      !Loop over torsions
      DO B=1,NIJKL(I)
         L = LTORIJKL(I,B)
         COUNTER = COUNTER +1
         IF(COUNTER .gt. MAXCONNECTIONS) THEN
            WRITE(1,*) 'ERROR, Number of connections greater than max connections'
            WRITE(*,*) 'ERROR, Number of connections greater than max connections'
            ISTOP = 1
            RETURN
         END IF
         CONNECTED_TO(I,COUNTER) = L
      END DO

      DO B=1,FNIJKL(I)
         L = FLTORIJKL(I,B)
         COUNTER = COUNTER +1
         IF(COUNTER .gt. MAXCONNECTIONS) THEN
            WRITE(1,*) 'ERROR, Number of connections greater than max connections'
            WRITE(*,*) 'ERROR, Number of connections greater than max connections'
            ISTOP = 1
            RETURN
         END IF
         CONNECTED_TO(I,COUNTER) = L
      END DO

      !Loop over OOP
      DO B=1,NOOPIJKL(I)
         L = LOOPIJKL(I,B)
         COUNTER = COUNTER +1
         IF(COUNTER .gt. MAXCONNECTIONS) THEN
            WRITE(1,*) 'ERROR, Number of connections greater than max connections'
            WRITE(*,*) 'ERROR, Number of connections greater than max connections'
            ISTOP = 1
            RETURN
         END IF
         CONNECTED_TO(I,COUNTER) = L
      END DO
   END IF

   !Connections() beyond NATOMS is 0 as VS have no connections
   CONNECTIONS(I) = COUNTER
END DO

RETURN

END SUBROUTINE BUILD_CONNECTIVITY

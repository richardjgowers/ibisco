SUBROUTINE BUILD_CONNECTIVITY(N,INDEX_LIST,EXCLUSION)

USE VAR

IMPLICIT NONE

INTEGER, INTENT(IN) :: N, INDEX_LIST(N), EXCLUSION
INTEGER :: A, B, COUNTER
INTEGER :: I, J, K, L
LOGICAL, DIMENSION(NATOMS,NATOMS) :: IS_CONNECTED

IS_CONNECTED = .FALSE.

DO I=1,NATOMS !Cycle over 1 to NATOMS, NATOMS+1 to NITEMS is 0 because they are VS
   !Loop over bonds
   DO B=1,NBONDS(I)
      J = JBOND(I,B)
      IS_CONNECTED(I,J) = .true.
      IS_CONNECTED(J,I) = .true.
   END DO

   !Loop over angles
   DO B=1,NIJK(I)
      J = JANGLEIJK(I,B)
      K = KANGLEIJK(I,B)
      IS_CONNECTED(I,J) = .true.
      IS_CONNECTED(J,I) = .true.
      IS_CONNECTED(I,K) = .true.
      IS_CONNECTED(K,I) = .true.      
   END DO

   DO B=1,NOANGLEIJK(I)
      J = NOJANGLEIJK(I,B)
      K = NOKANGLEIJK(I,B)
      IS_CONNECTED(I,J) = .true.
      IS_CONNECTED(J,I) = .true.
      IS_CONNECTED(I,K) = .true.
      IS_CONNECTED(K,I) = .true.   
   END DO

   IF(EXCLUSION .eq. 5) THEN !If 1-4 ARE included
      !Loop over torsions
      DO B=1,NIJKL(I)
         J = JTORIJKL(I,B)
         K = KTORIJKL(I,B)
         L = LTORIJKL(I,B)
         IS_CONNECTED(I,J) = .true.
         IS_CONNECTED(J,I) = .true.
         IS_CONNECTED(I,K) = .true.
         IS_CONNECTED(K,I) = .true.
         IS_CONNECTED(I,L) = .true.
         IS_CONNECTED(L,I) = .true.
      END DO

      DO B=1,FNIJKL(I)
         J = FJTORIJKL(I,B)
         K = FKTORIJKL(I,B)
         L = FLTORIJKL(I,B)
         IS_CONNECTED(I,J) = .true.
         IS_CONNECTED(J,I) = .true.
         IS_CONNECTED(I,K) = .true.
         IS_CONNECTED(K,I) = .true.
         IS_CONNECTED(I,L) = .true.
         IS_CONNECTED(L,I) = .true.
      END DO

      !Loop over OOP
      DO B=1,NOOPIJKL(I)
         J = JOOPIJKL(I,B)
         K = KOOPIJKL(I,B)
         L = LOOPIJKL(I,B)
         IS_CONNECTED(I,J) = .true.
         IS_CONNECTED(J,I) = .true.
         IS_CONNECTED(I,K) = .true.
         IS_CONNECTED(K,I) = .true.
         IS_CONNECTED(I,L) = .true.
         IS_CONNECTED(L,I) = .true.
      END DO
   END IF !IF EXCLUSION .eq. 5 
END DO

DO A=1,N
   COUNTER = 0
   I = INDEX_LIST(A)

   IF(I .gt. NATOMS) THEN !This is now a VS, as VS are appended after beads, all beads are done
      RETURN
   END IF

   DO J=1,NATOMS
      IF(IS_CONNECTED(I,J)) THEN
         COUNTER = COUNTER +1
         CONNECTED_TO(I,COUNTER) = J
      END IF
   END DO
   CONNECTIONS(I) = COUNTER
END DO



RETURN

END SUBROUTINE BUILD_CONNECTIVITY

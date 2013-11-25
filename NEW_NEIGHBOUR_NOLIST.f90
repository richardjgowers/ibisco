SUBROUTINE NEW_NEIGHBOUR_NOLIST(INDEX_LIST,POINT,CELL,LCLIST,N,RLIST,LIST,MAXNAB &
     , MAP, SIZEMAP, HEAD, MAXNUMCELL)

  USE VAR

  INTEGER, INTENT(IN) :: N, MAXNUMCELL, MAXNAB
  INTEGER, INTENT(IN) :: INDEX_LIST(N), CELL(N), LCLIST(N), SIZEMAP, MAP(SIZEMAP), HEAD(MAXNUMCELL)
  INTEGER, INTENT(INOUT) :: POINT(N+1), LIST(MAXNAB)
  INTEGER :: A, B, C, D, I, J, NLIST
  REAL(KIND=RKIND), INTENT(IN) :: RLIST
  REAL(KIND=RKIND) :: RXI, RYI, RZI
  REAL(KIND=RKIND) :: RIJSQ
  REAL(KIND=RKIND) :: RXIJ, RYIJ, RZIJ
  REAL(KIND=RKIND) :: RLISTSQ

  RLISTSQ = RLIST * RLIST

  NLIST = 0
  DO A=1,N-1
     I = INDEX_LIST(A)
     POINT(A) = NLIST+1
     RXI = RX(I)
     RYI = RY(I)
     RZI = RZ(I)

     DO B=A+1,N
        J = INDEX_LIST(B)
        !Check if connected
        NONBOND=1
        IF(CONNECTIONS(I) .gt. 0) THEN !If I has connections:
           DO C=1,CONNECTIONS(I)
              D = CONNECTED_TO(I,C)
              IF(J .eq. D) THEN !If J is C then they are connected
                 NONBOND=0 
                 EXIT
              END IF
           END DO
        END IF
        !Check if VS-VS
        IF(IBRDESCR .eq. 0) THEN
           IF(I .gt. NATOMS .and. J .gt. NATOMS) THEN
              NONBOND = 0
           END IF
        END IF

        IF(NONBOND .eq. 1) THEN
           RXIJ  = RXI - RX(J)
           RYIJ  = RYI - RY(J)
           RZIJ  = RZI - RZ(J)

           RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
           RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
           RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ
           RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

           IF (RIJSQ.LT.RLISTSQ) THEN
              NLIST = NLIST +1
              LIST(NLIST) = J              
              IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
           END IF
        END IF
        J = LCLIST(A)
     END DO

  END DO
  POINT(N) = NLIST +1
  POINT(N+1) = NLIST +1

RETURN 
END SUBROUTINE NEW_NEIGHBOUR_NOLIST

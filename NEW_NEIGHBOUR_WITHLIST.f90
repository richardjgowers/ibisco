SUBROUTINE NEW_NEIGHBOUR_WITHLIST(INDEX_LIST,POINT,CELL,LCLIST,N,RLIST,LIST,MAXNAB &
     , MAP, SIZEMAP, HEAD, MAXNUMCELL)

  USE VAR

  IMPLICIT	NONE

  INTEGER, INTENT(IN) :: N, MAXNUMCELL, MAXNAB
  INTEGER, INTENT(IN) :: INDEX_LIST(N), CELL(N), LCLIST(N), SIZEMAP, MAP(SIZEMAP), HEAD(MAXNUMCELL)
  INTEGER, INTENT(INOUT) :: POINT(N+1), LIST(MAXNAB)
  INTEGER :: A, B, C, D, I, J, NLIST, JCELL, JCELL0, NABOR
  REAL(KIND=RKIND), INTENT(IN) :: RLIST
  REAL(KIND=RKIND) :: RXI, RYI, RZI
  REAL(KIND=RKIND) :: RIJSQ
  REAL(KIND=RKIND) :: RXIJ, RYIJ, RZIJ
  REAL(KIND=RKIND) :: RLISTSQ


  RLISTSQ = RLIST * RLIST

  !***************************BUILD UP THE NEIGHBOUR LIST*********************
  !***********LOOP IN THE CELL WHERE ATOM I IS **************************
  NLIST = 0            !!!!!!!!!!EXCLUDE THE BONDED PAIRS FROM INTERACTION LIST**************
  DO A = 1,N
     I = INDEX_LIST(A)

     POINT(A) = NLIST + 1
     JCELL = CELL(A)
     RXI = RX(I)
     RYI = RY(I)
     RZI = RZ(I)
  
     B = LCLIST(A)
     DO
        IF(B .eq. 0) EXIT !If reached end of list
        J = INDEX_LIST(B) !Fetch real index of B

        !           if(abs(ip-jp) .le. contactA)then
        NONBOND=1
        IF(CONNECTIONS(I) .gt. 0) THEN !If I has connections:
           DO C=1,CONNECTIONS(I) !Connections is the number of connected atoms to I
              D = CONNECTED_TO(I,C) !Connected_to has the index of all connected atoms
              IF(J .eq. D) THEN !If J is D then they are connected
                 NONBOND=0 
                 EXIT
              END IF
           END DO
        END IF
        IF(IBRDESCR .eq. 0) THEN
           IF(I .gt. NATOMS .and. J .gt. NATOMS) THEN !Exclude VS-VS interactions
              NONBOND = 0
           END IF
        END IF

        IF (NONBOND.EQ.1) THEN
           RXIJ = RXI - RX(J)
           RYIJ = RYI - RY(J)
           RZIJ = RZI - RZ(J)

           RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
           RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
           RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
           RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2

           IF (RIJSQ.LT.RLISTSQ) THEN
              NLIST = NLIST + 1
              LIST(NLIST) = J
              IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
           ENDIF
        ENDIF
        B = LCLIST(B) !Error is here
     END DO

        !***********LOOP IN THE NEIGHBOURING CELLS**********************************
     JCELL0 = 13*(JCELL - 1)
     DO NABOR = 1,13
        JCELL = MAP(JCELL0+NABOR)
        B = HEAD(JCELL)

        DO
           IF(B .eq. 0) EXIT
           J = INDEX_LIST(B)
           !              if(abs(ip-jp) .le. contactA)then
           !              else
           !                 NONBOND=1
           !              end if
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
           IF(IBRDESCR .eq. 0) THEN
              IF(I .gt. NATOMS .and. J .gt. NATOMS) THEN
                 NONBOND = 0
              END IF
           END IF

           IF (NONBOND.EQ.1) THEN
              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J) 
              RZIJ = RZI - RZ(J)

              RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
              RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
              RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
              RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2

              IF (RIJSQ.LT.RLISTSQ) THEN
                 NLIST = NLIST + 1
                 LIST(NLIST) = J
                 IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
              ENDIF

           ENDIF
           B = LCLIST(B)
        END DO
     END DO

  END DO
  POINT(N+1) = NLIST + 1

  RETURN
END SUBROUTINE NEW_NEIGHBOUR_WITHLIST


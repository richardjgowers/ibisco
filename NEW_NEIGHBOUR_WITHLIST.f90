!> @file
!> @brief Neighbourlist creation with linked cell lists
!> @details Called from UPDATE_NEIGHBOURLIST.f90

SUBROUTINE NEW_NEIGHBOUR_WITHLIST(INDEX_LIST,CELL,LCLIST,N,RLIST,LIST,MAXNAB &
     , MAP, SIZEMAP, HEAD, MAXNUMCELL)

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N !! @param n Number of particles to process
  INTEGER, INTENT(IN) :: MAXNUMCELL, MAXNAB
  INTEGER, INTENT(IN) :: INDEX_LIST(N) !! Array which translates onto the master coordinate list
  INTEGER, INTENT(IN) :: CELL(N), LCLIST(N), SIZEMAP, MAP(SIZEMAP), HEAD(MAXNUMCELL)
  INTEGER, INTENT(INOUT) :: LIST(MAXNAB,N)
  INTEGER :: A, B, C, D, I, J, NLIST, JCELL, JCELL0, NABOR
  INTEGER :: NONBOND
  REAL(KIND=RKIND), INTENT(IN) :: RLIST
  REAL(KIND=RKIND) :: RXI, RYI, RZI
  REAL(KIND=RKIND) :: RIJSQ
  REAL(KIND=RKIND) :: RXIJ, RYIJ, RZIJ
  REAL(KIND=RKIND) :: RLISTSQ


  RLISTSQ = RLIST * RLIST
  LIST = 0 !Should set all unoccupied list slots to 0, can then use this to detect end of list
  !***************************BUILD UP THE NEIGHBOUR LIST*********************
  !***********LOOP IN THE CELL WHERE ATOM I IS **************************
  !$OMP PARALLEL DO DEFAULT(NONE), SCHEDULE(STATIC,1)&
  !$OMP& SHARED(N,INDEX_LIST,CELL,RX,RY,RZ,LCLIST,CONNECTIONS,CONNECTED_TO,IBRDESCR,NATOMS)&
  !$OMP& SHARED(BOXX,BOXY,BOXZ,BOXXINV,BOXYINV,BOXZINV,RLISTSQ,LIST,MAXNAB,MAP,HEAD)&
  !$OMP& PRIVATE(A,I,NLIST,JCELL,RXI,RYI,RZI,B,J,NONBOND,C,D,RXIJ,RYIJ,RZIJ,RIJSQ,JCELL0,NABOR)
  DO A = 1,N
     I = INDEX_LIST(A)
     NLIST = 0 !Counter for number of neighbours this atom has

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
           RIJSQ = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ

           IF (RIJSQ.LT.RLISTSQ) THEN
              NLIST = NLIST + 1
              LIST(NLIST,A) = J
              IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
           ENDIF
        ENDIF
        B = LCLIST(B)
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
              RIJSQ = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ

              IF (RIJSQ.LT.RLISTSQ) THEN
                 NLIST = NLIST + 1
                 LIST(NLIST,A) = J
                 IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
              ENDIF

           ENDIF
           B = LCLIST(B)
        END DO
     END DO

  END DO
  !$OMP END PARALLEL DO

  RETURN
END SUBROUTINE NEW_NEIGHBOUR_WITHLIST


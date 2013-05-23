SUBROUTINE NEW_NEIGHBOUR_WITHLIST (INDEX_LIST,POINT,CELL,LCLIST,N,RLIST,LIST,MAXNAB &
     , MAP, SIZEMAP, HEAD, MAXNUMCELL)

  USE VAR

  IMPLICIT	NONE

  INTEGER, INTENT(IN) :: N, MAXNUMCELL, MAXNAB
  INTEGER, INTENT(IN) :: INDEX_LIST(N), CELL(NITEMS), LCLIST(NITEMS), SIZEMAP, MAP(SIZEMAP), HEAD(MAXNUMCELL)
  INTEGER, INTENT(INOUT) :: POINT(N+1), LIST(MAXNAB)
  INTEGER :: A, B, C, I, J, NLIST, JCELL, JCELL0, NABOR
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
     JCELL = CELL(I)
     RXI = RX(I)
     RYI = RY(I)
     RZI = RZ(I)

     J = LCLIST(I)
     DO
        IF(J .eq. 0) EXIT
        !           if(abs(ip-jp) .le. contactA)then
        NONBOND=1
        IF(CONNECTIONS(I) .gt. 0) THEN !If I has connections:
           DO B=1,CONNECTIONS(I)
              C = CONNECTED_TO(I,B)
              IF(J .eq. C) THEN !If J is C then they are connected
                 NONBOND=0 
                 EXIT
              END IF
           END DO
        END IF
        !           else
        !              NONBOND=1
        !           end if
        !call NON_BOND_ARRAY(IP,JP)
        IF (NONBOND.EQ.1) THEN
           RXIJ = RXI - RX(J)
           RYIJ = RYI - RY(J)
           RZIJ = RZI - RZ(J)

           RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
           RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
           RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
           RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2

!if(type_label(i) .eq. 1 .and.type_label(j) .eq. 1 .and. i .eq. 2 )write(*,*)i,j,jcell, RIJSQ,RLISTSQ,RLIST

           IF (RIJSQ.LT.RLISTSQ) THEN
 !             if(type_label(i) .eq. 1 .and.type_label(j) .eq. 1 .and. i .eq. 2 )write(*,*)'giao'
              NLIST = NLIST + 1
              LIST(NLIST) = J
              IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
           ENDIF
        ENDIF
        J = LCLIST(J)
     END DO

        !***********LOOP IN THE NEIGHBOURING CELLS**********************************
     JCELL0 = 13*(JCELL - 1)
     DO NABOR = 1,13
        JCELL = MAP(JCELL0+NABOR)
        J = HEAD(JCELL)
        DO
           IF(J .eq. 0) EXIT
           !              if(abs(ip-jp) .le. contactA)then
           !              else
           !                 NONBOND=1
           !              end if

           NONBOND=1
           IF(CONNECTIONS(I) .gt. 0) THEN !If I has connections:
              DO B=1,CONNECTIONS(I)
                 C = CONNECTED_TO(I,B)
                 IF(J .eq. C) THEN !If J is C then they are connected
                    NONBOND=0 
                    EXIT
                 END IF
              END DO
           END IF

           IF (NONBOND.EQ.1) THEN
              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J) 
              RZIJ = RZI - RZ(J)

              RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
              RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
              RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
              RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2
!if(type_label(i) .eq. 1 .and.type_label(j) .eq. 1 .and. i .eq. 2 )write(*,*)i,j,jcell, RIJSQ,RLISTSQ,RLIST


              IF (RIJSQ.LT.RLISTSQ) THEN
!if(type_label(i) .eq. 1 .and.type_label(j) .eq. 1 .and. i .eq. 2 )write(*,*)'giao'
                 NLIST = NLIST + 1
                 LIST(NLIST) = J
                 IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
              ENDIF
           ENDIF
           J = LCLIST(J)
        END DO
     END DO

  END DO
  POINT(N+1) = NLIST + 1

  WRITE(*,*) POINT(1), POINT(2), POINT(3)

  RETURN
END SUBROUTINE NEW_NEIGHBOUR_WITHLIST


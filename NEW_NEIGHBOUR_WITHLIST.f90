SUBROUTINE NEW_NEIGHBOR_WITHLIST (INDEX_LIST,SIZEINDEX_LIST,POINT,CELL,LCLIST,N,RLIST,LIST,MAXNAB &
     , MAP, SIZEMAP)

  USE VAR

  IMPLICIT	NONE

  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN) :: INDEX_LIST(N), CELL(N), LCLIST(N), MAXNAB, SIZEMAP, MAP(SIZEMAP)
  REAL(KIND=RKIND), INTENT(IN) :: RLIST
  INTEGER, INTENT(INOUT) :: POINT(N+1), LIST(MAXNAB)

  INTEGER :: A, I, J, NLIST, JCELL0, NABOR
  REAL(KIND=RKIND) :: RLISTSQ

  INTEGER :: NLIST, JCELL
  REAL*8 :: RXI, RYI, RZI
  REAL*8 :: RIJSQ
  REAL*8 :: RXIJ, RYIJ, RZIJ

  RLISTSQ = RLIST * RLIST

  !***************************BUILD UP THE NEIGHBOUR LIST*********************
  !***********LOOP IN THE CELL WHERE ATOM I IS **************************
  NLIST = 0            !!!!!!!!!!EXCLUDE THE BONDED PAIRS FROM INTERACTION LIST**************
  DO A = 1,N
     I = INDEX_LIST(A)
     POINT(I) = NLIST + 1
     JCELL = CELL(I)
     RXI = RX(I)
     RYI = RY(I)
     RZI = RZ(I)

     J = LCLIST(I)
     DO
        IF(J .eq. 0) EXIT
        !           if(abs(ip-jp) .le. contactA)then
        CALL NON_BOND_ARRAY(I,J)
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

           IF (RIJSQ.LT.RLISTSQ) THEN
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
           CALL NON_BOND_ARRAY(I,J)
           !              else
           !                 NONBOND=1
           !              end if

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
           J = LCLIST(J)
        END DO
     END DO

  END DO
  POINT(N+1) = NLIST + 1

  RETURN
END SUBROUTINE NEW_NEIGHBOR_WITHLIST


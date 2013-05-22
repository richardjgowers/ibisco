SUBROUTINE LINKS(HEAD,MAXNUMCELL,INDEX_LIST,N,CELL,NCELLX, NCELLY, NCELLZ, LCLIST)
  USE VAR
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: MAXNUMCELL, N, INDEX_LIST(N)
  INTEGER, INTENT(IN) :: NCELLX, NCELLY, NCELLZ
  INTEGER, INTENT(INOUT) :: HEAD(MAXNUMCELL), CELL(N), LCLIST(N)
  INTEGER I,J,JCELL

  HEAD = 0
  DO I = 1,N
     J = INDEX_LIST(I)
     JCELL = 1 + INT((RX(J)*BOXXINV + 0.5D0)*NCELLX)*NCELLZ*NCELLY                     &
          + INT((RY(J)*BOXYINV + 0.5D0)*NCELLY)*NCELLZ              &
          + INT((RZ(J)*BOXZINV + 0.5D0)*NCELLZ)
     LCLIST(J) = HEAD(JCELL)
     HEAD(JCELL) = J
     CELL(J) = JCELL

  END DO
  RETURN
END SUBROUTINE LINKS

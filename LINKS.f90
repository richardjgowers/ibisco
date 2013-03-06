         SUBROUTINE LINKS()
         USE VAR
         IMPLICIT NONE
         INTEGER I,JCELL
       
         HEAD = 0
         DO I = 1,NATOMS
            JCELL = 1 + INT((RX(I)*BOXXINV + 0.5D0)*NCELLX)*NCELLZ*NCELLY                     &
                      + INT((RY(I)*BOXYINV + 0.5D0)*NCELLY)*NCELLZ              &
                      + INT((RZ(I)*BOXZINV + 0.5D0)*NCELLZ)
            LCLIST(I) = HEAD(JCELL)
            HEAD(JCELL) = I
            CELL(I) = JCELL
!if(i .eq. 3727)then
!write(666,*)i,'CELL',jcell
!write(666,*)'LCLIST(I)',LCLIST(I),'HEAD(JCELL)',HEAD(JCELL)
!end if
         END DO
         RETURN
         END

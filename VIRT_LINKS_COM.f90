      SUBROUTINE VIRT_LINKS_COM()
      USE VAR
      IMPLICIT NONE
      INTEGER I,JCELL

      VHEAD = 0
      DO I = 1,NVIRTA

            JCELL = 1 + INT((VIRTRX(I)*BOXXINV + 0.5D0)*NCELLX)*NCELLZ*NCELLY                     &
                      + INT((VIRTRY(I)*BOXYINV + 0.5D0)*NCELLY)*NCELLZ              &
                      + INT((VIRTRZ(I)*BOXZINV + 0.5D0)*NCELLZ)
            VLCLIST(I) = VHEAD(JCELL)
!          VHEAD(JCELL) = INDEX_VSITE(I)
            VHEAD(JCELL) = I
            VCELL(I) = JCELL  

      END DO

      RETURN

      END


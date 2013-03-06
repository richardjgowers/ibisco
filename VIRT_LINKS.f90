!
!	Created by Nicodemo
!	Dec 2010
!
!	JCELL = number of the cell in which the atom I is contained
!	
!	When prog find the cell in which atom I is contained it checks the vector VHEAD, if atom I is the first 
!	atom found in that cell VHEAD(JCELL) == 0 otherwise it contains the index of the last atom found in that
!	cell. After this, prog replaces the value of vector VHEAD at the point JCELL with the new atom found in 
!	that position. At the end the position of the atom I in stored in vector VCELL(I).   
!	In this way is possible to find all the position of the atom starting from the last atom found in a specific
!	cell.
!
!	Ex. :
!
!	i = 1
!	jcell = 4
!	vlclist(1) = vhead(4) = 0
!	vhead(4) = 1
!	vcell(1) = 4

!	i = 2
!	jcell = 4
!	vlclist(2) = vhead(4) = 1
!	vhead(4) = 2
!	vcell(2) = 4

!	i = 3
!	jcell = 4
!	vlclist(3) = vhead(4) = 2
!	vhead(4) = 3
!	vcell(3) = 4

!	i = 4
!	jcell = 5
!	vlclist(4) = vhead(5) = 0
!	vhead(5) = 4
!	vcell(4) = 5

!	i = 5
!	jcell = 4
!	vlclist(5) = vhead(4) = 3
!	vhead(4) = 5
!	vcell(3) = 4
!
!	ecc. ecc.
!


      SUBROUTINE VIRT_LINKS()
      USE VAR
      IMPLICIT NONE
      INTEGER I,JCELL

      VHEAD = 0

      DO I = 1,NVIRTA
            JCELL = 1 + INT((RX(INDEX_VSITE(I))*BOXXINV + 0.5D0)*NCELLX)*NCELLZ*NCELLY                     &
                      + INT((RY(INDEX_VSITE(I))*BOXYINV + 0.5D0)*NCELLY)*NCELLZ              &
                      + INT((RZ(INDEX_VSITE(I))*BOXZINV + 0.5D0)*NCELLZ)
            VLCLIST(I) = VHEAD(JCELL)
!          VHEAD(JCELL) = INDEX_VSITE(I)
            VHEAD(JCELL) = I
            VCELL(I) = JCELL  
      END DO

      RETURN

      END






















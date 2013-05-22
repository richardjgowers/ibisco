!	SETS UP A LIST OF THE TWENTY SIX NEIGHBOURING CELLS OF EACH OF
!	             THE SMALL CELLS IN THE CENTRAL BOX.

! There are 6 cells for each face of the cell, 12 cells for each edge and 8 cells for each vertex. 26 neighbour cells for each cell

! In the vector MAP is stored the position of the 13 neighbours cells for each cell, the counter (IMAP) advances in 13 of 13.
! e.g. For the cell 1
! ICELL(1,1,1) = 1
! IMAP = (ICELL ( IX, IY, IZ ) - 1 ) * 13 = 0
! so, the firt 13 positions in the vector MAP will store the index of the neighbour of cell 1
! For the cell 2
! ICELL(1,1,2) = 2
! IMAP = (ICELL ( IX, IY, IZ ) - 1 ) * 13 = 13
! so, the positions from 14 to 26 in the vector MAP will store the index of the neighbour of cell 2




SUBROUTINE MAPS (MAP,MAPSIZE,NCELLX,NCELLY,NCELLZ)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NCELLX, NCELLY, NCELLZ,MAPSIZE
  INTEGER, INTENT(INOUT) :: MAP(MAPSIZE)
  INTEGER :: IX, IY, IZ, IMAP, ICELL
  !       *******************************************************************

  !    ** FIND THE NEAREST NEIGHBOURS OF EACH CELL **

  DO IX = 1, NCELLX !50

     DO IY = 1, NCELLY !40

        DO IZ = 1, NCELLZ !30

           IMAP = ( ICELL ( IX, IY, IZ ,NCELLX, NCELLY, NCELLZ  ) - 1 ) * 13

           MAP( IMAP + 1  ) = ICELL( IX + 1, IY    , IZ     ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 2  ) = ICELL( IX + 1, IY + 1, IZ     ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 3  ) = ICELL( IX    , IY + 1, IZ     ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 4  ) = ICELL( IX - 1, IY + 1, IZ     ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 5  ) = ICELL( IX + 1, IY    , IZ - 1 ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 6  ) = ICELL( IX + 1, IY + 1, IZ - 1 ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 7  ) = ICELL( IX    , IY + 1, IZ - 1 ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 8  ) = ICELL( IX - 1, IY + 1, IZ - 1 ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 9  ) = ICELL( IX + 1, IY    , IZ + 1 ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 10 ) = ICELL( IX + 1, IY + 1, IZ + 1 ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 11 ) = ICELL( IX    , IY + 1, IZ + 1 ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 12 ) = ICELL( IX - 1, IY + 1, IZ + 1 ,NCELLX, NCELLY, NCELLZ  )
           MAP( IMAP + 13 ) = ICELL( IX    , IY    , IZ + 1 ,NCELLX, NCELLY, NCELLZ  )

        END DO! 30            CONTINUE

     END DO! 40         CONTINUE

  END DO!50      CONTINUE

        RETURN
        END

!       *******************************************************************
!    ** STATEMENT FUNCTION TO GIVE CELL INDEX **
      FUNCTION ICELL ( IX, IY, IZ,NCELLX, NCELLY, NCELLZ )

!      USE VAR
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: NCELLX, NCELLY, NCELLZ
      INTEGER :: IX, IY, IZ, ICELL
!       *******************************************************************

      ICELL = 1 + MOD ( IZ - 1 + NCELLZ, NCELLZ )     &
                  + MOD ( IY - 1 + NCELLY, NCELLY ) * NCELLZ &
                  + MOD ( IX - 1 + NCELLX, NCELLX ) * NCELLZ * NCELLY

      END


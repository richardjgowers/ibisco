SUBROUTINE UPDATE_NEIGHBOURLIST()

USE VAR

IMPLICIT NONE

INTEGER :: A, I

IF (ENSEMBLE == 2) THEN !If NPT
   !Update the size of cells and corresponding map of cells.
   NCELLX_ATOM = BOXX / RLIST_ATOM
   NCELLY_ATOM = BOXY / RLIST_ATOM
   NCELLZ_ATOM = BOXZ / RLIST_ATOM

   NUMCELL_ATOM = NCELLX_ATOM * NCELLY_ATOM * NCELLZ_ATOM 

   CALL MAPS (MAP_ATOM,MAPSIZE_ATOM &
        , NCELLX_ATOM, NCELLY_ATOM, NCELLZ_ATOM) 


   IF(IBRDESCR .eq. 0) THEN

      NCELLX_BEAD = BOXX / RLIST_BEAD
      NCELLY_BEAD = BOXY / RLIST_BEAD
      NCELLZ_BEAD = BOXZ / RLIST_BEAD

      NUMCELL_BEAD = NCELLX_BEAD * NCELLY_BEAD * NCELLZ_BEAD 

      CALL MAPS (MAP_BEAD,MAPSIZE_BEAD &
           , NCELLX_BEAD, NCELLY_BEAD, NCELLZ_BEAD) 
   END IF
END IF

write(*,*) 'LINKS'

CALL LINKS (HEAD_ATOM,MAXNUMCELL_ATOM,ATOM,NUMATOMS,CELL_ATOM &
     , NCELLX_ATOM, NCELLY_ATOM, NCELLZ_ATOM, LCLIST_ATOM)

IF(IBRDESCR .eq. 0) THEN
   CALL LINKS (HEAD_BEAD,MAXNUMCELL_BEAD,BEAD,NCOARSE,CELL_BEAD &
        , NCELLX_BEAD, NCELLY_BEAD, NCELLZ_BEAD, LCLIST_BEAD)
END IF

write(*,*) 'Neighbour list'

CALL NEW_NEIGHBOUR_WITHLIST(ATOM,POINT_ATOM,CELL_ATOM,LCLIST_ATOM,NUMATOMS,RLIST_ATOM &
     , LIST_ATOM, MAXNAB_ATOM, MAP_ATOM, MAXMAPSIZE_ATOM, HEAD_ATOM, MAXNUMCELL_ATOM)

IF(IBRDESCR .eq. 0) THEN
   CALL NEW_NEIGHBOUR_WITHLIST(BEAD,POINT_BEAD,CELL_BEAD,LCLIST_BEAD,NCOARSE,RLIST_BEAD &
        , LIST_BEAD, MAXNAB_BEAD, MAP_BEAD, MAXMAPSIZE_BEAD, HEAD_BEAD, MAXNUMCELL_BEAD)
END IF

RETURN
END SUBROUTINE UPDATE_NEIGHBOURLIST

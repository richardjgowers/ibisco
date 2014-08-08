!> @file
!> @brief Generates data structures for hybrid simulations
!!
!> @details Hybrid simulations work on the basis of having 2 pointer arrays which point to the master
!! coordinate list.  These two pointer lists are called ATOM and BEAD. These are integer arrays 
!! which are the length of the number of atoms or beads.  ATOM(1) returns the index of the first atom
!! in the master arrays, eg RXYZ(1, ATOM(1)) is the x coordinate of the first atom.
!!
!! In this way, all atoms or beads can easily be cycled through.
!! 
!! called from main.f90

SUBROUTINE MAKE_LISTS()

USE VAR

IMPLICIT NONE

INTEGER :: I, NUMATOM, NUMBEADVS

NUMATOMS =0
NUMBEADS =0

DO I=1,NATOMS
   IF(TYPE_LABEL(I) .eq. 1) THEN !1 is atom
      NUMATOMS = NUMATOMS +1
   ELSE IF(TYPE_LABEL(I) .eq. 2) THEN!type_label = 2 is bead
      NUMBEADS = NUMBEADS +1
   ELSE
      WRITE(1,*) 'ATOM ',I,' HAS INVALID TYPE LABEL, CHECK interaction FILE'
      WRITE(*,*) 'ATOM ',I,' HAS INVALID TYPE LABEL, CHECK interaction FILE'
      ISTOP = 1
      RETURN
   END IF
END DO

NCOARSE = NUMBEADS + NVIRTA
ALLOCATE( ATOM(NUMATOMS) &
     , BEAD(NCOARSE) )

ATOM = 0
BEAD = 0
NUMATOM =0 
NUMBEADVS =0

!Make atom list
DO I=1,NATOMS
   IF(type_label(I) .eq. 1) then !If atom
      NUMATOM = NUMATOM +1      
      ATOM(NUMATOM) = I
   ELSE
      NUMBEADVS = NUMBEADVS +1
      BEAD(NUMBEADVS) = I
   END IF
END DO

DO I=1,NVIRTA
   NUMBEADVS = NUMBEADVS +1
   BEAD(NUMBEADVS) = NATOMS+I
END DO

RETURN

END SUBROUTINE MAKE_LISTS

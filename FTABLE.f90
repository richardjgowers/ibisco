!> @file
!> @brief Calculates the forcefield tables
!!
!> @details Forcefield input is tabulated potential against distance.  To turn this into
!! force information, the derivative at each point is taken.  
!!
!! Called from main.f90

SUBROUTINE FTABLE ()

  USE VAR
  IMPLICIT NONE

  INTEGER :: I, J
  INTEGER :: N_ENTRY, E
  REAL*8 :: RAMP_SIZE, SCALE
  !      *******************************************************************

  !	TABLE FORCE FOR BONDED ATOMS AND NON-BONDED INRTRACTIONS WILL BE DERIVATIVE
  !	OF POTENTIALS DIVIDED BY DISTANCE

  !	IF YOU WANT TO USE TABLE FOR BOND AND BEND POTENTIALS,
  !	MAKE TABLE FOR FORCE FOR THESE INTERACTIONS

  IF (INTERACT == 1) THEN

     !	MAKE TABLE FORCE FOR BOND
     DO I = 1, NBTYPE
        if(.not. typeBond(i))then
           DO J = 1, NDATB(I) - 1 
              BOND_FORCE(I,J) = -(BOND_POT(I,J+1) - BOND_POT(I,J-1))/BINB(I)/RBOND(I,J)/2.0D0
           END DO
           BOND_FORCE(I, 0)     = BOND_FORCE(I, 1)
           BOND_FORCE(I, NDATB) = BOND_FORCE(I, NDATB(I)-1)
        end if
     END DO

     !	MAKE TABLE FORCE FOR BEND
     DO I = 1, NATYPE

        DO J = 1, NDATAN(I) - 1
           BEND_FORCE(I,J) = -(BEND_POT(I,J+1) - BEND_POT(I,J-1))/BINA(I)/D2R/2.0D0
        END DO

        BEND_FORCE(I, 0) = BEND_FORCE(I, 1)
        BEND_FORCE(I, NDATAN(I)) = BEND_FORCE(I, NDATAN(I)-1)
     END DO
  END IF

  !	MAKE TABLE FORCE FOR TORSION
  DO I = 1, NTTYPE
     DO J = 1, NDATT(I) - 1
        TOR_FORCE(I,J) = -(TOR_POT(I,J+1) - TOR_POT(I,J-1))/BINT(I)/D2R/2.0D0
     END DO

     TOR_FORCE(I, 0) = TOR_FORCE(I, 1)
     TOR_FORCE(I, NDATT(I)) = TOR_FORCE(I, 0)

  END DO

  !	MAKE TABLE FORCE FOR NON-BONDED INTERACTIONS

  DO I = 1, NNBTYPE
     NBOND_FORCE(1,I) = -(NBOND_POT(2,I) - NBOND_POT(1,I))/BINNB(I)/RNBOND(1,I)
     NBOND_FORCE(0,I) = NBOND_FORCE(1,I)
     DO J = 2, NDATNB(I) - 1 
        NBOND_FORCE(J,I) = -(NBOND_POT(J+1,I) - NBOND_POT(J-1,I))/BINNB(I)/RNBOND(J,I)/2.0D0
     END DO
     NBOND_FORCE(NDATNB(I),I) = 0.0D0

     !Shift forces towards end of cutoff to make smooth
!     RAMP_SIZE = 0.1 !Length of potential to shift
!     N_ENTRY = INT(RAMP_SIZE/BINNB(I)) !Number of entries to modify at tail of this table
!     DO J=1,N_ENTRY
!        SCALE = REAL(J)/N_ENTRY !Amount to shift potentials towards zero by
!        E = NDATNB(I) - (N_ENTRY - J) !Entry to modify
!        NBOND_FORCE(E,I) = NBOND_FORCE(E,I) - NBOND_FORCE(E,I)*SCALE
!     END DO
  END DO

  !	MAKE TABLE FORCE FOR IMPROPER TORSION
  DO I = 1, NOTYPE

     DO J = 1, NDATO(I) - 1
        OOP_FORCE(I,J) = -(OOP_POT(I,J+1) - OOP_POT(I,J-1))/BINO(I)/D2R/2.0D0
     END DO

     OOP_FORCE(I, 0) = OOP_FORCE(I, 1)
     OOP_FORCE(I, NDATO(I)) = OOP_FORCE(I, 0)
  END DO



  RETURN
END SUBROUTINE FTABLE

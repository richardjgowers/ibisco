!> @file
!> @brief Reads the file virtbonds
!!
!> @details virtbonds defines VS and beads which would be connected in a fully CG model, but are 
!! separated by many atoms.  This is to avoid beads which are too close from trying to interact with
!! eachother through the nonbonded forcefield.
!!
!! The file 'virtbonds' is a simple list of all the connected sites, eg:
!!
!! @verbatim
!! 1    8185
!! 1    9
!! 8185 9
!! etc
!! @endverbatim
!!
!! The index of virtual sites is the (number of virtual site) + (number of atoms), so in our example
!! 8185 refers to the first virtual site in a system with 8184 atoms.
!!

SUBROUTINE RDVIRTBONDS()

  USE VAR

  IMPLICIT NONE

  INTEGER :: IOS = 0
  INTEGER :: I, J

  !Reads the file 'virtbonds' if it exists and adds to the nonbonded exclusions
  !
  !'virtbonds' is a list of pairs of atoms that need to be exluded from nonbonded interactions
  !the index given must be the global index of the atoms

  OPEN(4, IOSTAT=IOS, FILE='virtbonds', STATUS='OLD')

  IF(IOS .NE. 0) THEN !File is not compulsory 
     CLOSE(4)
     RETURN 
  END IF
  
  DO !Read infinitely until EOF
     READ(4,*,IOSTAT=IOS) I, J
     IF(IOS .NE. 0) EXIT

     CONNECTIONS(I) = CONNECTIONS(I) + 1
     IF(CONNECTIONS(I) .GT. MAXCONNECTIONS) THEN
        WRITE(*,*) 'Max number of connections exceeded with atom ',I
        WRITE(1,*) 'Max number of connections exceeded with atom ',I
        ISTOP = 1
        RETURN
     END IF
     CONNECTED_TO(I,CONNECTIONS(I)) = J

     CONNECTIONS(J) = CONNECTIONS(J) + 1
     IF(CONNECTIONS(J) .GT. MAXCONNECTIONS) THEN
        WRITE(*,*) 'Max number of connections exceeded with atom ',J
        WRITE(1,*) 'Max number of connections exceeded with atom ',J
        ISTOP = 1
        RETURN
     END IF
     CONNECTED_TO(J,CONNECTIONS(J)) = I

  END DO

  CLOSE(4)

  RETURN

END SUBROUTINE RDVIRTBONDS

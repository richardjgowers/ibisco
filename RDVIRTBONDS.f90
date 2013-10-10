SUBROUTINE RDVIRTBONDS()

  USE VAR

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

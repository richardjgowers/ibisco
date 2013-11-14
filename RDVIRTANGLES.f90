SUBROUTINE RDVIRTANGLES()

  USE VAR

  IMPLICIT NONE

  INTEGER :: I, J, K, t_ijk

  OPEN(4,IOSTAT=IOS, FILE='virtangles', STATUS='old')

  IF(IOS .ne. 0) THEN !File not compulsory
     CLOSE(4)
     RETURN
  END IF

  DO 
     READ(4,*,IOSTAT=IOS) I, J, K, t_ijk
     IF(IOS .ne. 0) EXIT

     NIJK(I) = NIJK(I) + 1
     JANGLEIJK(I,NIJK(I)) = J
     KANGLEIJK(I,NIJK(I)) = K

  END DO

  CLOSE(4)

  RETURN

END SUBROUTINE RDVIRTANGLES

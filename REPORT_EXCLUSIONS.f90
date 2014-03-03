SUBROUTINE REPORT_EXCLUSIONS()
  ! Write a record of what nonbonded partners were excluded.
  USE VAR

  IMPLICIT NONE

  INTEGER :: I, J

  OPEN(120, FILE='exclusions.rep')

  DO I=1,NITEMS
     WRITE(120, *) I
     WRITE(120, *) (CONNECTED_TO(I,J), J=1,CONNECTIONS(I))
  END DO

  CLOSE(120)

END SUBROUTINE REPORT_EXCLUSIONS

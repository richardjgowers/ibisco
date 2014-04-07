!> @brief Debug module to record of what nonbonded partners were excluded.
!> @details Some nonbonded exclusions are detected automatically, to check that this has
!!          been done properly the debug option REPORT_EXCLUSIONS can be enabled in "ibi-preprocess"
!!          to write a list of these to a file, exclusions.rep
!> @author Rich

SUBROUTINE REPORT_EXCLUSIONS()
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

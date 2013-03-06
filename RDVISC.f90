!********************************************************************************
	SUBROUTINE RDVISC()
      USE VAR 
      USE RNEMD

      IMPLICIT NONE

       OPEN (40,IOSTAT=IOS, FILE='VISC',STATUS='OLD')

      IF (IOS.NE.0) THEN
         WRITE (1,*) ' **** FATAL ERROR! FILE VISC DOES NOT EXIST ****'
         WRITE (*,*) ' **** FATAL ERROR! FILE VISC DOES NOT EXIST ****'
         ISTOP = 1
         RETURN
          WRITE (1, *) "ISTOP?", ISTOP
      END IF
       
       READ (40,*)
       READ (40,*) NUMSLAB,NEXCH,NEMDPROF,NEMDTRAJ

      CLOSE (40)

      END SUBROUTINE

!> @file
!> @brief Reads angles between VS and beads
!> @details The file 'virtangles' allows angular potentials to be defined between a mixture of
!! beads and atoms. 
!! 
!! The format of the file virtangles is as below
!!
!! @verbatim
!! 1 8185 9 	 10 
!! 8185 9 10 	 11 
!! 9 10 8186 	 11 
!! 10 8186 18 	 12 
!! 8186 18 8187	 13 
!! 18 8187 26 	 12 
!! 8187 26 27 	 11 
!! 26 27 8188 	 11 
!! 27 8188 35 	 12 
!! 8188 35 8189	 13
!! @endverbatim
!!
!! Each line represents an angle.  The first three numbers represent the indices of the three 
!! involved particles. The index of a virtual site is (number of atoms in system) + (index of VS)
!! The last number on the line refers to the index of the angular potential, as defined in 
!! 'interaction'.

SUBROUTINE RDVIRTANGLES()

  USE VAR

  IMPLICIT NONE

  INTEGER :: IOS = 0
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

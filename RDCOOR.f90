!> @file
!> @brief Reads the file coordinate
!!
!> @details The coordinate file must be in a very specific format, detailed below
!!
!! The header needs to follow this format.  From this header, the time, box dimensions and 
!! number of molecules in the system are read
!!
!! @verbatim
!! **********************************************************************
!!  ********** Time    1.0000000000000000      (ps)
!!  ********** Box Length(X, Y, Z in nanometres)                **********
!!   5.5422830581665039        5.5422830581665039        5.5422830581665039     
!!  **********************************************************************
!!  ********** Record for each atom is in the form:-            **********
!!  ********** Index Atom_Type No._of_bonds X Y Z (coords.in nm)   *******
!!  ********** Vx Vy Vz (in nm/ps) Indices_of_bonded_atoms         *******
!!  **********************************************************************
!! num_of_molecules:           24
!! @endverbatim
!!
!! For each molecule, the following format is followed
!!
!! @verbatim
!! 341 atoms_in_molecule_no.           1  polyamide
!!     1  7  1   0.67832373D+01  -0.11701208D+02   0.20280079D+01   
!!              -0.25225833D+00   0.25166263D+00   0.16474343D+00       2
!!     2  3  4   0.69946531D+01  -0.11558007D+02   0.19942064D+01   
!!               0.38742163D+00   0.10000559D+00  -0.62110003D-01       1     3     4     5
!!     3  1  1   0.70425797D+01  -0.11580311D+02   0.20895295D+01   
!!               0.19891964D+01   0.22775315D+01  -0.28066809D+00       2
!! etc etc..
!! @endverbatim

SUBROUTINE RDCOOR()

  USE MODULEPARSING
  USE VAR

  IMPLICIT NONE

  INTEGER :: IOS
  INTEGER :: I, J, L=0, LL, K, TT
  character(len=20) :: text

  OPEN (3, IOSTAT=IOS, FILE='coordinate', STATUS='OLD')

  IF (IOS.NE.0) THEN
     WRITE (1,*) ' **** FATAL ERROR! File coordinate does not exist ****'
     WRITE (*,*) ' **** FATAL ERROR! File coordinate does not exist ****'
     ISTOP=1
     RETURN
  END IF

  READ (3, '(A80)') LINE
  CALL PARSE ()
  READ (STRNGS(1),*) TITLE

  READ (3, '(A80)') LINE
  CALL PARSE ()
  IF (STRNGS(1) == 'Time') THEN
     READ (STRNGS(2),*) INITIME
  END IF
  READ(3,*)
  READ(3,*)BOXX, BOXY, BOXZ
  BOXXINV = 1.0D0 / BOXX
  BOXYINV = 1.0D0 / BOXY
  BOXZINV = 1.0D0 / BOXZ

  BOXX2 = BOXX / 2.0D0
  BOXY2 = BOXY / 2.0D0
  BOXZ2 = BOXZ / 2.0D0

  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)text, NMOL

  ALLOCATE(NATM(NMOL))
  allocate(name_mol(nmol))

  ALLOCATE(TYPE_LABEL(NITEMS))
  TYPE_LABEL = 1 !All things are labelled as atoms by default

  DO I = 1, NMOL !101

     READ(3,*,iostat=ios)NATM(I),text,tt,name_mol(i)
     if(ios .ne. 0)then
        call error_noname()
        istop=1
        return
     end if

     DO J = 1, NATM(I)  !102
        L = L + 1
        IF ( NATOMS .lt. L ) THEN
           WRITE (1,*)' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
           WRITE (*,*)' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
           ISTOP=1
           RETURN
        END IF
        READ (3,*) LL,ITYPE(L),NBONDS(L),SX(L),SY(L),SZ(L)
        READ (3,*) VX(L),VY(L),VZ(L),(JBOND(L,K),K=1,NBONDS(L))
        !READ (3,*)(JBOND(L,K),K=1,NBONDS(L))

        MOL(L) = I !Which molecule atom L belongs to

        IF(IBRDESCR .eq. 0) THEN
           IF(name_label(itype(l)) .EQ. 'A' .OR. name_label(itype(l)) .EQ. 'a') THEN
              TYPE_LABEL(L) = 1	
           ELSE IF (name_label(itype(l)) .EQ. 'B' .OR. name_label(itype(l)) .EQ. 'b') THEN
              TYPE_LABEL(L) = 2
           ELSE
              WRITE(*,*) '**** FATAL ERROR, the label must be equal to A or B ****'
              WRITE(*,*) '**** in atom number: ',L,' *****'
              WRITE(1,*) '**** FATAL ERROR, the label must be equal to A or B ****'
              WRITE(1,*) '**** in atom number: ',L,' *****'
              ISTOP = 1			
              RETURN			
           END IF
        END IF

        VX(L) = VX(L)*1.e3 / VSCALE
        VY(L) = VY(L)*1.e3 / VSCALE 
        VZ(L) = VZ(L)*1.e3 / VSCALE 

     END DO !102 		CONTINUE
  END DO !101     CONTINUE

  IF ( NATOMS /= L ) THEN
     WRITE (1,*)' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
     WRITE (*,*)' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
     ISTOP=1
     RETURN
  END IF
  CLOSE (3)

  RETURN
END SUBROUTINE RDCOOR

subroutine error_noname
		WRITE (*,*)
		WRITE (*,*) '************************ FATAL ERROR ***************************'
		WRITE (*,*) '*           Some molecules name missing. Check your            *'
		WRITE (*,*) '*                      coordinate file                         *'
		WRITE (*,*) '****************************************************************'
		WRITE (*,*)
end subroutine

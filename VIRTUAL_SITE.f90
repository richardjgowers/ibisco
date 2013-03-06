! 	Created by Nicodemo
!	Dec 2010
!
!	This subroutine reads the index of the real atom which are virtual site also and
!	and the index of the atoms for each virtual bead
!
!
!	If VIRTSITE=0 the program compute the virtual site as the COM (Center Of Mass) of the 
!	group of atoms.
!	NVIRTA = the total number of virtual site in the system
!	NB = number of different kind of bead in the system
!	VIRTRX(or Y,Z) = coordinate of the virtual site. If VIRTSITE == 1 this represent the coordinate of 
!		a real atoms. Otherwise it represents the value of the center of mass of a group of atoms
!	COMPCOM = vector that contains the index of the atom for each bead (only for COM)
!
!	***************************************************************************
	SUBROUTINE VIRTUAL_SITE ()
	
	USE MODULEPARSING
	USE VAR

	IMPLICIT NONE
    INTEGER :: iost,I,NV,TOTATM=0,J,H, NBEAD
	INTEGER :: L
	REAL(KIND=RKIND) :: SUMTOTBMASS, SUMTOTX,SUMTOTY,SUMTOTZ
	REAL(KIND=RKIND) :: SPX,SPY,SPZ
	REAL(KIND=RKIND) :: RPX,RPY,RPZ
	CHARACTER(len=80)::TEXT_TITLE
	CHARACTER(len=10):: ENDCHECK
        
!	**************************************************************************


! NBEAD represents the maximum number of atoms in a Bed. I choose the value of 15 but it can ben increased 
! or decreased as necessary. 15 should be a good compromise value between the real occuring case and 
! saving the use of the memory. It is used only in the computation of the virtual site through a function

NBEAD = 15

OPEN(UNIT=10,FILE=name_file_virt,STATUS='old', IOSTAT=iost)

    IF (iost .NE. 0) THEN
		WRITE(*,*) '*** FATAL ERROR! File virtual site does not exist ****'
		WRITE(1,*) '*** FATAL ERROR! File virtual site does not exist ****'
		ISTOP = 1
		RETURN
	END IF

READ(10,*) TEXT_TITLE
READ(10,*)
READ(10,*) NVIRTA
READ(10,'(A80)') LINE
	CALL PARSE()
IF (STRNGS(1) == 'bead_type') THEN
		READ (STRNGS(2),*) NB
END IF
READ(10,*)

ALLOCATE(VIRTRX(NVIRTA))
ALLOCATE(VIRTRY(NVIRTA))
ALLOCATE(VIRTRZ(NVIRTA))
ALLOCATE(INDEX_VSITE(NVIRTA))
ALLOCATE(TYPE_BEAD(NB))
ALLOCATE(NUMATOM(NB))
ALLOCATE(TOTBMASS(NVIRTA))
ALLOCATE(INVTOTBMASS(NVIRTA))
ALLOCATE(VITYPE(NVIRTA))
allocate(init_numbcomp(nvirta))

TOTBMASS = 0
INVTOTBMASS = 0

DO I=1,NB
	READ (10, '(A80)') LINE
		CALL PARSE ()
		IF (STRNGS(1) .EQ. 'index') THEN
						ISTOP = 1
		WRITE(*,*) '             ******************* FATAL ERROR **********************'
		WRITE(*,*) '             * Number of bead_types different from Number of Bead *'
		WRITE(*,*) '             *                Check Virtual File                  *'
		WRITE(*,*) '             ******************************************************'
			RETURN 
		END IF
	READ (STRNGS(1),*) TYPE_BEAD(I)
	READ (STRNGS(2),*) NUMATOM(I)
!	READ(10,*)TYPE_BEAD(I),NUMATOM(I)
END DO

READ(10,*)
READ(10,*)	

 ! VS Defined as one atom

ALLOCATE(INDX_ATM(NVIRTA,NATOMS))
ALLOCATE(virtual_center(NATOMS))
INDX_ATM = 0
virtual_center = 0

	DO I=1,NVIRTA
	    READ(10,*,IOSTAT=iost)INDEX_VSITE(I),VITYPE(I)
        virtual_center(INDEX_VSITE(I)) = I
	END DO 	

READ(10,*) ENDCHECK

	DO I=1,NVIRTA
		DO J = 1,NB	
			IF(VITYPE(I) .EQ. TYPE_BEAD(J)) THEN
                init_numbcomp(i) = NUMATOM(J)
				DO H = 1,NUMATOM(J)
					READ(10,*,IOSTAT=iost)INDX_ATM(I,H)
!	If iost NE 0 means that IBIsCO reaches the end of file "virtual" but it is trying to read again, 
!	so there is a mismatch beetwen the total number of atoms in virtual beads and the number of index 
!	in this file  
					IF(iost .NE. 0) THEN
		WRITE(*,*)
		WRITE(*,*) '**** FATAL ERROR! Number of atom different from number of atom in bead ****'
		WRITE(*,*) '****                    Check virtual file                             ****'
		WRITE(*,*)
		WRITE(1,*) '**** FATAL ERROR! Number of atom different from number of atom in bead ****'
		WRITE(1,*) '****                    Check virtual file                             ****'
		WRITE(1,*)		
		ISTOP = 1
		RETURN
		END IF
				END DO
			END IF
		END DO
	END DO

!	After the reading if another reading give iost .EQ. 0 means that there are other index in the file 
!	"virtual" not readed so there is a mismatch beetwen the total number of atoms in virtual beads and 
!	the number of index  in this file  

	READ(10,*,IOSTAT=iost)
	IF(iost .EQ. 0) THEN
		WRITE(*,*)
		WRITE(*,*) '*** FATAL ERROR! Number of atom different from number of atom in bead *****'
		WRITE(*,*) '****                    Check virtual file                             ****'
		WRITE(*,*)
		WRITE(1,*) '*** FATAL ERROR! Number of atom different from number of atom in bead *****'
		WRITE(1,*) '****                    Check virtual file                             ****'
		WRITE(1,*)
		ISTOP = 1
		RETURN
	END IF

RETURN
END
!	**************************************************************************

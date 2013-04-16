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
	SUBROUTINE VIRTUAL_SITE_COM ()
	
	USE MODULEPARSING
	USE VAR

	IMPLICIT NONE
        INTEGER :: iost,I,NV,TOTATM=0,J,H, NBEAD
	INTEGER :: L,count_VS,jj,num_indx
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
ALLOCATE(virtNmol(nmol))
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

ALLOCATE(COMVECT(NVIRTA))
ALLOCATE(TYPECOM(NVIRTA))
ALLOCATE(COMPCOM(NVIRTA,NBEAD))

! Initializer for counter of VS in each molecule
 count_VS = 0
 jj=1
 num_indx = natm(1)

 DO L=1,NVIRTA

    SUMTOTBMASS = 0

    READ(10,*,IOSTAT=iost) COMPCOM(L,1), VITYPE(L)
    IF (iost .NE. 0) THEN
       ISTOP = 1
       WRITE(*,*)
       WRITE(*,*) '**** FATAL ERROR! Number of atom different from number of atom in bead ****'
       WRITE(*,*) '****                    Check virtual file                             ****'
       WRITE(*,*)
       WRITE(1,*) '**** FATAL ERROR! Number of atom different from number of atom in bead ****'
       WRITE(1,*) '****                    Check virtual file                             ****'
       WRITE(1,*)		
       RETURN
    END IF
    if(type_label(COMPCOM(L,1)) .ne. 1) then
       WRITE(*,*)
       WRITE(*,*) '**** FATAL ERROR! Molecule nr',COMPCOM(L,1), ' is a bead!!! ****'
       WRITE(*,*)
       STOP
    end if

    ! Find how many VS we have for each molecule
    if(COMPCOM(L,1) .lt. num_indx)then
       count_VS = count_VS+1
    else
       virtNmol(jj)=count_VS
       jj=jj+1
       num_indx = num_indx+natm(jj)
       count_VS = 1
    end if
    if(l .eq. nvirta)virtNmol(nmol)=count_VS
    !--------------------------------------------


    DO I = 1,NB
       IF (VITYPE(L) .EQ. TYPE_BEAD(I)) THEN
          init_numbcomp(l) = NUMATOM(I)
          DO J=2,NUMATOM(I)
             READ(10,*,IOSTAT=iost) COMPCOM(L,J)
             IF (iost .NE. 0) THEN
                ISTOP = 1
		WRITE(*,*)
		WRITE(*,*) '**** FATAL ERROR! Number of atom different from number of atom in bead ****'
		WRITE(*,*) '****                    Check virtual file                             ****'
		WRITE(*,*)
		WRITE(1,*) '**** FATAL ERROR! Number of atom different from number of atom in bead ****'
		WRITE(1,*) '****                    Check virtual file                             ****'
		WRITE(1,*)		
                RETURN
             END IF
             if(type_label(COMPCOM(L,J)) .ne. 1) then
		WRITE(*,*)
		WRITE(*,*) '**** FATAL ERROR! Molecule nr',COMPCOM(L,J), ' is a bead!!! ****'
		WRITE(*,*)
		STOP
             end if
          END DO
       END IF
    END DO
 END DO

READ(10,*,IOSTAT=iost)
IF (iost .EQ. 0) THEN
   ISTOP = 1
   WRITE(*,*)
   WRITE(*,*) '**** FATAL ERROR! Number of atom different from number of atom in bead ****'
   WRITE(*,*) '****                    Check virtual file                             ****'
   WRITE(*,*)
   WRITE(1,*) '**** FATAL ERROR! Number of atom different from number of atom in bead ****'
   WRITE(1,*) '****                    Check virtual file                             ****'
   WRITE(1,*)		
   RETURN
END IF

 DO I=1,NVIRTA
    SUMTOTBMASS = 0
    DO J=1,NB
       IF (VITYPE(I) .EQ. TYPE_BEAD(J)) THEN
          init_numbcomp(I) = NUMATOM(J)
          DO H = 1,NUMATOM(J)
             SUMTOTBMASS = SUMTOTBMASS+MASS(ITYPE(COMPCOM(I,H)))
          END DO
       END IF
    END DO
    TOTBMASS(I) = SUMTOTBMASS
    INVTOTBMASS(I) = 1/TOTBMASS(I)
!Calculate mass coefficient for virtual site
    DO H=1,init_numbcomp(I)
       masscoeff(I,H) = MASS(ITYPE(COMPCOM(I,H)))*INVTOTBMASS(I)
    END DO
 END DO

    DO I=1,NVIRTA

       SUMTOTX = 0
       SUMTOTY = 0
       SUMTOTZ = 0

       DO J=1,NB
          IF (VITYPE(I) .EQ. TYPE_BEAD(J)) THEN
             DO H = 1,NUMATOM(J)
                SUMTOTX = SUMTOTX+MASS(ITYPE(COMPCOM(I,H)))*SX(COMPCOM(I,H))
                SUMTOTY = SUMTOTY+MASS(ITYPE(COMPCOM(I,H)))*SY(COMPCOM(I,H))
                SUMTOTZ = SUMTOTZ+MASS(ITYPE(COMPCOM(I,H)))*SZ(COMPCOM(I,H))
             END DO
          END IF
       END DO

       VIRTRX(I) = SUMTOTX*INVTOTBMASS(I)
       VIRTRY(I) = SUMTOTY*INVTOTBMASS(I)
       VIRTRZ(I) = SUMTOTZ*INVTOTBMASS(I)

       IF ((VIRTRX(I)> BOXX2).OR.(VIRTRX(I) < -BOXX2).OR.(VIRTRY(I)> BOXY2) &
            .OR.(VIRTRY(I)< -BOXY2).OR.(VIRTRZ(I)> BOXZ2).OR.(VIRTRZ(I)<-BOXZ2)) THEN

          RPX = VIRTRX(I)/BOXX2
          RPY = VIRTRY(I)/BOXY2
          RPZ = VIRTRZ(I)/BOXZ2

          SPX = RPX - 2.0*INT(RPX/2.0)
          SPY = RPY - 2.0*INT(RPY/2.0)
          SPZ = RPZ - 2.0*INT(RPZ/2.0)

          RPX = (SPX - 2.0*INT(SPX))*BOXX2
          RPY = (SPY - 2.0*INT(SPY))*BOXY2
          RPZ = (SPZ - 2.0*INT(SPZ))*BOXZ2

          VIRTRX(I) = RPX
          VIRTRY(I) = RPY
          VIRTRZ(I) = RPZ 

       END IF

    END DO

RETURN
END
!	**************************************************************************

!     Modified by Nicodemo Di Pasquale
!           June 2011
!
! mol_counter = molecules counter. Counter for the number of atoms/bead in one molecule of the system
!
! init_change(i) says if that group of atoms have to be changed with the corresponding bead
!     init_change(i) = 1 ---> atoms have to be changed
!     init_change(i) = 0 ---> atoms remains atoms
!
!

	SUBROUTINE RDCOOR()

	USE MODULEPARSING
	USE VAR
        IMPLICIT NONE
!	character(20)  text
	INTEGER :: I, J, L, LL, K,q,kk,qq,p=0,jj=0, lv,TT
      integer :: shifted,ql
      integer :: nbon
      INTEGER :: PP=0, JP, NOBOUND
	INTEGER :: count_VS,num_indx
!      integer :: init=20, fin=27
      integer :: num_mol,counter=1,Btype, recover, atom_in_bead, mol_counter, &
                  num_bond,tot_virt_site=0!, original_atom=1
      integer,pointer :: init_change(:)
      real(kind=rkind) :: comrx,comry,comrz,totmass_ad, vec_pcom,x,vcomrx,vcomry,vcomrz
      real(kind=rkind) :: Fcomrx,Fcomry,Fcomrz
      real(kind=rkind) :: xf,yf,zf
      character(len=20) :: text
!      character(len=20) :: flush
!       *******************************************************************

	OPEN (3, IOSTAT=IOS, FILE='coordinate', STATUS='OLD')

	IF (IOS.NE.0) THEN
	        WRITE (1,*)	&
		   ' **** FATAL ERROR! File coordinate does not exist ****'
		WRITE (*,*)	&
		   ' **** FATAL ERROR! File coordinate does not exist ****'
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
	L = 0
	ALLOCATE(NATM(NMOL))
    allocate(name_mol(nmol))

IF (IBRDESCR .EQ. 0) THEN

            ALLOCATE(TYPE_LABEL(NATOMS))

	TYPE_LABEL = 0

	DO 101 I = 1, NMOL

		READ(3,*,iostat=ios)NATM(I),text,tt,name_mol(i)
        if(ios .ne. 0)then
            call error_noname()
            istop=1
            return
        end if	
    

		DO 102 J = 1, NATM(I)
            	L = L + 1
	IF ( NATOMS .lt. L ) THEN
	        WRITE (1,*)	&
		   ' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
	        WRITE (*,*)	&
		   ' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
        	ISTOP=1
         	RETURN
      	END IF
            	READ (3,*) LL,ITYPE(L),NBONDS(L),SX(L),SY(L),SZ(L)
            	READ (3,*) VX(L),VY(L),VZ(L),(JBOND(L,K),K=1,NBONDS(L))
		!READ (3,*)(JBOND(L,K),K=1,NBONDS(L))

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

		VX(L) = VX(L)*1.e3 / VSCALE
		VY(L) = VY(L)*1.e3 / VSCALE 
		VZ(L) = VZ(L)*1.e3 / VSCALE 

102 		CONTINUE
101     CONTINUE


    if(virtsite .eq. 1)then
ALLOCATE(virtNmol(nmol))
! Find how many VS we have for each molecule
! Initializer for counter of VS in each molecule
 count_VS = 0
 jj=1
 num_indx = natm(1)

	DO I=1,NVIRTA
        if(INDEX_VSITE(I) .lt. num_indx)then
            count_VS = count_VS+1
        else
            virtNmol(jj)=count_VS
            jj=jj+1
            num_indx = num_indx+natm(jj)
            count_VS = 1
        end if
        if(l .eq. nvirta)virtNmol(nmol)=count_VS
	END DO 
end if

!--------------------------------------------

! Divide beads and atoms in two different part to be processed separately

NUM_BEAD = 0
NUM_VS   = 0
NUM_BA  = 0

if(MTS_CHECK .eq. 0)then
    if(virtsite .eq. 0)then

! From 1 to NUM_BEAD we have Bead
! From NUM_BEAD+1 to NUM_BA we have beads connected to atoms (not considered in MTS) 
! From NUM_BA + 1 to NATOMS we have atoms 

        DO I = 1,NATOMS
            NOBOUND = 0
            IF(TYPE_LABEL(I) .EQ. 2) THEN
                DO JP = 1,NBONDS(I)
                   IF(TYPE_LABEL(JBOND(I,JP)) .EQ. 1) THEN
                      NOBOUND = 1                              
                   END IF
                END DO
               IF(NOBOUND .EQ. 0) THEN
                  NUM_BEAD = NUM_BEAD+1
                  INDEX_AB(NUM_BEAD) = I
               END IF
            END IF
        END DO

        PP = NUM_BEAD

        DO I = 1,NATOMS
            NOBOUND = 0
            IF (TYPE_LABEL(I) .EQ. 2) THEN
                DO JP = 1,NBONDS(I)
                    IF(TYPE_LABEL(JBOND(I,JP)) .EQ. 1) THEN
                        NOBOUND = 1                              
                    END IF
                END DO
                IF(NOBOUND .EQ. 1) THEN
                    PP = PP+1
                    INDEX_AB(PP) = I
                END IF
            END IF
        END DO

        NUM_BA = PP
    
        DO I = 1,NATOMS
            IF (TYPE_LABEL(I) .EQ. 1) THEN
                PP = PP + 1
                INDEX_AB(PP) = I
            END IF
        END DO
    else

! From 1 to NUM_BEAD we have Bead
! From NUM_BEAD+1 to NUM_BA we have beads connected to atoms (not considered in MTS) 
! From NUM_BA + 1 to NUM_VS we have Virtual Site
! From NUM_VS+1 to NATOMS we have atoms (no VS) 

        DO I = 1,NATOMS
            NOBOUND = 0
            IF(TYPE_LABEL(I) .EQ. 2) THEN
                DO JP = 1,NBONDS(I)
                   IF(TYPE_LABEL(JBOND(I,JP)) .EQ. 1) THEN
                      NOBOUND = 1                              
                   END IF
                END DO
               IF(NOBOUND .EQ. 0) THEN
                  NUM_BEAD = NUM_BEAD+1
                  INDEX_AB(NUM_BEAD) = I
               END IF
            END IF
        END DO

        NUM_BA = NUM_BEAD

        DO I = 1,NATOMS
            NOBOUND = 0
            IF (TYPE_LABEL(I) .EQ. 2) THEN
                DO JP = 1,NBONDS(I)
                    IF(TYPE_LABEL(JBOND(I,JP)) .EQ. 1) THEN
                        NOBOUND = 1                              
                    END IF
                END DO
                IF(NOBOUND .EQ. 1) THEN
                    NUM_BA = NUM_BA+1
                    INDEX_AB(NUM_BA) = I
                END IF
            END IF
        END DO

        NUM_VS = NUM_BA

        DO I = 1,NATOMS
            IF (TYPE_LABEL(I) .EQ. 1 .and.  virtual_center(I) .ne. 0) THEN
                NUM_VS = NUM_VS + 1
                INDEX_AB(NUM_VS) = I
            END IF
        END DO

        PP = NUM_VS

        DO I = 1,NATOMS
            IF (TYPE_LABEL(I) .EQ. 1  .and.  virtual_center(I) .eq. 0) THEN
                PP = PP + 1
                INDEX_AB(PP) = I
            END IF
        END DO

    end if !if(virtsite .eq. 0)
else

    if(virtsite .eq. 0)then

! From 1 to NUM_BEAD we have Bead
! From NUM_BEAD+1 to NATOMS we have atoms 

        DO I = 1,NATOMS
            IF(TYPE_LABEL(I) .EQ. 2) THEN
                NUM_BEAD = NUM_BEAD+1
                INDEX_AB(NUM_BEAD) = I
            END IF
        END DO
    
        PP = NUM_BEAD
    
        DO I = 1,NATOMS
            IF (TYPE_LABEL(I) .EQ. 1) THEN
                PP = PP + 1
                INDEX_AB(PP) = I
            END IF
        END DO

    else

! From 1 to NUM_BEAD we have Bead
! From NUM_BEAD+1 to NUM_VS we have virtual site 
! From NUM_VS+1 to NATOMS we have atoms 

        DO I = 1,NATOMS
            IF(TYPE_LABEL(I) .EQ. 2) THEN
                NUM_BEAD = NUM_BEAD+1
                INDEX_AB(NUM_BEAD) = I
            END IF
        END DO
        NUM_VS = NUM_BEAD
        DO I = 1,NATOMS
            IF (TYPE_LABEL(I) .EQ. 1 .and.  virtual_center(I) .ne. 0) THEN
                NUM_VS = NUM_VS + 1
                INDEX_AB(NUM_VS) = I
            END IF
        END DO
        PP = NUM_VS
        DO I = 1,NATOMS
            IF (TYPE_LABEL(I) .EQ. 1  .and.  virtual_center(I) .eq. 0) THEN
                PP = PP + 1
                INDEX_AB(PP) = I
            END IF
        END DO

    end if ! VIRTSITE
end if ! MTS


	IF ( NATOMS /= L ) THEN
	        WRITE (1,*)	&
		   ' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
	        WRITE (*,*)	&
		   ' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
        	ISTOP=1
         	RETURN
      	END IF
        CLOSE (3)



ELSE

	DO I = 1, NMOL

		READ(3,*,iostat=ios)NATM(I),text,tt,name_mol(i)
        if(ios .ne. 0)then
            call error_noname()
            istop=1
            return
        end if	
		
		DO J = 1, NATM(I)
            	L = L + 1
            	READ (3,*) LL,ITYPE(L),NBONDS(L),SX(L),SY(L),SZ(L)
            	READ (3,*) VX(L),VY(L),VZ(L),(JBOND(L,K),K=1,NBONDS(L))
		!READ (3,*)(JBOND(L,K),K=1,NBONDS(L))

		VX(L) = VX(L)*1.e3 / VSCALE
		VY(L) = VY(L)*1.e3 / VSCALE 
		VZ(L) = VZ(L)*1.e3 / VSCALE 

 	END DO
     END DO

	IF ( NATOMS /= L ) THEN
	        WRITE (1,*)	&
		   ' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
	        WRITE (*,*)	&
		   ' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
        	ISTOP=1
         	RETURN
      	END IF
        CLOSE (3)

END IF

	RETURN
	END

subroutine error_noname

		WRITE (*,*)
		WRITE (*,*) '************************ FATAL ERROR ***************************'
		WRITE (*,*) '*           Some molecules name missing. Check your            *'
		WRITE (*,*) '*                      coordinate file                         *'
		WRITE (*,*) '****************************************************************'
		WRITE (*,*)

end subroutine

!	*********************************************************************

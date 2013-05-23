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

  INTEGER :: count_VS,num_indx
  !      integer :: init=20, fin=27
  integer :: num_mol,counter=1,Btype, recover, atom_in_bead, mol_counter, &
       num_bond,tot_virt_site=0!, original_atom=1
  integer,pointer :: init_change(:)
  real(kind=rkind) :: comrx,comry,comrz,totmass_ad, vec_pcom,x,vcomrx,vcomry,vcomrz
  real(kind=rkind) :: Fcomrx,Fcomry,Fcomrz
  real(kind=rkind) :: xf,yf,zf
  character(len=20) :: text

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
     WRITE(*,*) INITIME
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

  IF (IBRDESCR .EQ. 0) THEN !Hybrid RDCOOR
     ALLOCATE(TYPE_LABEL(NITEMS))
  END IF

  TYPE_LABEL = 0

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
     WRITE (1,*)	&
          ' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
     WRITE (*,*)	&
          ' **** FATAL ERROR! No. Atoms in coordinate does not equal to NATOMS ****'
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

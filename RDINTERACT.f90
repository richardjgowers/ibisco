SUBROUTINE RDINTERACT()
  USE MODULEPARSING
  USE VAR

  IMPLICIT NONE

  INTEGER		I, J, L, LL, IB, JB, IA, JA, KA, IT, JT, KT, LT, K,ITA,h
  INTEGER		INB, JNB, TYPEI
  INTEGER :: iosIN = 0
  REAL*8		R
  CHARACTER(len=80) :: LINE2
  logical :: alloc=.true.

  OPEN (4, IOSTAT=IOS, FILE='interaction', STATUS='OLD')

  IF (IOS.NE.0) THEN
     WRITE (1,*)	&
          ' **** FATAL ERROR! File interaction does not exist ****'
     WRITE (*,*)	&
          ' **** FATAL ERROR! File interaction does not exist ****'
     ISTOP=1
     RETURN
  END IF

  !*******************************************************
  !****************reading atoms parameters***************
  !*******************************************************
  !	READ(4,*)
  READ (4, '(A80)') LINE
  CALL PARSE ()
  IF (STRNGS(1) == 'atom_types') THEN	
     READ (STRNGS(2),*) NTYPE	
  END IF

  ALLOCATE(LABEL(NTYPE))
  ALLOCATE(MASS0(NTYPE))
  ALLOCATE(MASS(NTYPE))
  ALLOCATE(INVMASS(NTYPE))

  READ(4,*)
  DO I = 1, NTYPE
     READ (4, '(A80)') LINE
     CALL PARSE ()
     IF (STRNGS(1) .EQ. 'bond_types') THEN
        ISTOP = 1
        WRITE(*,*) '             ******************* FATAL ERROR **********************'
        WRITE(*,*) '             * Number of atom_types different from Number of Atom *'
        WRITE(*,*) '             *                Check Interaction File              *'
        WRITE(*,*) '             ******************************************************'
        RETURN 
     END IF
     if(ibrdescr .eq. 0)then
        READ (STRNGS(1),*) TYPEI
        READ (STRNGS(2),*) LABEL(TYPEI)
        READ (STRNGS(3),*) MASS0(TYPEI)
        READ (STRNGS(4),*,iostat=ios) name_label(TYPEI)
        if(ios .ne. 0)then
           call error_inter ()
           ISTOP = 1
           return
        end if
        if(name_label(TYPEI) .ne. 'A' .and. name_label(TYPEI) .ne. 'a')then
           if(name_label(TYPEI) .ne. 'B' .and. name_label(TYPEI) .ne. 'b')then  
              call error_inter ()
              ISTOP = 1
              return
           end if
        end if
     else
        READ (STRNGS(1),*) TYPEI
        READ (STRNGS(2),*) LABEL(TYPEI)
        READ (STRNGS(3),*) MASS0(TYPEI)
     end if
     MASS(I) = MASS0(I)/NA/MASSSCALE!/1000.0
     INVMASS(I) = 1.0D0 / MASS(I)
  END DO

  ALLOCATE(IBONDT(NTYPE, NTYPE ))
  ALLOCATE(INBONDT(NTYPE, NTYPE ))
  !**** Zero Pointer Arrays for Bonds, Bends and Torsions ****
  IBONDT = 0
  INBONDT = 0 
  !hj	IANGT = 0
  !hj	ITORT = 0
  !hj	IOOPT = 0
  !*******************************************************
  !****************reading bonds parameters***************
  !*******************************************************

  READ (4, '(A80)') LINE
  CALL PARSE ()
  IF (STRNGS(1) == 'bond_types') THEN
     READ (STRNGS(2),*) NBTYPE
  ELSE
     ISTOP = 1
     WRITE(*,*) '             ******************* FATAL ERROR **********************'
     WRITE(*,*) '             * Number of atom_types different from Number of Atom *'
     WRITE(*,*) '             *                Check Interaction File              *'
     WRITE(*,*) '             ******************************************************'
     RETURN 
  END IF
  READ(4,*)
  IF (NBTYPE.GT.0) THEN
     ALLOCATE(RBOND(NBTYPE,0:MAXINPUT))
     ALLOCATE(BOND_FORCE(NBTYPE,0:MAXINPUT))
     ALLOCATE(BOND_POT(NBTYPE,0:MAXINPUT))
     ALLOCATE(BINB(NBTYPE))
     ALLOCATE(NDATB(NBTYPE))
     allocate(typeBond(nbtype))
     typeBond = .false.
     BOND_FORCE = 0.0D0
     BOND_POT = 0.0D0

     h=0
     DO I = 1, NBTYPE
        READ (4, '(A80)') LINE
        CALL PARSE ()
        IF (STRNGS(1) .EQ. 'angle_types') THEN
           ISTOP = 1
           WRITE(*,*) '             ******************* FATAL ERROR **********************'
           WRITE(*,*) '             * Number of bond_types different from Number of Bond *'
           WRITE(*,*) '             *                Check Interaction File              *'
           WRITE(*,*) '             ******************************************************'
           RETURN 
        END IF
        READ (STRNGS(1),*) IB
        READ (STRNGS(2),*) JB
        !write(*,*)STRNGS(1),STRNGS(2),STRNGS(3),STRNGS(4)

        IF(IBONDT(IB,JB).EQ.0) THEN
           IBONDT(IB,JB) = I
           IBONDT(JB,IB) = I
        ELSE
           WRITE (1,*)	&
                '**** FATAL ERROR: Duplicate entry in interaction file for bonds ****'
           WRITE (*,*)	&
                '**** FATAL ERROR: Duplicate entry in interaction file for bonds ****'
           ISTOP=1
           RETURN
        ENDIF

        ! If you choose to use table for potentials INTERACT==1, else if you want to use gaussian function INTERACT==0

        IF (INTERACT == 1) THEN

           if(STRNGS(3) .eq. 'constraint')then

              !   If there are some constraints, it allocates the corresponding vectors        
              if(alloc)then
                 allocate(constr(nbtype))
                 alloc = .false.
              end if
              typeBond(i) = .true.
              shakeOK = .true.
              READ (STRNGS(4),*) constr(i)
              h=h+1
           else

              shakeOK=.false.
              OPEN (11, IOSTAT=IOS, FILE=STRNGS(3), STATUS='OLD')
              IF (IOS.NE.0) THEN
                 WRITE (1,*)	&
                      ' **** FATAL ERROR! File ', STRNGS(3),' does not exist ****'
                 WRITE (*,*)	&
                      ' **** FATAL ERROR! File ', STRNGS(3),' does not exist ****'
                 ISTOP=1
                 RETURN
              END IF

              K = 0
              DO WHILE (.TRUE.)
                 READ(11,*,IOSTAT=IOS2) RBOND(I,K), BOND_POT(I,K)
                 IF(IOS2 .ne. 0) EXIT
                 K = K+1
              END DO
              NDATB(I) = K - 1

              CLOSE(11)

              !Rescale energies
              DO J =0,NDATB(I)
                 BOND_POT(I,J) = BOND_POT(I,J) *1000.0 /NA /ESCALE
              END DO

              BINB(I) = RBOND(I,1)-RBOND(I,0)

           END IF

        end if

     END DO

     if(h .eq. NBTYPE)then
        NoBond=.true. ! if there are only constraints there is no need for the subroutine FPBOND.f90
     else
        NoBond=.false. 
     end if


  END IF !If NBTYPE gt 0

  !**********************************
  !**** Angle-bending parameters ****
  !**********************************

  READ (4, '(A80)') LINE
  CALL PARSE ()
  IF (STRNGS(1) == 'angle_types') THEN
     READ (STRNGS(2),*) NATYPE
  ELSE
     ISTOP = 1
     WRITE(*,*) '             ******************* FATAL ERROR **********************'
     WRITE(*,*) '             * Number of bond_types different from Number of Bond *'
     WRITE(*,*) '             *                Check Interaction File              *'
     WRITE(*,*) '             ******************************************************'
     RETURN 
  END IF
  READ(4,*)
  IF (NATYPE.GT.0) THEN
     ALLOCATE(ANGLE(NATYPE, 0:MAXINPUT))
     ALLOCATE(BEND_FORCE(NATYPE, 0:MAXINPUT))
     ALLOCATE(BEND_POT(NATYPE, 0:MAXINPUT))
     ALLOCATE(BINA(NATYPE))
     ALLOCATE(NDATAN(NATYPE))

     ALLOCATE(JANGLEIJK(NITEMS,10))
     ALLOCATE(KANGLEIJK(NITEMS,10))

     ALLOCATE(NOJANGLEIJK(NITEMS,10))
     ALLOCATE(NOKANGLEIJK(NITEMS,10))

     ALLOCATE(IANGT(NTYPE,NTYPE,NTYPE))
     ALLOCATE(ATIANG(NATYPE*3))

     IANGT = 0

     BEND_FORCE = 0.0D0
     BEND_POT = 0.0D0

     DO I = 1, NATYPE
        ITA = 3*(I-1)
        READ (4, '(A80)') LINE
        CALL PARSE ()
        IF (STRNGS(1) .EQ. 'torsion_types') THEN
           ISTOP = 1
           WRITE(*,*) '             ******************* FATAL ERROR ************************'
           WRITE(*,*) '             * Number of angle_types different from Number of Angle *'
           WRITE(*,*) '             *                 Check Interaction File               *'
           WRITE(*,*) '             ********************************************************'
           RETURN 
        END IF
        READ (STRNGS(1),*) IA
        READ (STRNGS(2),*) JA
        READ (STRNGS(3),*) KA

        IF(IANGT(IA,JA,KA).EQ.0) THEN
           IANGT(IA,JA,KA) = I
           IANGT(KA,JA,IA) = I
           ATIANG(ITA+1) = IA
           ATIANG(ITA+2) = JA
           ATIANG(ITA+3) = KA
        ELSE
           WRITE (1,*)	&
                '**** FATAL ERROR: Duplicate entry in interaction file for angles ****'
           WRITE (*,*)	&
                '**** FATAL ERROR: Duplicate entry in interaction file for angles ****'
           ISTOP=1
           RETURN
        ENDIF

        IF (INTERACT == 1) THEN

           OPEN (11, IOSTAT=IOS, FILE=STRNGS(4), STATUS='OLD')
           IF (IOS.NE.0) THEN
              WRITE (1,*)	&
                   ' **** FATAL ERROR! File ', STRNGS(4),' does not exist ****'

              WRITE (*,*)	&
                   ' **** FATAL ERROR! File ', STRNGS(4),' does not exist ****'
              ISTOP=1
              RETURN
           END IF

           READ(11,*,IOSTAT=IOS2)ANGLE(I,0), BEND_POT(I,0)
           READ(11,*,IOSTAT=IOS2)ANGLE(I,1), BEND_POT(I,1)
           BEND_POT(I,0) = BEND_POT(I,0)*1000.0/ NA / ESCALE
           BEND_POT(I,1) = BEND_POT(I,1)*1000.0/ NA / ESCALE

           BINA(I) = ANGLE(I,1)-ANGLE(I,0)

           IF(ANGLE(I, 0) /= 0.0) THEN
              ANGLE(I,0) = 0.0	
              ANGLE(I,2) = ANGLE(I,1)
              BEND_POT(I,2) = BEND_POT(I,1)
              ANGLE(I,1) = BINA(I)
              BEND_POT(I,1) = BEND_POT(I,0)
              K = 3
           ELSE
              K = 2
           END IF

           DO WHILE (.TRUE.)
              READ(11,*,IOSTAT=IOS2)ANGLE(I, K), BEND_POT(I,K)
              BEND_POT(I,K) = BEND_POT(I,K)*1000.0/ NA / ESCALE
              IF (IOS2 /= 0) EXIT
              K = K + 1 
           END DO
           CLOSE (11)
           NDATAN(I) = K - 1
        END IF

     END DO
  END IF
  !**********************************
  !**** Torsion angle parameters ****
  !**********************************

  READ (4, '(A80)') LINE
  CALL PARSE ()
  IF (STRNGS(1) == 'torsion_types') THEN
     READ (STRNGS(2),*) NTTYPE
  ELSE
     ISTOP = 1
     WRITE(*,*) '             ******************** FATAL ERROR ***********************'
     WRITE(*,*) '             * Number of angle_types different from Number of Angle *'
     WRITE(*,*) '             *                 Check Interaction File               *'
     WRITE(*,*) '             ********************************************************'
     RETURN 
  END IF
  READ(4,*)

  IF (NTTYPE.GT.0) THEN
     ALLOCATE(ANGLE_TOR(NTTYPE, 0:MAXINPUT))
     ALLOCATE(TOR_FORCE(NTTYPE, 0:MAXINPUT))
     ALLOCATE(TOR_POT(NTTYPE, 0:MAXINPUT))
     ALLOCATE(BINT(NTTYPE))
     ALLOCATE(NDATT(NTTYPE))

     ALLOCATE(JTORIJKL(NATOMS,10))
     ALLOCATE(KTORIJKL(NATOMS,10))
     ALLOCATE(LTORIJKL(NATOMS,10))

     ALLOCATE(FJTORIJKL(NATOMS,10))
     ALLOCATE(FKTORIJKL(NATOMS,10))
     ALLOCATE(FLTORIJKL(NATOMS,10))

     ALLOCATE(ITORT(NTYPE,NTYPE,NTYPE,NTYPE))
     ITORT = 0
     TOR_FORCE = 0.0D0
     TOR_POT = 0.0D0

     DO I = 1, NTTYPE
        READ (4, '(A80)') LINE
        CALL PARSE ()
        IF (STRNGS(1) .EQ. 'non_bonded_interaction_types') THEN
           ISTOP = 1
           WRITE(*,*) '             ********************* FATAL ERROR **************************'
           WRITE(*,*) '             * Number of torsion_types different from Number of Torsion *'
           WRITE(*,*) '             *                   Check Interaction File                 *'
           WRITE(*,*) '             ************************************************************'
           RETURN 
        END IF
        READ (STRNGS(1),*) IT
        READ (STRNGS(2),*) JT
        READ (STRNGS(3),*) KT
        READ (STRNGS(4),*) LT

        IF(ITORT(IT,JT,KT,LT).EQ.0) THEN
           ITORT(IT,JT,KT,LT) = I
           ITORT(LT,KT,JT,IT) = I
        ELSE
           WRITE (1,*)	&
                '**** FATAL ERROR: Duplicate entry in interaction file for torsions  ****'
           WRITE (*,*)	&
                '**** FATAL ERROR: Duplicate entry in interaction file for torsions  ****'
           WRITE (*,*)	&
                IT,JT,KT,LT
           ISTOP=1
           RETURN
        ENDIF

        OPEN (11, IOSTAT=IOS, FILE=STRNGS(5), STATUS='OLD')
        IF (IOS.NE.0) THEN
           WRITE (1,*)	&
                ' **** FATAL ERROR! File ', STRNGS(5),' does not exist ****'

           WRITE (*,*)	&
                ' **** FATAL ERROR! File ', STRNGS(5),' does not exist ****'
           ISTOP=1
           RETURN
        END IF

        READ(11,*,IOSTAT=IOS2)ANGLE_TOR(I,0), TOR_POT(I,0)
        READ(11,*,IOSTAT=IOS2)ANGLE_TOR(I,1), TOR_POT(I,1)
        TOR_POT(I,0) = TOR_POT(I,0)*1000.0/ NA / ESCALE
        TOR_POT(I,1) = TOR_POT(I,1)*1000.0/ NA / ESCALE

        BINT(I) = ANGLE_TOR(I,1)-ANGLE_TOR(I,0)

        IF(ANGLE_TOR(I, 0) /= 0.0) THEN
           ANGLE_TOR(I,0) = 0.0	
           ANGLE_TOR(I,2) = ANGLE_TOR(I,1)
           TOR_POT(I,2) = TOR_POT(I,1)
           ANGLE_TOR(I,1) = BINT(I)
           TOR_POT(I,1) = TOR_POT(I,0)
           K = 3
        ELSE
           K = 2
        END IF

        DO WHILE (.TRUE.)
           READ(11,*,IOSTAT=IOS2)ANGLE_TOR(I, K), TOR_POT(I,K)
           TOR_POT(I,K) = TOR_POT(I,K)*1000.0/ NA / ESCALE
           IF (IOS2 /= 0) EXIT
           K = K + 1 
        END DO
        NDATT(I) = K - 1
        CLOSE (11)


     END DO
  END IF

  !**********************************
  !**** Non Bonded parameters *******
  !**********************************

  READ (4, '(A80)') LINE
  CALL PARSE ()
  IF (STRNGS(1) == 'non_bonded_interaction_types') THEN
     READ (STRNGS(2),*) NNBTYPE
  ELSE
     ISTOP = 1
     WRITE(*,*) '           ********************* FATAL ERROR **************************'
     WRITE(*,*) '           * Number of torsion_types different from Number of Torsion *'
     WRITE(*,*) '           *                   Check Interaction File                 *'
     WRITE(*,*) '           ************************************************************'
     RETURN
  END IF
  READ(4,*)
  IF (NNBTYPE.GT.0) THEN
     ALLOCATE(RNBOND(0:MAXINPUT, NNBTYPE))
     ALLOCATE(NBOND_FORCE(0:MAXINPUT,NNBTYPE))
     ALLOCATE(NBOND_POT(0:MAXINPUT,NNBTYPE))
     ALLOCATE(BINNB(NNBTYPE))
     ALLOCATE(NDATNB(NNBTYPE))
     NBOND_FORCE = 0.0D0
     NBOND_POT = 0.0D0

     DO I = 1, NNBTYPE

        READ (4, '(A80)') LINE
        CALL PARSE ()
        IF (STRNGS(1) .EQ. 'Out_of_Planes') THEN
           ISTOP = 1
           WRITE(*,*) '******************************** FATAL ERROR **********************************'
           WRITE(*,*) '* Number of non_bonded_interaction_types different from Number NB interaction *'
           WRITE(*,*) '*                          Check Interaction File                             *'
           WRITE(*,*) '*******************************************************************************'
           RETURN 
        END IF
        READ (STRNGS(1),*) INB
        READ (STRNGS(2),*) JNB

        IF(INBONDT(INB,JNB).EQ.0) THEN
           INBONDT(INB,JNB) = I
           INBONDT(JNB,INB) = I
        ELSE
           WRITE (1,*)	&
                '**** FATAL ERROR: Duplicate entry in interaction file for non-bonded interactions ****'
           WRITE (*,*)	&
                '**** FATAL ERROR: Duplicate entry in interaction file for non-bonded interactions ****'
           WRITE (*,*)	&
                IT,JT
           ISTOP=1
           RETURN
        ENDIF

        OPEN (11, IOSTAT=IOS, FILE=STRNGS(3), STATUS='OLD')
        IF (IOS.NE.0) THEN
           WRITE (1,*)	&
                ' **** FATAL ERROR! File ', STRNGS(3),' does not exist ****'

           WRITE (*,*)	&
                ' **** FATAL ERROR! File ', STRNGS(3),' does not exist ****'
           ISTOP=1
           RETURN
        END IF

        READ(11,*,IOSTAT=IOS2)RNBOND(0, I), NBOND_POT(0,I)
        READ(11,*,IOSTAT=IOS2)RNBOND(1, I), NBOND_POT(1,I)
        NBOND_POT(0,I) = NBOND_POT(0,I)*1000.0/ NA / ESCALE
        NBOND_POT(1,I) = NBOND_POT(1,I)*1000.0/ NA / ESCALE

        BINNB(I) = RNBOND(1,I) - RNBOND(0,I)

        IF(RNBOND(0,I) /= 0.0) THEN
           RNBOND(0,I) = 0.0	
           RNBOND(2,I) = RNBOND(1,I)
           NBOND_POT(2,I) = NBOND_POT(1,I)
           RNBOND(1,I) = BINNB(I)
           NBOND_POT(1,I) = NBOND_POT(0,I)
           K = 3
        ELSE
           K = 2
        END IF

        DO WHILE (.TRUE.)
           READ(11,*,IOSTAT=IOS2)RNBOND(K,I), NBOND_POT(K,I)
           IF (RNBOND(K,I) .NE. RNBOND(K-1,I) .AND. IOS2 == 0) THEN
              NBOND_POT(K,I) = NBOND_POT(K,I)*1000.0/ NA / ESCALE
              IF (IOS2 /= 0) EXIT
              K = K + 1 
           ELSE
              EXIT
           END IF
        END DO

        NDATNB(I) = K - 1
        CLOSE (11)

     END DO
  END IF

  !**********************************
  !**** Out of Planes parameters ****
  !**********************************

  READ (4, '(A80)') LINE
  CALL PARSE ()
  IF (STRNGS(1) == 'Out_of_Planes') THEN
     READ (STRNGS(2),*) NOTYPE
  ELSE
     ISTOP = 1
     WRITE(*,*) '******************************** FATAL ERROR **********************************'
     WRITE(*,*) '* Number of non_bonded_interaction_types different from Number NB interaction *'
     WRITE(*,*) '*                          Check Interaction File                             *'
     WRITE(*,*) '*******************************************************************************'
     RETURN
  END IF
  READ(4,*)

  IF (NOTYPE.GT.0) THEN
     ALLOCATE(ANGLE_OOP(NOTYPE, 0:MAXINPUT))
     ALLOCATE(OOP_FORCE(NOTYPE, 0:MAXINPUT))
     ALLOCATE(OOP_POT(NOTYPE, 0:MAXINPUT))
     ALLOCATE(BINO(NOTYPE))
     ALLOCATE(NDATO(NOTYPE))

     ALLOCATE(JOOPIJKL(NATOMS,20))
     ALLOCATE(KOOPIJKL(NATOMS,20))
     ALLOCATE(LOOPIJKL(NATOMS,20))

     ALLOCATE(IOOPT(NTYPE,NTYPE,NTYPE,NTYPE))
     IOOPT = 0
     OOP_FORCE = 0.0D0
     OOP_POT = 0.0D0

     DO I = 1, NOTYPE
        READ (4, '(A80)',IOSTAT = iosIN) LINE
        IF (iosIN .NE. 0) THEN
           ISTOP = 1
           WRITE(*,*) '**************************** FATAL ERROR **********************************'
           WRITE(*,*) '* Number of Out_of_Planes different from Number Out of Planes interaction *'
           WRITE(*,*) '*                           Check Interaction File                        *'
           WRITE(*,*) '***************************************************************************'
           RETURN 
        END IF
        CALL PARSE ()
        READ (STRNGS(1),*) IT
        READ (STRNGS(2),*) JT
        READ (STRNGS(3),*) KT
        READ (STRNGS(4),*) LT

        IF (I .EQ. NOTYPE) THEN
           READ (4, '(A80)',IOSTAT = iosIN) LINE2
           IF(iosIN .EQ. 0) THEN
              ISTOP = 1
              WRITE(*,*) '**************************** FATAL ERROR **********************************'
              WRITE(*,*) '* Number of Out_of_Planes different from Number Out of Planes interaction *'
              WRITE(*,*) '*                           Check Interaction File                        *'
              WRITE(*,*) '***************************************************************************'
              RETURN 
           END IF
        END IF

        IF(IOOPT(IT,JT,KT,LT).EQ.0) THEN
           !                 IOOPT(IT,JT,KT,LT) = I
           IOOPT(IT,JT,KT,LT) = I
           !                  IOOPT(LT,KT,JT,IT) = I
           !               IOOPT(IT,JT,LT,KT) = I
           !		IOOPT(IT,LT,JT,KT) = I
           !		IOOPT(IT,LT,KT,JT) = I
           !		IOOPT(IT,KT,JT,LT) = I
           !		IOOPT(IT,KT,LT,JT) = I
        ELSE
           WRITE (1,*)	&
                '**** FATAL ERROR: Duplicate entry in interaction file for out of planes  ****'
           WRITE (*,*)	&
                '**** FATAL ERROR: Duplicate entry in interaction file for  out of planes   ****'
           WRITE (*,*)	&
                IT,JT,KT,LT
           ISTOP=1
           RETURN
        ENDIF

        OPEN (11, IOSTAT=IOS, FILE=STRNGS(5), STATUS='OLD')
        IF (IOS.NE.0) THEN
           WRITE (1,*)	&
                ' **** FATAL ERROR! File', STRNGS(5),' does not exist ****'

           WRITE (*,*)	&
                ' **** FATAL ERROR! File', STRNGS(5),' does not exist ****'
           ISTOP=1
           RETURN
        END IF


        READ(11,*,IOSTAT=IOS2)ANGLE_OOP(I,0), OOP_POT(I,0)
        READ(11,*,IOSTAT=IOS2)ANGLE_OOP(I,1), OOP_POT(I,1)
        OOP_POT(I,0) = OOP_POT(I,0)*1000.0/ NA / ESCALE
        OOP_POT(I,1) = OOP_POT(I,1)*1000.0/ NA / ESCALE

        if(ANGLE_OOP(I,1) .lt. ANGLE_OOP(I,0))then
           BINO(I) = ANGLE_OOP(I,0)-ANGLE_OOP(I,1)
        else
           BINO(I) = ANGLE_OOP(I,1)-ANGLE_OOP(I,0)
        end if

        IF(ANGLE_OOP(I, 0) /= 0.0) THEN
           ANGLE_OOP(I,0) = 0.0	
           ANGLE_OOP(I,2) = ANGLE_OOP(I,1)
           OOP_POT(I,2) = OOP_POT(I,1)
           ANGLE_OOP(I,1) = BINO(I)
           OOP_POT(I,1) = OOP_POT(I,0)
           K = 3
        ELSE
           K = 2
        END IF

        DO WHILE (.TRUE.)
           READ(11,*,IOSTAT=IOS2)ANGLE_OOP(I, K), OOP_POT(I,K)
           OOP_POT(I,K) = OOP_POT(I,K) *1000.0/ NA / ESCALE
           IF (IOS2 /= 0) EXIT
           K = K + 1 
        END DO
        NDATO(I) = K - 1
        CLOSE (11)


     END DO


  END IF

  CLOSE (4)
  RETURN
END SUBROUTINE RDINTERACT

subroutine error_inter

  WRITE (*,*)
  WRITE (*,*) '************************ FATAL ERROR ***************************'
  WRITE (*,*) '*           the label for hybrid model must be equal to:       *'
  WRITE (*,*) '*           A or a:  atoms                                     *'
  WRITE (*,*) '*           B or b:  beads                                     *'
  WRITE (*,*) '****************************************************************'
  WRITE (*,*)

end subroutine error_inter

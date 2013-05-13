
!	update 16/1/2008 15:57:50
!	********************************************************************************
SUBROUTINE SETLIS()

  USE VAR
  IMPLICIT NONE
  INTEGER		N, M, I, J, K, J1, K1, L, L1,Z,P,O,ll
  INTEGER		IT, JT, KT, LT, SWITCH,ITA, ZT,PT,OT,tmp
  INTEGER	::	H,Q=0,X=0,NOWRITE=0,COUNTER=0
  INTEGER,DIMENSION(500):: ATOM_TORS = 0, ATOM_ANG = 0
  INTEGER, DIMENSION(10,10,10,10) :: NOWRITE2 = 0
  !       *******************************************************************

  !       SWITCH = 0
  IF (NBTYPE.GT.0) THEN
     I = 0

     NIJK = 0
     NIJKL = 0
     FNIJKL = 0
     NOOPIJKL = 0
     NOANGLEIJK = 0

     DO 140 N = 1, NMOL
        DO 130 M = 1, NATM(N)
           I = I + 1
           IT = ITYPE(I)

           !**** Search for possible improper torsion ****

           IF (NOTYPE.GT.0) THEN
              IF (NBONDS(I).EQ.3) THEN

                do ll = 1,6
                    if(ll .eq. 1)then ! questa era quella che c'era
                        Z = JBOND(I,1)
                        P = JBOND(I,2)
                        O = JBOND(I,3)
                    elseif(ll .eq. 2)then
                        Z = JBOND(I,3)
                        P = JBOND(I,1)
                        O = JBOND(I,2)
                    elseif(ll .eq. 3)then
                        Z = JBOND(I,2)
                        P = JBOND(I,1)
                        O = JBOND(I,3)
                    elseif(ll .eq. 4)then ! questa era quella che c'era con (ll .eq. 2)
                        Z = JBOND(I,1)
                        P = JBOND(I,3)
                        O = JBOND(I,2)
                   elseif(ll .eq. 5)then
                        Z = JBOND(I,2)
                        P = JBOND(I,3)
                        O = JBOND(I,1)
                   elseif(ll .eq. 6)then
                        Z = JBOND(I,3)
                        P = JBOND(I,2)
                        O = JBOND(I,1)
                   end if
                   ZT = ITYPE(Z)
                   PT = ITYPE(P)
                   OT = ITYPE(O)

                    IF (IOOPT(IT,ZT,PT,OT).NE.0 .and. NOOPIJKL(I) .eq. 0) THEN
                       !		        WRITE(*,*) 'found',i,z,p,o,ll
                       NOWRITE2(IT,ZT,PT,OT) = NOWRITE2(IT,ZT,PT,OT)+ 1
                       IF (NOWRITE2(IT,ZT,PT,OT) .EQ. 1) THEN
                          WRITE(*,*) '             ********************* WARNING *********************'
                          WRITE(*,*) '             *     You defined an improper harmonic torsion    *'
                          WRITE(*,*) '             *               for the type atoms                *'
                          WRITE(*,700) IT,ZT,PT,OT  
700                       FORMAT('              *               ',4(I4,2X),'          *')
                          WRITE(*,*) '             ***************************************************'
                          WRITE(*,*)		
                          WRITE(1,*)
                          WRITE(1,*) '             ********************* WARNING *********************'
                          WRITE(1,*) '             *     You defined an improper harmonic torsion    *'
                          WRITE(1,*) '             *               for the type atoms                *'
                          WRITE(1,700) IT,ZT,PT,OT  
                          WRITE(1,*) '             ***************************************************'
                          WRITE(1,*)
                       END IF

                       !**** We have a valid OOP ****

                       NOOPIJKL(I) = NOOPIJKL(I) + 1
                       JOOPIJKL(I,NOOPIJKL(I)) = Z
                       KOOPIJKL(I,NOOPIJKL(I)) = P
                       LOOPIJKL(I,NOOPIJKL(I)) = O

                    END IF
                 end do

              END IF !IF NBONDS == 3
           END IF !If NOTYPE > 1

           DO 120 J1 = 1,NBONDS(I)
              J = JBOND(I,J1)
              JT = ITYPE(J)

              !**** Check for bonds ****

              IF (IBONDT(IT,JT).EQ.0) THEN
                 !**** No parameters defined for this bond type! ****
                 WRITE (1,*) ' '
                 WRITE (1,*) &
                      ' **** FATAL ERROR! FATAL ERROR! ****'
                 WRITE (1,*) ' No parameters found for the bond'
                 WRITE (1,*) ' between  atoms',I,' & ',J
                 WRITE (1,*) ' of atom types ',IT,' & ',JT
                 WRITE (1,*) ' **** Check interaction file ****'
                 WRITE (*,*) ' '
                 WRITE (*,*) &
                      ' **** FATAL ERROR! FATAL ERROR! ****'
                 WRITE (*,*) ' No parameters found for the bond'
                 WRITE (*,*) ' between  atoms',I,' & ',J
                 WRITE (*,*) ' of atom types ',IT,' & ',JT
                 WRITE (*,*) ' **** Check interaction file ****'
                 ISTOP = 1
                 RETURN

              END IF

              !**** Search for possible bond angles? ****

              IF (NATYPE.GT.0) THEN

                 DO 110 K1 = 1,NBONDS(J)
                    K = JBOND(J,K1)
                    KT = ITYPE(K)
                    SWITCH = 0

                    DO ITA = 1,NATYPE*3-2,3
                       IF (ATIANG(ITA).EQ.IT) THEN
                          IF(ATIANG(ITA+1).EQ.JT) THEN
                             IF(ATIANG(ITA+2).EQ.KT) THEN 
                                SWITCH = 1
                             ENDIF
                          ENDIF
                       ENDIF
                    ENDDO
                    DO ITA = 3,NATYPE*3,3
                       IF (ATIANG(ITA).EQ.IT) THEN
                          IF(ATIANG(ITA-1).EQ.JT) THEN
                             IF(ATIANG(ITA-2).EQ.KT) THEN
                                SWITCH = 1
                             ENDIF
                          ENDIF
                       ENDIF
                    ENDDO
                    IF (K.EQ.I) GO TO 110
                    !                      IF (K.GT.I.and. SWITCH.EQ.1) THEN
                    IF (K.GT.I) THEN


                       IF (IANGT(IT,JT,KT).EQ.0) THEN

                          NOWRITE = 0

                          DO H = 1,X,3
                             IF ((IT+JT+KT) .EQ. (ATOM_ANG(H)+ATOM_ANG(H+1)+ATOM_ANG(H+2))) THEN
				IF (IT .EQ. ATOM_ANG(H)) THEN
                                   IF (JT .EQ. ATOM_ANG(H+1)) THEN
                                      IF (KT .EQ. ATOM_ANG(H+2)) THEN
                                         NOWRITE = 1
                                         EXIT
                                      END IF
                                   END IF
				END IF
                             END IF
                             IF ((IT+JT+KT) .EQ. (ATOM_ANG(H+1)+ATOM_ANG(H+2)+ATOM_ANG(H))) THEN
				IF (IT .EQ. ATOM_ANG(H+2)) THEN
                                   IF (JT .EQ. ATOM_ANG(H+1)) THEN
                                      IF (KT .EQ. ATOM_ANG(H)) THEN
                                         NOWRITE = 1
                                         EXIT
                                      END IF
                                   END IF
				END IF
                             END IF
                          END DO

                          IF (NOWRITE .EQ. 0) THEN

                             X = X+3
                             ATOM_ANG(X-2) = IT
                             ATOM_ANG(X-1) = JT
                             ATOM_ANG(X)   = KT

                             !			WRITE (*,*) &
                             !                        ' **** WARNING ****'
                             !                      WRITE (*,*) ' No parameters found for the angle'
                             !                     WRITE (*,570) IT,JT,KT
                             !				570 FORMAT(1X, 'of atom types:',3(3X,I3))

                          END IF


                          WRITE (1,*) ' No parameters found for the angle'
                          WRITE (1,*) ' between  atoms',I,' & ',J,' & ',K
                          WRITE (1,*) ' of atom types ',IT,' & ',JT,' & ',KT
                          WRITE (1,*) ' **** Check interaction file ****'

                          NOANGLE = 1

                          !**** No parameters defined for this angle type! ****
                          !                        WRITE (1,*) ' '
                          !                        WRITE (1,*) &
                          !                         ' **** FATAL ERROR! FATAL ERROR! ****'
                          !                        WRITE (1,*) ' No parameters found for the angle'
                          !                        WRITE (1,*) ' between  atoms',I,' & ',J,' & ',K
                          !                        WRITE (1,*) ' of atom types ',IT,' & ',JT,' & ',KT
                          !                        WRITE (1,*) ' **** Check interaction file ****'
                          !                        WRITE (*,*) ' '
                          !                        WRITE (*,*) &
                          !                         ' **** FATAL ERROR! FATAL ERROR! ****'
                          !                        WRITE (*,*) ' No parameters found for the angle'
                          !                        WRITE (*,*) ' between  atoms',I,' & ',J,' & ',K
                          !                        WRITE (*,*) ' of atom types ',IT,' & ',JT,' & ',KT
                          !                        WRITE (*,*) ' **** Check interaction file ****'
                          !                        ISTOP = 1
                          !                        RETURN
                          !	                END IF

                          NOANGLEIJK(I) = NOANGLEIJK(I) + 1
                          NOJANGLEIJK(I,NOANGLEIJK(I)) = J
                          NOKANGLEIJK(I,NOANGLEIJK(I)) = K

                          !**** We have found a valence angle so add it to the list ****

                       ELSE

                          NIJK(I) = NIJK(I) + 1
                          JANGLEIJK(I,NIJK(I)) = J
                          KANGLEIJK(I,NIJK(I)) = K

                       END IF
                    END IF
                    !**** Search for possible torsions? ****

                    IF (NTTYPE.GT.0) THEN
                       DO 100 L1 = 1,NBONDS(K)
                          L = JBOND(K,L1)
                          LT = ITYPE(L)

                          IF (L.EQ.J) GO TO 100

                          IF (L.GT.I) THEN
                             check_TORS:IF (ITORT(IT,JT,KT,LT).EQ.0) THEN
                                !**** No parameters defined for this torsion type! ****	

                                check_YES_NOTYPE:IF (NOTYPE.GT.0) THEN
                                   checkOOP:IF(IOOPT(IT,JT,KT,LT) .EQ. 0) THEN

                                      NOWRITE = 0

                                      DO H = 1,Q,4
                                         IF ((IT+JT+KT+LT) .EQ. (ATOM_TORS(H+1)+ATOM_TORS(H+2)+ATOM_TORS(H+3)+ATOM_TORS(H))) THEN
                                            IF (IT .EQ. ATOM_TORS(H)) THEN
                                               IF (JT .EQ. ATOM_TORS(H+1)) THEN
                                                  IF (KT .EQ. ATOM_TORS(H+2)) THEN
                                                     IF (LT .EQ. ATOM_TORS(H+3)) THEN
                                                        NOWRITE = 1
                                                        EXIT
                                                     END IF
                                                  END IF
                                               END IF
                                            END IF
                                         END IF
                                         IF ((IT+JT+KT+LT) .EQ. (ATOM_TORS(H+1)+ATOM_TORS(H+2)+ATOM_TORS(H+3)+ATOM_TORS(H))) THEN
                                            IF (IT .EQ. ATOM_TORS(H+3)) THEN
                                               IF (JT .EQ. ATOM_TORS(H+2)) THEN
                                                  IF (KT .EQ. ATOM_TORS(H+1)) THEN
                                                     IF (LT .EQ. ATOM_TORS(H)) THEN
                                                        NOWRITE = 1
                                                        EXIT
                                                     END IF
                                                  END IF
                                               END IF
                                            END IF
                                         END IF
                                      END DO

                                      IF (NOWRITE .EQ. 0) THEN

                                         Q = Q+4
                                         ATOM_TORS(Q-3) = IT
                                         ATOM_TORS(Q-2) = JT
                                         ATOM_TORS(Q-1) = KT
                                         ATOM_TORS(Q)   = LT				

                                         !				WRITE (*,*) &
                                         !               	         ' **** WARNING ****'
                                         !              	        WRITE (*,*) ' No parameters found for the torsion'
                                         !               	         WRITE (*,500) I,J,K,L
                                         !					500 FORMAT(1X, 'between atoms:',4(3X,I3))
                                         !                	        WRITE (*,510) IT,JT,KT,LT
                                         !					510 FORMAT(1X, 'of atom types:',4(3X,I3))
                                         !              	         WRITE (*,*) ' **** Check interaction file ****'
                                         !				WRITE (*,*)

                                      END IF

                                      WRITE (1,*) ' '
                                      WRITE (1,*) &
                                           ' **** FATAL ERROR! FATAL ERROR! ****'
                                      WRITE (1,*) ' No parameters found for the torsion'
                                      WRITE (1,*) ' between  atoms',I,' & ',J,' & ',K,' & ',L
                                      WRITE (1,*) ' of atom types ',IT,' & ',JT,' & ',KT,' & ',LT
                                      WRITE (1,*) ' **** Check interaction file ****'
                                      !                        WRITE (*,*) ' '

                                      !                        ISTOP = 1
                                      !                        RETURN

                                      NOTORS = 1

                                      FNIJKL(I) = FNIJKL(I) + 1
                                      FJTORIJKL(I,FNIJKL(I)) = J
                                      FKTORIJKL(I,FNIJKL(I)) = K
                                      FLTORIJKL(I,FNIJKL(I)) = L

                                   ELSE

                                      NOOPIJKL(I) = NOOPIJKL(I) + 1
                                      JOOPIJKL(I,NOOPIJKL(I)) = J
                                      KOOPIJKL(I,NOOPIJKL(I)) = K
                                      LOOPIJKL(I,NOOPIJKL(I)) = L

                                   END IF checkOOP

                                ELSE

                                   NOWRITE = 0

                                   DO H = 1,Q,4
                                      IF ((IT+JT+KT+LT) .EQ. (ATOM_TORS(H+1)+ATOM_TORS(H+2)+ATOM_TORS(H+3)+ATOM_TORS(H))) THEN
                                         IF (IT .EQ. ATOM_TORS(H)) THEN
                                            IF (JT .EQ. ATOM_TORS(H+1)) THEN
                                               IF (KT .EQ. ATOM_TORS(H+2)) THEN
                                                  IF (LT .EQ. ATOM_TORS(H+3)) THEN
                                                     NOWRITE = 1
                                                     EXIT
                                                  END IF
                                               END IF
                                            END IF
                                         END IF
                                      END IF
                                   END DO

                                   IF (NOWRITE .EQ. 0) THEN

                                      Q = Q+4

                                      ATOM_TORS(Q-3) = IT
                                      ATOM_TORS(Q-2) = JT
                                      ATOM_TORS(Q-1) = KT
                                      ATOM_TORS(Q)   = LT

                                      !			WRITE (*,*) &
                                      !                        ' **** WARNING ****'
                                      !                      WRITE (*,*) ' No parameters found for the torsion'
                                      !                        WRITE (*,500) I,J,K,L
                                      !				500 FORMAT(1X, 'between atoms:',4(3X,I3))
                                      !                     WRITE (*,550) IT,JT,KT,LT
                                      !			550 FORMAT(1X, 'of atom types:',4(3X,I3))
                                      !                       WRITE (*,*) ' **** Check interaction file ****'
                                      !		WRITE (*,*)

                                   END IF

                                   WRITE (1,*) ' '
                                   WRITE (1,*) &
                                        ' **** FATAL ERROR! FATAL ERROR! ****'
                                   WRITE (1,*) ' No parameters found for the torsion'
                                   WRITE (1,*) ' between  atoms',I,' & ',J,' & ',K,' & ',L
                                   WRITE (1,*) ' of atom types ',IT,' & ',JT,' & ',KT,' & ',LT
                                   WRITE (1,*) ' **** Check interaction file ****'
                                   !                        WRITE (*,*) ' '

                                   !                        ISTOP = 1
                                   !                        RETURN

                                   NOTORS = 1

                                   FNIJKL(I) = FNIJKL(I) + 1
                                   FJTORIJKL(I,FNIJKL(I)) = J
                                   FKTORIJKL(I,FNIJKL(I)) = K
                                   FLTORIJKL(I,FNIJKL(I)) = L


                                END IF check_YES_NOTYPE

                             ELSE

                                !**** We have found a torsion angle so add it to the list ****		
                                NIJKL(I) = NIJKL(I) + 1
                                JTORIJKL(I,NIJKL(I)) = J
                                KTORIJKL(I,NIJKL(I)) = K
                                LTORIJKL(I,NIJKL(I)) = L


                             END IF check_TORS

                          END IF

100                       CONTINUE !Loop over bonds on K (Torsion search)
                       END IF

110                    CONTINUE !Loop over all angles
                    END IF

120                 CONTINUE !Loop over all bonds

130                 CONTINUE !Loop over all atoms within molecules
140                 CONTINUE !Loop over all molecules
                 END IF

                 IF (NOANGLE .EQ. 1) THEN

                    WRITE(*,*) '             ********************* WARNING *********************'
                    WRITE(*,*) '             *       No parameters found for the angle(s)      *'
                    DO H=3,X,3
                       WRITE (*,570) ATOM_ANG(H-2),ATOM_ANG(H-1),ATOM_ANG(H)
570                    FORMAT(2X, '            *         of atom types:',3(3X,I3),'        *')
                    END DO
                    WRITE(*,*) '             ***************************************************'
                    WRITE(*,*)

                 END IF

                 IF (NOTORS .EQ. 1) THEN

                    WRITE(*,*) '             ********************* WARNING *********************'
                    WRITE(*,*) '             *       No parameters found for the torsion(s)    *'
                    DO H=4,Q,4
                       WRITE (*,510) ATOM_TORS(H-3),ATOM_TORS(H-2),ATOM_TORS(H-1),ATOM_TORS(H)
510                    FORMAT(2X, '            *        of atom types:',4(3X,I3),'   *')
                    END DO
                    WRITE(*,*) '             ***************************************************'
                    WRITE(*,*)

                 END IF


                 ! Find the highest connectivity

                 if(ibrdescr .eq. 0)then
                    contactA = 0
                    contactB = 0
                    tmp=0
                    do i=1,natoms-1 
                       do j=i+1,natoms
                          if(type_label(i) .eq. 1 .and. type_label(j) .eq. 1)then
                             call non_bond_array(i,j)
                             if(nonbond .ne. 1)then
                                tmp = abs(j-i)+1
                                if(contactA .lt. tmp)contactA=tmp
                             end if
                          elseif(type_label(i) .eq. 2 .and. type_label(j) .eq. 2)then
                             call non_bond_array_bead(i,j)
                             if(nonbond .ne. 1)then
                                tmp = abs(j-i)+1
                                if(contactB .lt. tmp)contactB=tmp
                             end if
                          end if
                       end do
                    end do
                 else
                    contactA = 0
                    tmp=0
                    do i=1,natoms 
                       do j=i+1,natoms
                          call non_bond_array(i,j)
                          if(nonbond .ne. 1)then
                             tmp = abs(j-i)
                             if(contactA .lt. tmp)contactA=tmp+1
                          end if
                       end do
                    end do
                 end if


                 RETURN
               END SUBROUTINE SETLIS

               !	***************************************************************************

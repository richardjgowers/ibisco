!> @file
!> @brief Neighbourlist creation with linked cell lists
!> @details Called from UPDATE_NEIGHBOURLIST.f90

SUBROUTINE NEW_NEIGHBOUR_WITHLIST(INDEX_LIST,CELL,LCLIST,N,RLIST,LIST,MAXNAB &
     , MAP, SIZEMAP, HEAD, MAXNUMCELL)

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: N !! @param n Number of particles to process
  INTEGER, INTENT(IN) :: MAXNUMCELL, MAXNAB
  INTEGER, INTENT(IN) :: INDEX_LIST(N) !! Array which translates onto the master coordinate list
  INTEGER, INTENT(IN) :: CELL(N), LCLIST(N), SIZEMAP, MAP(SIZEMAP), HEAD(MAXNUMCELL)
  INTEGER, INTENT(INOUT) :: LIST(MAXNAB,N)
  INTEGER :: A, B, C, D, I, J, NLIST, JCELL, JCELL0, NABOR
  INTEGER :: NONBOND
  REAL*4, INTENT(IN) :: RLIST
  REAL*4, DIMENSION(3) :: RXYZ_I, RXYZ_IJ
  REAL*4 :: RIJSQ
  REAL*4 :: RLISTSQ


  RLISTSQ = RLIST * RLIST
  LIST = 0 !Should set all unoccupied list slots to 0, can then use this to detect end of list
  !***************************BUILD UP THE NEIGHBOUR LIST*********************
  !***********LOOP IN THE CELL WHERE ATOM I IS **************************
  !$OMP PARALLEL DO DEFAULT(NONE), SCHEDULE(STATIC,1)&
  !$OMP& SHARED(N,INDEX_LIST,CELL,RXYZ,LCLIST,CONNECTIONS,CONNECTED_TO,IBRDESCR,NATOMS)&
  !$OMP& SHARED(BOX,BOXINV,RLISTSQ,LIST,MAXNAB,MAP,HEAD, NNEBS)&
  !$OMP& PRIVATE(A,I,NLIST,JCELL,RXYZ_I,B,J,NONBOND,C,D,RXYZ_IJ,RIJSQ,JCELL0,NABOR)
  DO A = 1,N
     I = INDEX_LIST(A)
     NLIST = 0 !Counter for number of neighbours this atom has

     JCELL = CELL(A)

     RXYZ_I(:) = RXYZ(:,I)
  
     B = LCLIST(A)
     DO
        IF(B .eq. 0) EXIT !If reached end of list
        J = INDEX_LIST(B) !Fetch real index of B

        NONBOND=1
        IF(CONNECTIONS(I) .gt. 0) THEN !If I has connections:
           DO C=1,CONNECTIONS(I) !Connections is the number of connected atoms to I
              D = CONNECTED_TO(I,C) !Connected_to has the index of all connected atoms
              IF(J .eq. D) THEN !If J is D then they are connected
                 NONBOND=0 
                 EXIT
              END IF
           END DO
        END IF
        IF(IBRDESCR .eq. 0) THEN
           IF(I .gt. NATOMS .and. J .gt. NATOMS) THEN !Exclude VS-VS interactions
              NONBOND = 0
           END IF
        END IF

        IF (NONBOND.EQ.1) THEN
           RXYZ_IJ(:) = RXYZ_I(:) - RXYZ(:,J)

           RXYZ_IJ(:) = RXYZ_IJ(:) - ANINT(RXYZ_IJ(:) * BOXINV(:)) * BOX(:)
           RIJSQ = SUM(RXYZ_IJ(:) * RXYZ_IJ(:))

           IF (RIJSQ.LT.RLISTSQ) THEN
              NLIST = NLIST + 1
              LIST(NLIST,A) = J
              IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
           ENDIF
        ENDIF
        B = LCLIST(B)
     END DO

        !***********LOOP IN THE NEIGHBOURING CELLS**********************************
     JCELL0 = 13*(JCELL - 1)
     DO NABOR = 1,13
        JCELL = MAP(JCELL0+NABOR)
        B = HEAD(JCELL)

        DO
           IF(B .eq. 0) EXIT
           J = INDEX_LIST(B)

           NONBOND=1
           IF(CONNECTIONS(I) .gt. 0) THEN !If I has connections:
              DO C=1,CONNECTIONS(I)
                 D = CONNECTED_TO(I,C)
                 IF(J .eq. D) THEN !If J is C then they are connected
                    NONBOND=0 
                    EXIT
                 END IF
              END DO
           END IF
           IF(IBRDESCR .eq. 0) THEN
              IF(I .gt. NATOMS .and. J .gt. NATOMS) THEN
                 NONBOND = 0
              END IF
           END IF

           IF (NONBOND.EQ.1) THEN
              RXYZ_IJ(:) = RXYZ_I(:) - RXYZ(:,J)

              RXYZ_IJ(:) = RXYZ_IJ(:) - ANINT(RXYZ_IJ(:) * BOXINV(:)) * BOX(:)
              RIJSQ = SUM(RXYZ_IJ(:) * RXYZ_IJ(:))

              IF (RIJSQ.LT.RLISTSQ) THEN
                 NLIST = NLIST + 1
                 LIST(NLIST,A) = J
                 IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
              ENDIF

           ENDIF
           B = LCLIST(B)
        END DO
     END DO

     NNEBS(I) = NLIST ! Save number of neighbours this particle has

  END DO
  !$OMP END PARALLEL DO

  RETURN
END SUBROUTINE NEW_NEIGHBOUR_WITHLIST


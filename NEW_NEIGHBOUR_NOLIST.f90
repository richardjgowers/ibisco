!> @file
!> @brief Neighbour list creation without linked cell lists, generally the slowest option
!> @note Not currently used, kept in for debugging purposes
!> @todo Parallelise this if it's ever used
SUBROUTINE NEW_NEIGHBOUR_NOLIST(N, INDEX_LIST, MAXNAB, LIST, RLIST)

  USE VAR

  INTEGER, INTENT(IN) :: N !< Number of particles in index list
  INTEGER, INTENT(IN) :: INDEX_LIST(N) !< Gives global index of each particle
  INTEGER, INTENT(IN) :: MAXNAB !< Maximum possible number neighbours per particle
  INTEGER, INTENT(INOUT) :: LIST(MAXNAB, N) !< The neighbour list
  REAL*4, INTENT(IN) :: RLIST !< Cutoff radius for the neighbour list
  INTEGER :: A, B, C, D, I, J, NLIST
  REAL*4, DIMENSION(3) :: RXYZ_I
  REAL*4 :: RIJSQ
  REAL*4, DIMENSION(3) :: RXYZ_IJ
  REAL*4 :: RLISTSQ

  RLISTSQ = RLIST * RLIST
  LIST = 0 ! Set neighbour list to 0, this indicates end of list

  DO A=1,N-1
     I = INDEX_LIST(A)
     NLIST = 0

     RXYZ_I(:) = RXYZ(:,I)

     DO B=A+1,N
        J = INDEX_LIST(B)
        !Check if connected
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
        !Check if VS-VS
        IF(IBRDESCR .eq. 0) THEN
           IF(I .gt. NATOMS .and. J .gt. NATOMS) THEN
              NONBOND = 0
           END IF
        END IF

        IF(NONBOND .eq. 1) THEN
           RXYZ_IJ(:)  = RXYZ_I(:) - RXYZ(:,J)

           RXYZ_IJ(:) = RXYZ_IJ(:) - ANINT(RXYZ_IJ(:) * BOXINV(:)) * BOX(:)

           RIJSQ = SUM(RXYZ_IJ(:) * RXYZ_IJ(:))

           IF (RIJSQ.LT.RLISTSQ) THEN
              NLIST = NLIST +1
              LIST(NLIST, A) = J              
              IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
           END IF
        END IF
     END DO
  END DO

RETURN 
END SUBROUTINE NEW_NEIGHBOUR_NOLIST

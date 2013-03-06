!
!	Written by Nicodemo
!
!	Dec 2010
!
! 	Structure of the subroutine:
!
!
!		DO (a loop among all the particles of the system)
!
!			IF(the type of the particle is 'BEAD')	
!				In this case it searchs among the bead-bead and among the bead-virtual site using the 
!				cut off distance specified for the bead (that in general is higher than that for the atoms
!			ELSE(the type of the particle is 'ATOM')
!				In this case the search is made only among the pairs atom-atom
!			END IF
!
!	      END DO
!
!

SUBROUTINE VIRT_NEIGHBOR_WITHLIST_COM ()

USE VAR
IMPLICIT NONE

integer :: NLIST, I, J, K, JCELL0, JCELL, NABOR,IP,JP,NLISTDPD,VNLIST
integer :: VNLISTDPD, VNLISTSEC
real(kind=rkind) :: RXI, RYI, RZI
real(kind=rkind) :: RIJSQ
real(kind=rkind) :: RXIJ, RYIJ, RZIJ


!       *******************************************************************

!		SETS UP A LIST OF 13 NEIGHBOURING
!		CELLS OF EACH OF THE SMALL CELLS IN THE CENTRAL BOX.

IF (ENSEMBLE == 2) THEN
    NCELLX = BOXX / RLIST
    NCELLY = BOXY / RLIST
    NCELLZ = BOXZ / RLIST

    NUMCELL = NCELLX * NCELLY * NCELLZ 
    CALL MAPS ()
ENDIF

 CALL LINKS ()
 CALL VIRT_LINKS_COM ()
            
    NLIST = 0
    VNLIST = 0  

    DO IP = 1,NATOMS

        POINT(IP) = NLIST + 1
        VIRT_POINT(IP) = VNLIST + 1
        JCELL = CELL(IP)

!IF IP is a Bead the list for the virtual atom is written 

        IF(TYPE_LABEL(IP) .EQ. 2) THEN
            RXI = RX(IP)
            RYI = RY(IP)
            RZI = RZ(IP)          
            JP = LCLIST(IP)

 600       IF (JP.GT.0) THEN

        if(type_label(jp) .eq. 2)then ! if jp is not a Bead it has to be discarded form Nlist 
            if(abs(ip-jp) .le. contactB)then
               CALL NON_BOND_ARRAY_BEAD(IP,JP)
            else
               NONBOND=1
            end if
        else
            NONBOND=0
        end if

           IF (NONBOND.EQ.1) THEN
		RXIJ = RXI - RX(JP)
		RYIJ = RYI - RY(JP)
		RZIJ = RZI - RZ(JP)
           
		RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
		RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
		RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
		RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2
		IF (RIJSQ.LT.RLISTIBRSQ) THEN
			NLIST = NLIST + 1
			LIST(NLIST) = JP
			IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
		ENDIF

           ENDIF

           JP = LCLIST(JP)

           GOTO 600

           ENDIF          ! 600

!*********NEIGHBOUR LIST, VIRTUAL ATOM, SAME CELL


	JCELL = CELL(IP)
    JP = VHEAD(JCELL)

	  IF (JP .GT. 0) THEN
!		JP = VLCLIST(JP)
	  	DO
		IF (JP.EQ.0) EXIT

			RXIJ = RXI - VIRTRX(JP)
			RYIJ = RYI - VIRTRY(JP)
			RZIJ = RZI - VIRTRZ(JP)
      	     
			RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
			RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
			RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ

			RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2

			IF (RIJSQ.LT.RLISTIBRSQ) THEN

				VNLIST = VNLIST + 1
				VLIST(VNLIST) = JP

				IF (VNLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'

			END IF

 !           end if

			JP = VLCLIST(JP)
	  	END DO
	  END IF

!*********NEIGHBOUR LIST, VIRTUAL ATOM, NEIGHBOURING CELLS

	           JCELL0 = 13*(JCELL - 1)
      DO NABOR = 1,13
		JCELL = MAP(JCELL0+NABOR)
		JP = VHEAD(JCELL)

	IF (JP .GT. 0) THEN
	    DO
	        IF (JP.EQ. 0) EXIT

  	              RXIJ = RXI - VIRTRX(JP)
        	      RYIJ = RYI - VIRTRY(JP) 
        	      RZIJ = RZI - VIRTRZ(JP)

        	      RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
          	      RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
            RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
      	    RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2

			IF (RIJSQ.LT.RLISTIBRSQ) THEN

				VNLIST = VNLIST + 1
				VLIST(VNLIST) = JP

                 	      IF (VNLIST.EQ.MAXNAB) STOP 'VIRT LIST TOO SMALL'
                  END IF


!                end if

                JP = VLCLIST(JP)
            END DO
	END IF

     END DO

!*********************************************************************************
!***********	LOOP IN THE NEIGHBOURING CELLS	**********************************
!*********************************************************************************

           JCELL0 = 13*(CELL(IP) - 1)
           DO NABOR = 1,13

              JCELL = MAP(JCELL0+NABOR)
              JP = HEAD (JCELL)

 400          IF (JP.GT.0) THEN
        if(type_label(jp) .eq. 2)then ! if jp is not a Bead it has to be discarded form Nlist 
            if(abs(ip-jp) .le. contactB)then
                CALL NON_BOND_ARRAY_BEAD(IP,JP)
            else
                NONBOND=1
            end if
        else
                NONBOND=0
        end if
              IF (NONBOND.EQ.1) THEN

        	      RXIJ = RXI - RX(JP)
        	      RYIJ = RYI - RY(JP) 
        	      RZIJ = RZI - RZ(JP)

         	     RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
         	     RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
         	     RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ

          	     RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2
              IF (RIJSQ.LT.RLISTIBRSQ) THEN

                 NLIST = NLIST + 1
                 LIST(NLIST) = JP
                 IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'

              ENDIF

              ENDIF
              JP = LCLIST(JP)
              GOTO 400
              ENDIF  !400
           END DO

ELSE

           RXI = RX(IP)
           RYI = RY(IP)
           RZI = RZ(IP)
           
           JP = LCLIST(IP)

 200       IF (JP.GT.0) THEN

        if(type_label(jp) .eq. 1)then ! if jp is not an Atom it has to be discarded form Nlist 
            if(abs(ip-jp) .le. contactA)then
                CALL NON_BOND_ARRAY(IP,JP)
            else
                NONBOND=1
            end if
        else
                NONBOND=0
        end if
                !CALL NON_BOND_ARRAY(IP,JP)
           IF (NONBOND.EQ.1) THEN
		RXIJ = RXI - RX(JP)
		RYIJ = RYI - RY(JP)
		RZIJ = RZI - RZ(JP)
           
		RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
		RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
		RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
		RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2

		IF (RIJSQ.LT.RLISTSQ) THEN
			NLIST = NLIST + 1
			LIST(NLIST) = JP
			IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
		ENDIF

           ENDIF

           JP = LCLIST(JP)

           GOTO 200

           ENDIF          ! 200

!*********************************************************************************
!***********	LOOP IN THE NEIGHBOURING CELLS	**********************************
!*********************************************************************************

           JCELL0 = 13*(CELL(IP) - 1)
           DO NABOR = 1,13

              JCELL = MAP(JCELL0+NABOR)
              JP = HEAD (JCELL)

 300          IF (JP.GT.0) THEN

        if(type_label(jp) .eq. 1)then ! if jp is not an Atom it has to be discarded form Nlist 
            if(abs(ip-jp) .le. contactA)then
                CALL NON_BOND_ARRAY(IP,JP)
            else
                NONBOND=1
            end if
        else
                NONBOND=0
        end if
              !CALL NON_BOND_ARRAY(IP,JP)

              IF (NONBOND.EQ.1) THEN

        	      RXIJ = RXI - RX(JP)
        	      RYIJ = RYI - RY(JP) 
        	      RZIJ = RZI - RZ(JP)

         	     RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
         	     RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
         	     RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ

          	     RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2

              IF (RIJSQ.LT.RLISTSQ) THEN

                 NLIST = NLIST + 1
                 LIST(NLIST) = JP

                 IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'

              ENDIF

              ENDIF
              JP = LCLIST(JP)
              GOTO 300
              ENDIF  !300
           END DO

END IF
        END DO
        POINT(NATOMS+1) = NLIST + 1
	VIRT_POINT(NATOMS+1) = VNLIST + 1

!     The previous part of the subroutine does not found all the virtual site near a particular bead, 
!     so the virtual site neighbour list must be updated with this loop in which the other virtual site
!     are found

!     *********************************************************************************
!     ***********	LOOP IN THE NEIGHBOURING CELLS	****************************
!     *********************************************************************************

    VNLISTSEC = 0

    DO IP = 1,NVIRTA
    
        VIRT_POINT_SEC(IP) = VNLISTSEC + 1 
    
        RXI = VIRTRX(IP)
        RYI = VIRTRY(IP)
        RZI = VIRTRZ(IP)
    
        JCELL = VCELL(IP)
        JCELL0 = 13*(VCELL(IP) - 1)
    
        DO NABOR = 1,13
    
            JCELL = MAP(JCELL0+NABOR)
            JP = HEAD (JCELL)
    
            DO
                IF (JP.EQ.0) EXIT
                IF(TYPE_LABEL(JP) .EQ. 2) THEN
                    IF (JP.EQ. 0) EXIT
                    RXIJ = RXI - RX(JP)
                    RYIJ = RYI - RY(JP) 
                    RZIJ = RZI - RZ(JP)
                    RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
                    RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
                    RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
                    RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2
    
                    IF (RIJSQ.LT.RLISTIBRSQ) THEN
                        VNLISTSEC = VNLISTSEC + 1
                        VLIST_SEC(VNLISTSEC) = JP
                        IF (VNLISTSEC.EQ.MAXNAB) STOP 'VIRT LIST TOO SMALL'
                    ENDIF
                ENDIF
                JP = LCLIST(JP)
            END DO
        END DO
    END DO
    
    VIRT_POINT_SEC(NVIRTA+1) = VNLISTSEC + 1 


      RETURN
      END
 

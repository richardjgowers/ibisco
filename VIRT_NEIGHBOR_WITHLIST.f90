!
!	Written by Nicodemo
!
!	June 2012
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

SUBROUTINE VIRT_NEIGHBOR_WITHLIST()

  USE VAR
  IMPLICIT NONE

integer :: NLIST, I, J, K, JCELL0, JCELL, NABOR,IP,JP,VNLIST
integer :: VNLISTSEC
real(kind=rkind) :: RXI, RYI, RZI
real(kind=rkind) :: RIJSQ
real(kind=rkind) :: RXIJ, RYIJ, RZIJ

!       *******************************************************************
!

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
              
NLIST = 0
VNLIST = 0

    DO IP = 1,NATOMS
        POINT(IP) = NLIST + 1
        JCELL = CELL(IP)

!IF IP is a Bead the list for the virtual atom is written 
        IF(TYPE_LABEL(IP) .EQ. 2) THEN
           RXI = RX(IP)
           RYI = RY(IP)
           RZI = RZ(IP)
           JP = LCLIST(IP)

 600        IF (JP.GT.0) THEN
                do
                if(type_label(jp) .eq. 2)then ! If jp is a Bead it checks if within the cut-off and possibly add it in the list
                    if(abs(ip-jp) .le. contactB)then
                        CALL NON_BOND_ARRAY_BEAD(IP,JP)
                    else
                        NONBOND=1
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
                else    ! If jp is an Atoms it checks if it is a VS and then if is within the cut-off and possibly add it in the Vlist
                     if(VIRT_VS_IND(jp) .ne. 0)then
		                RXIJ = RXI - RX(JP)
		                RYIJ = RYI - RY(JP)
		                RZIJ = RZI - RZ(JP)
                		RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
                		RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
		                RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
		                RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2
	      		        IF (RIJSQ.LT.RLISTIBRSQ) THEN
	      			        NLIST = NLIST + 1
	      			        LIST(NLIST) = jp
	      			        IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
	      		        END IF
                    end if
                end if
                JP = LCLIST(JP)   
                if(JP .eq. 0) exit         
                end do
           end if

!*********************************************************************************
!***********	LOOP IN THE NEIGHBOURING CELLS	**********************************
!*********************************************************************************

        JCELL0 = 13*(CELL(IP) - 1)
        DO NABOR = 1,13
            JCELL = MAP(JCELL0+NABOR)
            JP = HEAD (JCELL)
 400        IF (JP.GT.0) THEN
                do             
                    if(type_label(jp) .eq. 2)then ! If jp is a Bead it checks if within the cut-off and possibly add it in the list
                        if(abs(ip-jp) .le. contactB)then
                            CALL NON_BOND_ARRAY_BEAD(IP,JP)
                        else
                            NONBOND=1
                        end if
!                        CALL NON_BOND_ARRAY_BEAD(IP,JP)
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
                    else
                        if(VIRT_VS_IND(jp) .ne. 0)then
		                    RXIJ = RXI - RX(JP)
		                    RYIJ = RYI - RY(JP)
		                    RZIJ = RZI - RZ(JP)
                    		RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
                    		RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
		                    RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
		                    RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2
	      	    	        IF (RIJSQ.LT.RLISTIBRSQ) THEN
	      			            NLIST = NLIST + 1
	      			            LIST(NLIST) = jp
	      	    		        IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
	      	    	        END IF
                        end if
                    end if
                    JP = LCLIST(JP)
                    if(JP .eq. 0) exit         
                end do
            end if
        END DO

    ELSE ! Neighbour list for the atoms

        RXI = RX(IP)
        RYI = RY(IP)
        RZI = RZ(IP)           
        JP = LCLIST(IP)

        if(VIRT_VS_IND(ip) .ne. 0)then

! If the atom in a VS also the neighbour list has to include beads also

 201    IF (JP .GT. 0) THEN
            do  
		        RXIJ = RXI - RX(JP)
		        RYIJ = RYI - RY(JP)
		        RZIJ = RZI - RZ(JP)
              	RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
		        RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
		        RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
		        RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2
                if(type_label(JP) .eq. 1)then
                    if(abs(ip-jp) .le. contactA)then
                        CALL NON_BOND_ARRAY(IP,JP)
                    else
                        NONBOND=1
                    end if    
!                    CALL NON_BOND_ARRAY(IP,JP)
                    IF (NONBOND.EQ.1) THEN
                        IF (RIJSQ.LT.RLISTSQ) THEN
			                NLIST = NLIST + 1
			                LIST(NLIST) = JP
			                IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
		                ENDIF
	                ENDIF
                else
                    IF (RIJSQ.LT.RLISTIBRSQ) THEN
                        NLIST = NLIST + 1
	                    LIST(NLIST) = JP
		                IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
	                ENDIF
                end if! if(type_label(JP) .eq. 1)then
                 JP = LCLIST(JP)
                 if(JP .eq. 0) exit 
             end do
           ENDIF          ! 200

!*********************************************************************************
!***********	LOOP IN THE NEIGHBOURING CELLS for VS	**************************
!*********************************************************************************

        JCELL0 = 13*(CELL(IP) - 1)
        DO NABOR = 1,13
           JCELL = MAP(JCELL0+NABOR)
           JP = HEAD (JCELL)
 301        IF (JP .GT. 0 ) THEN
                do
            	    RXIJ = RXI - RX(JP)
            	    RYIJ = RYI - RY(JP) 
            	    RZIJ = RZI - RZ(JP)
             	    RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
             	    RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
             	    RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
                    RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2
                    if(type_label(JP) .eq. 1)then
                    if(abs(ip-jp) .le. contactA)then
                        CALL NON_BOND_ARRAY(IP,JP)
                    else
                        NONBOND=1
                    end if
!                        CALL NON_BOND_ARRAY(IP,JP)
                        if(nonbond .eq. 1)then
                            IF (RIJSQ.LT.RLISTSQ) THEN
			                    NLIST = NLIST + 1
			                    LIST(NLIST) = JP
			                    IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
		                    ENDIF
                        end if
                    else
                        IF (RIJSQ.LT.RLISTIBRSQ) THEN
			                NLIST = NLIST + 1
			                LIST(NLIST) = JP
			                IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
		                ENDIF
                    end if! if(type_label(JP) .eq. 1)then
                    JP = LCLIST(JP)
                    if(JP .eq. 0) exit     
                end do
              ENDIF  !300
           END DO

        else !if(VIRT_VS_IND(ip) .ne. 0)

 200    IF (JP .GT. 0) THEN
            do
                if(type_label(jp) .eq. 1)then  
                    if(abs(ip-jp) .le. contactA)then
                        CALL NON_BOND_ARRAY(IP,JP)
                    else
                        NONBOND=1
                    end if      
!                    CALL NON_BOND_ARRAY(IP,JP)
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
                end if
                JP = LCLIST(JP)
                if(JP .eq. 0) exit 
             end do
           ENDIF          ! 200

!*********************************************************************************
!***********	LOOP IN THE NEIGHBOURING CELLS	**********************************
!*********************************************************************************

        JCELL0 = 13*(CELL(IP) - 1)
        DO NABOR = 1,13
           JCELL = MAP(JCELL0+NABOR)
           JP = HEAD (JCELL)
 300        IF (JP .GT. 0 ) THEN
                do
                    if(type_label(jp) .eq. 1)then
                    if(abs(ip-jp) .le. contactA)then
                        CALL NON_BOND_ARRAY(IP,JP)
                    else
                        NONBOND=1
                    end if    
!                        CALL NON_BOND_ARRAY(IP,JP)
                        IF (NONBOND .EQ. 1) THEN
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
                    end if
                    JP = LCLIST(JP)
                    if(JP .eq. 0) exit     
                end do
              ENDIF  !300
           END DO
    end if !if(VIRT_VS_IND(ip) .ne. 0)
END IF ! (TYPE_LABEL(IP) .EQ. 2) THEN
        END DO ! do ip=1,natoms

    POINT(NATOMS+1) = NLIST + 1

RETURN
END
 

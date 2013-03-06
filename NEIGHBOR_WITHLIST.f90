        SUBROUTINE NEIGHBOR_WITHLIST ()

	USE VAR
	IMPLICIT	NONE
	INTEGER		NLIST, I, J, K, JCELL0, JCELL, NABOR,IP,JP,NLISTDPD
	REAL*8		RXI, RYI, RZI
	REAL*8		RIJSQ
	REAL*8		RXIJ, RYIJ, RZIJ
!       *******************************************************************
	IF (ENSEMBLE == 2) THEN
!		SETS UP A LIST OF 13 NEIGHBOURING
!		CELLS OF EACH OF THE SMALL CELLS IN THE CENTRAL BOX.
		NCELLX = BOXX / RLIST
		NCELLY = BOXY / RLIST
		NCELLZ = BOXZ / RLIST

        	NUMCELL = NCELLX * NCELLY * NCELLZ 
		CALL MAPS ()
        ENDIF
                CALL LINKS ()
!***************************BUILD UP THE NEIGHBOUR LIST*********************
!***********LOOP IN THE CELL WHERE ATOM I IS **************************
        IF (DPDINPUT.EQ.1.OR.LAINPUT.EQ.1) THEN
        NLIST = 0
        NLISTDPD = 0 
        DO IP = 1,NATOMS
           POINT(IP) = NLIST + 1
           DPDPOINT(IP) = NLISTDPD + 1
           JCELL = CELL(IP)
           RXI = RX(IP)
           RYI = RY(IP)
           RZI = RZ(IP)
           
           JP = LCLIST(IP)
 210       IF (JP.GT.0) THEN
           RXIJ = RXI - RX(JP)
           RYIJ = RYI - RY(JP)
           RZIJ = RZI - RZ(JP)
           
           RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
           RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
           RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
           RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2
           IF (RIJSQ.LT.RLISTSQ) THEN
             call NON_BOND_ARRAY(IP,JP)
             IF (NONBOND.EQ.1) THEN
                NLIST = NLIST + 1
                LIST(NLIST) = JP
                IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
             ELSE
                NLISTDPD = NLISTDPD + 1
                DPDLIST(NLISTDPD) = JP
                IF (NLISTDPD.EQ.NATOMS*50) STOP 'ADDITIONAL DPD LIST (NATOMS*50) TOO SMALL'
             ENDIF 
           ENDIF
           JP = LCLIST(JP)
           GOTO 210
           ENDIF          ! 210
!***********LOOP IN THE NEIGHBOURING CELLS**********************************
           JCELL0 = 13*(JCELL - 1)
           DO NABOR = 1,13
              JCELL = MAP(JCELL0+NABOR)
              JP = HEAD (JCELL)
 310          IF (JP.GT.0) THEN
              RXIJ = RXI - RX(JP)
              RYIJ = RYI - RY(JP) 
              RZIJ = RZI - RZ(JP)

              RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
              RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
              RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
              RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2
              IF (RIJSQ.LT.RLISTSQ) THEN
                 call NON_BOND_ARRAY(IP,JP)
                 IF (NONBOND.EQ.1) THEN
                    NLIST = NLIST + 1
                    LIST(NLIST) = JP
                    IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
                 ELSE
                    NLISTDPD = NLISTDPD + 1
                    DPDLIST(NLISTDPD) = JP
                    IF (NLISTDPD.EQ.NATOMS*50) STOP 'ADDITIONAL DPD LIST (NATOMS*50) TOO SMALL'
                 ENDIF
              ENDIF
              JP = LCLIST(JP)
              GOTO 310
              ENDIF  !310
           ENDDO

        ENDDO
        POINT(NATOMS+1) = NLIST + 1
        DPDPOINT(NATOMS+1) = NLISTDPD + 1
!***************************************************************************
        ELSE                 !!!!!!!!!!DPDINPUT = 0, NO DDP THERMORSTAT****************************
        NLIST = 0            !!!!!!!!!!EXCLUDE THE BONDED PAIRS FROM INTERACTION LIST**************
        DO IP = 1,NATOMS
           POINT(IP) = NLIST + 1
           JCELL = CELL(IP)
           RXI = RX(IP)
           RYI = RY(IP)
           RZI = RZ(IP)
           
           JP = LCLIST(IP)
 200       IF (JP.GT.0) THEN
            if(abs(ip-jp) .le. contactA)then
               CALL NON_BOND_ARRAY(IP,JP)
            else
             NONBOND=1
            end if
           !call NON_BOND_ARRAY(IP,JP)
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
!***********LOOP IN THE NEIGHBOURING CELLS**********************************
           JCELL0 = 13*(JCELL - 1)
           DO NABOR = 1,13
              JCELL = MAP(JCELL0+NABOR)
              JP = HEAD (JCELL)
 300          IF (JP.GT.0) THEN
                if(abs(ip-jp) .le. contactA)then
                   CALL NON_BOND_ARRAY(IP,JP)
                else
                    NONBOND=1
                end if
!              call NON_BOND_ARRAY(IP,JP)
              IF (NONBOND.EQ.1) THEN
              RXIJ = RXI - RX(JP)
              RYIJ = RYI - RY(JP) 
              RZIJ = RZI - RZ(JP)

              RXIJ = RXIJ - ANINT(RXIJ*BOXXINV)*BOXX
              RYIJ = RYIJ - ANINT(RYIJ*BOXYINV)*BOXY
              RZIJ = RZIJ - ANINT(RZIJ*BOXZINV)*BOXZ
              RIJSQ = RXIJ**2 + RYIJ**2 + RZIJ**2

!if (ip .eq. 440) then
!write(2000,*) ip,jp
!write(2000,*) jcell, jcell0
!write(2000,*)
!end if
              IF (RIJSQ.LT.RLISTSQ) THEN
                 NLIST = NLIST + 1
                 LIST(NLIST) = JP
                 IF (NLIST.EQ.MAXNAB) STOP 'LIST TOO SMALL'
              ENDIF
              ENDIF
              JP = LCLIST(JP)
              GOTO 300
              ENDIF  !300
           ENDDO

        ENDDO
        POINT(NATOMS+1) = NLIST + 1
        ENDIF

     RETURN
     END
 

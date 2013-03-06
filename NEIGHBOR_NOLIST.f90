
        SUBROUTINE NEIGHBOR_NOLIST ()

	USE VAR
	IMPLICIT	NONE
	INTEGER		NLIST, I, J, K, JCELL0, JCELL, NABOR,NLISTDPD
	REAL*8		RXI, RYI, RZI
	REAL*8		RIJSQ 
	REAL*8		RXIJ, RYIJ, RZIJ
!       *******************************************************************

        IF (DPDINPUT.EQ.1.OR.LAINPUT.EQ.1) THEN
	NLIST = 0
        NLISTDPD = 0
        DO 200 I = 1, NATOMS - 1

              POINT(I) = NLIST + 1
              DPDPOINT(I) = NLISTDPD + 1

              RXI      = RX(I)
              RYI      = RY(I)
              RZI      = RZ(I)

              DO 88 J = I + 1, NATOMS
                 RXIJ  = RXI - RX(J)
                 RYIJ  = RYI - RY(J)
                 RZIJ  = RZI - RZ(J)

                 RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
                 RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
                 RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ
                 RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

                 IF ( RIJSQ .LT. RLISTSQ ) THEN
                    
                  call NON_BOND_ARRAY(I,J)
                  IF (NONBOND == 1) THEN
                    NLIST = NLIST + 1
                    LIST(NLIST) = J

!                ** REMOVE THIS CHECK IF MAXNAB IS APPROPRIATE **

                    IF ( NLIST .EQ. MAXNAB ) STOP 'LIST TOO SMALL'
                  ELSE
                    NLISTDPD = NLISTDPD + 1
                    DPDLIST(NLISTDPD) = J
                    IF (NLISTDPD.EQ.NATOMS*50) STOP 'ADDITIONAL DPD LIST (NATOMS*50) TOO SMALL'
                  ENDIF
                 ENDIF
88            CONTINUE
		
200     CONTINUE
	POINT(NATOMS) = NLIST + 1
        POINT(NATOMS+1) = NLIST + 1
        DPDPOINT(NATOMS) = NLISTDPD + 1
        DPDPOINT(NATOMS+1) = NLISTDPD + 1

        ELSE
        NLIST = 0
        DO 100 I = 1, NATOMS - 1

              POINT(I) = NLIST + 1

              RXI      = RX(I)
              RYI      = RY(I)
              RZI      = RZ(I)

              DO 99 J = I + 1, NATOMS
              call NON_BOND_ARRAY(I,J)
              IF (NONBOND == 1) THEN

                 RXIJ  = RXI - RX(J)
                 RYIJ  = RYI - RY(J)
                 RZIJ  = RZI - RZ(J)

                 RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
                 RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
                 RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ
                 RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ

                 IF ( RIJSQ .LT. RLISTSQ ) THEN

                    NLIST = NLIST + 1
                    LIST(NLIST) = J

!                ** REMOVE THIS CHECK IF MAXNAB IS APPROPRIATE **

                    IF ( NLIST .EQ. MAXNAB ) STOP 'LIST TOO SMALL'

                 ENDIF
              END IF
99            CONTINUE

100     CONTINUE
        POINT(NATOMS) = NLIST + 1
        POINT(NATOMS+1) = NLIST + 1

        ENDIF

        RETURN
        END

!	*********************************************************************************************

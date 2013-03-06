        SUBROUTINE BONDTABLE ( )

	USE VAR
	IMPLICIT	NONE
	INTEGER		M, I, K
	REAL(KIND=RKIND) ::		BIN, A, OMEGA, XC, GP, GF, VIJ, &
& R, TMP, CON, BINBINV
	PARAMETER (CON = 1.2533141)
!    *******************************************************************
!	MAKE TABLE FOR STRETCHING POTENTIAL AND FORCE BY BOLTZMAN INVERSION

	BIN = RCUT / NDATA
	BINB = BIN
	BINBINV  = 1.0/BIN
	R = 0

	DO K = 1, NDATA

		R = K * BIN
		DO I = 1, NBTYPE
		RBOND(K,I) = R

		M=LB(I)
        	GP = 0.0
        	GF = 0.0

        	DO WHILE (M.LT.LB(I+1))
!		CONVERT TOTAL AREA (Ai), CENTER (XCi) AND WIDTH OF GAUSSIAN 
!		FUNCTIONS TO REDUCED UNIT (nanometer). DISTRIBUTION SHOULD BE NORMALIZED

	        A = GB(M,I)/10.0
       		OMEGA = GB(M+1,I)/10.0
       		XC = GB(M+2,I)/10.0
	
        	TMP = EXP(-2.0*((R-XC)/OMEGA)**2.0)
        	TMP = (A/(OMEGA*CON))*TMP
	
        	GP = GP + TMP
        	GF = GF + TMP*((R-XC)/OMEGA**2.0)
        	M = M + 3
        	END DO

		IF (GP==0.0) GP=1.0E-10
	       	VIJ = - LOG(GP)*TIN

!		TABLE FOR FORCE WILL BE FORCE DIVIDED BY DISTANCE
!#	FIJ = - 4.0*TIN*GF/(R*GP)
!#	BOND_FORCE(K,I) = FIJ
		BOND_POT(K,I)   = VIJ
		
		END DO
	END DO

	DO I = 1, NBTYPE

		RBOND(0,I) = 0.0
!#	BOND_FORCE(0,I) = BOND_FORCE(1,I)
		BOND_POT(0,I) = BOND_POT(1,I)
	END DO

        RETURN
        END

!	*********************************************************************************************

        SUBROUTINE ANGLETABLE ( )
	USE VAR
	
	IMPLICIT	NONE
	INTEGER		M, I, K
	REAL(KIND = RKIND) ::		BIN, A, OMEGA, XC, GP, GF, VIJ,  R, TMP, CON, BINAINV
	PARAMETER (CON = 1.2533141)
!	***************************************************************
!	MAKE TABLE FOR BENDING POTENTIAL AND FORCE BY BOLTZMAN INVERSION

	BIN = 180.0 / NDATA
	BINA = BIN
	BINAINV = 1.0/BIN

	R = 0
	DO K = 1, NDATA

		R = K * BIN
		DO I = 1, NATYPE
		ANGLE(K,I) = R

		M=LA(I)
        	GP = 0.0
        	GF = 0.0
		
        	DO WHILE (M.LT.LA(I+1))
	        A = GA(M,I)
       		OMEGA = GA(M+1,I)
       		XC = GA(M+2,I)
	
        	TMP = EXP(-2.0*((R-XC)/OMEGA)**2.0)
        	TMP = (A/(OMEGA*CON))*TMP
        	GP = GP + TMP
        	GF = GF + TMP*((R-XC)/OMEGA**2.0)
	
        	M = M + 3
        	END DO
        	VIJ = - LOG(GP)*TIN

!		WE HAVE TO DIVIDE THE FORCE D2R, BECAUSE DERIVATIVE OF THE POTENTIAL 
!		SHOULD BE RESPECT TO RADIAN
!#		FIJ = - 4.0*TIN*GF/GP/D2R
!#	BEND_FORCE(K,I) = FIJ
		BEND_POT(K,I)   = VIJ
		END DO
	END DO

	DO I = 1, NATYPE

		ANGLE(0,I) = 0.0
!#		BEND_FORCE(0,I) = BEND_FORCE(1,I)
		BEND_POT(0,I) = BEND_POT(1,I)
	END DO
        RETURN
        END
!	*********************************************************************************************

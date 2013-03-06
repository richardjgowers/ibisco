



SUBROUTINE FPANGLE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)
USE VAR
IMPLICIT	NONE
        INTEGER       I, J, K, M, L 
        REAL*8        RCUTSQ, W, ALPHA
        REAL*8        RXI, RYI, RZI, FXI, FYI, FZI, CT
        REAL*8        RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
	REAL*8        RXKJ, RYKJ, RZKJ, RKJSQ, THETA
	REAL*8        RXM, RYM, RZM, RM, RXN, RYN, RZN, RJI, RKI
	REAL*8	      RXS, RYS, RZS
        REAL*8        VIJ, WIJ, FIJ, RIJ, FIJK, RKJ, FXIJK, FYIJK, FZIJK
	REAL*8	      FXIJKL, FYIJKL, FZIJKL, FIJKL, PHI, COST, SIGNT, COSM, SINM
	REAL*8	      RXJK, RYJK, RZJK, RXKL, RYKL, RZKL, RJK, RKL, RN
	REAL*8	      COTANN, COTANM, COSN, SINN
	REAL*8	      RXIK, RYIK, RZIK, RIK, VIJK
	INTEGER       TI, TJ, TIJ, NI, TK, TIJK
!      *******************************************************************

	DO 298 L = 1, NIJK(I)

		J = JANGLEIJK(I, L)
		K = KANGLEIJK(I, L)

		TJ = ITYPE(J)
		TK = ITYPE(K)
	      	TIJK = IANGT(TI, TJ, TK)

		RXIJ = RXI - SX(J)
              	RYIJ = RYI - SY(J)
              	RZIJ = RZI - SZ(J)

		RXKJ = SX(K) - SX(J)
              	RYKJ = SY(K) - SY(J)
              	RZKJ = SZ(K) - SZ(J)

		RIJ = SQRT(RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0)
		RKJ = SQRT(RXKJ ** 2.0D0 + RYKJ ** 2.0D0 + RZKJ ** 2.0D0)
		CT = (RXIJ*RXKJ + RYIJ*RYKJ + RZIJ*RZKJ)/RIJ/RKJ
		THETA = DACOS(CT) * R2D

!	****** vector M cross product -IJxKJ	
		RXM = - ( RYIJ*RZKJ - RZIJ*RYKJ)
		RYM = - (-RXIJ*RZKJ + RZIJ*RXKJ)
		RZM = - ( RXIJ*RYKJ - RYIJ*RXKJ)
		RM = SQRT(RXM ** 2.0D0 + RYM ** 2.0D0 + RZM ** 2.0D0)

		RXM = RXM / RM
		RYM = RYM / RM
		RZM = RZM / RM

!	***** vector N cross product MxIJ
		RXN = ( RYM*RZIJ - RZM*RYIJ )/RIJ
		RYN = (-RXM*RZIJ + RZM*RXIJ )/RIJ
		RZN = ( RXM*RYIJ - RYM*RXIJ )/RIJ

!	***** vector S cross product -MxKJ
		RXS = -( RYM*RZKJ - RZM*RYKJ )/RKJ
		RYS = -(-RXM*RZKJ + RZM*RXKJ )/RKJ
		RZS = -( RXM*RYKJ - RYM*RXKJ )/RKJ

		NI = INT (THETA / BINA(TIJK))
		IF (THETA == 180.0) THEN
		FIJK = BEND_FORCE(TIJK,NI)
		VIJK = BEND_POT(TIJK,NI)
		ELSE
		ALPHA=(THETA-ANGLE(TIJK,NI))/BINA(TIJK)

		FIJK = BEND_FORCE(TIJK,NI)*(1.0D0-ALPHA) &
		+ ALPHA*BEND_FORCE(TIJK,NI+1) 

		VIJK = BEND_POT(TIJK,NI)*(1.0D0-ALPHA) &
		+ ALPHA*BEND_POT(TIJK,NI+1) 
		END IF 
		VANGLE = VANGLE + VIJK


		FXIJK  = FIJK * RXN/RIJ
                FYIJK  = FIJK * RYN/RIJ
                FZIJK  = FIJK * RZN/RIJ

                if(i .le. 765) write(9301,*) i, j, k, ti, tj, tk, fxijk, theta 

                FXI   = FXI + FXIJK
                FYI   = FYI + FYIJK
                FZI   = FZI + FZIJK

                FX(K) = FX(K) + FIJK* RXS/RKJ
                FY(K) = FY(K) + FIJK* RYS/RKJ
                FZ(K) = FZ(K) + FIJK* RZS/RKJ

                FX(J) = FX(J) - FIJK* RXS/RKJ - FXIJK
                FY(J) = FY(J) - FIJK* RYS/RKJ - FYIJK
                FZ(J) = FZ(J) - FIJK* RZS/RKJ - FZIJK

298	CONTINUE

        RETURN
        END

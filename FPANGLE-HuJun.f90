!      
!	Modified by Nicodemo
!	Jan. 2011
!
!	Added: lines 69-73    

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

		RXKJ = SX(J) - SX(K)
            RYKJ = SY(J) - SY(K)
            RZKJ = SZ(J) - SZ(K)

		RIJ = SQRT(RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0)
		RKJ = SQRT(RXKJ ** 2.0D0 + RYKJ ** 2.0D0 + RZKJ ** 2.0D0)
		CT = (RXIJ*RXKJ + RYIJ*RYKJ + RZIJ*RZKJ)/RIJ/RKJ

            if(CT .lt. -1.0D0)then
            write(3000,*)'CT PICCOLO',ct,i,j,k,timestepcheck
                  CT=-1.0D0
                  THETA=180.D0
            else if(CT .gt. 1.0D0) then
            write(3000,*)'CT GRANDE',ct,i,j,k,timestepcheck
                  CT=1.0D0
                  THETA=0.D0
            else
		      THETA = DACOS(CT) * R2D
            end if

!	****** vector M cross product -IJxKJ	
		RXM = - ( RYIJ*RZKJ - RZIJ*RYKJ)
		RYM = - (-RXIJ*RZKJ + RZIJ*RXKJ)
		RZM = - ( RXIJ*RYKJ - RYIJ*RXKJ)
		RM = SQRT(RXM ** 2.0D0 + RYM ** 2.0D0 + RZM ** 2.0D0)

      if(RM .eq. 0.) then

            RXIJ = RXI*0.99D0 - SX(J)
		RXM = - ( RYIJ*RZKJ - RZIJ*RYKJ)
		RYM = - (-RXIJ*RZKJ + RZIJ*RXKJ)
		RZM = - ( RXIJ*RYKJ - RYIJ*RXKJ)
		RM = SQRT(RXM ** 2.0D0 + RYM ** 2.0D0 + RZM ** 2.0D0)

            RXM = RXM / RM
		RYM = RYM / RM
		RZM = RZM / RM

      else

		RXM = RXM / RM
		RYM = RYM / RM
		RZM = RZM / RM

      end if

!	***** vector N cross product MxIJ
		RXN = ( RYM*RZIJ - RZM*RYIJ )/RIJ
		RYN = (-RXM*RZIJ + RZM*RXIJ )/RIJ
		RZN = ( RXM*RYIJ - RYM*RXIJ )/RIJ

!	***** vector S cross product -MxKJ
		RXS = -( RYM*RZKJ - RZM*RYKJ )/RKJ
		RYS = -(-RXM*RZKJ + RZM*RXKJ )/RKJ
		RZS = -( RXM*RYKJ - RYM*RXKJ )/RKJ

		NI = INT (THETA / BINA(TIJK))

		IF(NI .GT. NDATAN(TIJK)) THEN
			WRITE(*,*)'FATAL ERROR: Entry in angle table', TIJK,' does not exist'
			WRITE(1,*)'FATAL ERROR: Entry in angle table', TIJK,' does not exist'		
			STOP
		END IF

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
!		VANGLE = VANGLE + VIJK

                    if(type_label(i) .eq. 1)then
                        if(type_label(j) .eq. 2)then
                            VANGLE_MIX = VANGLE_MIX + VIJK              
                        else
                            if(type_label(k) .eq. 2)then
                                VANGLE_MIX = VANGLE_MIX + VIJK 
                            else
                                VANGLE = VANGLE + VIJK
                            end if
                        end if
                    elseif(type_label(i) .eq. 2)then
                        if(type_label(j) .eq. 1)then
                            VANGLE_MIX = VANGLE_MIX + VIJK                
                        else
                            if(type_label(k) .eq. 1 )then
                                VANGLE_MIX = VANGLE_MIX + VIJK
                            else
                                VANGLE_CG = VANGLE_CG + VIJK
                            end if
                        end if
                    end if 


		FXIJK  = FIJK * RXN/RIJ
                FYIJK  = FIJK * RYN/RIJ
                FZIJK  = FIJK * RZN/RIJ

                FXI   = FXI + FXIJK
                FYI   = FYI + FYIJK
                FZI   = FZI + FZIJK

                FX(K) = FX(K) - FIJK* RXS/RKJ
                FY(K) = FY(K) - FIJK* RYS/RKJ
                FZ(K) = FZ(K) - FIJK* RZS/RKJ

                FX(J) = FX(J) + FIJK* RXS/RKJ - FXIJK
                FY(J) = FY(J) + FIJK* RYS/RKJ - FYIJK
                FZ(J) = FZ(J) + FIJK* RZS/RKJ - FZIJK

298	CONTINUE

        RETURN
        END

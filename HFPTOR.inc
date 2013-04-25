
DO  M = 1, NIJKL(I) !296
		
    J = JTORIJKL(I, M)
    K = KTORIJKL(I, M)
    L = LTORIJKL(I, M)

    TJ = ITYPE(J)
    TK = ITYPE(K)
    TL = ITYPE(L)
    TIJKL = ITORT(TI, TJ, TK, TL)

    RXIJ = RXI - SX(J)
    RYIJ = RYI - SY(J)
    RZIJ = RZI - SZ(J)

    RXJK = SX(J) - SX(K)
    RYJK = SY(J) - SY(K)
    RZJK = SZ(J) - SZ(K)

    RXKL = SX(K) - SX(L)
    RYKL = SY(K) - SY(L)
    RZKL = SZ(K) - SZ(L)

!		RXKL = SX(L) - SX(K)
!              	RYKL = SY(L) - SY(K)
!              	RZKL = SZ(L) - SZ(K)

    RIJ = SQRT(RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0)
    RJK = SQRT(RXJK ** 2.0D0 + RYJK ** 2.0D0 + RZJK ** 2.0D0)
    RKL = SQRT(RXKL ** 2.0D0 + RYKL ** 2.0D0 + RZKL ** 2.0D0)

!	****** vector M cross product IJxJK
    RXM = RYIJ*RZJK - RZIJ*RYJK
    RYM =-RXIJ*RZJK + RZIJ*RXJK
    RZM = RXIJ*RYJK - RYIJ*RXJK
    RM = SQRT(RXM ** 2.0D0 + RYM ** 2.0D0 + RZM ** 2.0D0)

!	if RM = 0 IJ & JK are parallel

    if (RM.GT.0) THEN
        RXM = RXM / RM
        RYM = RYM / RM
        RZM = RZM / RM                
        COSM = (RXIJ*RXJK+RYIJ*RYJK+RZIJ*RZJK)/RIJ/RJK
        SINM = SQRT(1-COSM*COSM)
        COTANM = COSM/SINM
    elseif (RM .EQ. 0) then                
        RXIJ = RXI*0.99d0 - SX(J)       
        RIJ = SQRT(RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0)              
        RXM = RYIJ*RZJK - RZIJ*RYJK
        RYM =-RXIJ*RZJK + RZIJ*RXJK
        RZM = RXIJ*RYJK - RYIJ*RXJK
        RM = SQRT(RXM ** 2.0D0 + RYM ** 2.0D0 + RZM ** 2.0D0)
        RXM = RXM / RM
        RYM = RYM / RM
        RZM = RZM / RM
        COSM = (RXIJ*RXJK+RYIJ*RYJK+RZIJ*RZJK)/RIJ/RJK
        SINM = SQRT(1-COSM*COSM)
        COTANM = COSM/SINM         
    end if

!	****** vector N cross product JKxKL
    RXN = RYJK*RZKL - RZJK*RYKL
    RYN =-RXJK*RZKL + RZJK*RXKL
    RZN = RXJK*RYKL - RYJK*RXKL
    RN = SQRT(RXN ** 2.0D0 + RYN ** 2.0D0 + RZN ** 2.0D0)

    if (RN .GT. 0) then
        RXN = RXN / RN
        RYN = RYN / RN
        RZN = RZN / RN
        COSN = (RXJK*RXKL+RYJK*RYKL+RZJK*RZKL)/RJK/RKL
        SINN = SQRT(1-COSN*COSN)
        COTANN = COSN/SINN
    elseif (RN .EQ. 0)then
        RXKL = SX(K) - SX(L)*0.99d0
        RKL = SQRT(RXKL ** 2.0D0 + RYKL ** 2.0D0 + RZKL ** 2.0D0)                
        RXN = RYJK*RZKL - RZJK*RYKL
        RYN =-RXJK*RZKL + RZJK*RXKL
        RZN = RXJK*RYKL - RYJK*RXKL
        RN = SQRT(RXN ** 2.0D0 + RYN ** 2.0D0 + RZN ** 2.0D0)
        RXN = RXN / RN
        RYN = RYN / RN
        RZN = RZN / RN
        COSN = (RXJK*RXKL+RYJK*RYKL+RZJK*RZKL)/RJK/RKL
        SINN = SQRT(1-COSN*COSN)
        COTANN = COSN/SINN               
    endif

    COST = RXM*RXN + RYM*RYN + RZM*RZN
    SIGNT = -(RXM*RXKL + RYM*RYKL + RZM*RZKL)

    if(COST .GT. 1.0) then
        COST = 1.0
    else if (COST .LT. -1.0) THEN
        COST = -1.0
    end if

    PHI = DACOS(COST)*R2D

    IF (SIGNT.LT.0) PHI=360.D0 - PHI

    NI = INT (PHI / BINT(TIJKL))

    IF(NI .GT. NDATT(TIJKL)) THEN
        WRITE(*,*)'FATAL ERROR: Entry in torsion table', TIJKL,' does not exist'
        WRITE(1,*)'FATAL ERROR: Entry in torsion table', TIJKL,' does not exist'		
        STOP
    END IF

    IF (PHI == 360.0) THEN
        FIJKL = TOR_FORCE(TIJKL,NI)
        VIJKL = TOR_POT(TIJKL,NI)
    ELSE
        ALPHA=(PHI-ANGLE_TOR(TIJKL,NI))/BINT(TIJKL)
        FIJKL = TOR_FORCE(TIJKL,NI)*(1.0D0-ALPHA) &
                + ALPHA*TOR_FORCE(TIJKL,NI+1)
        VIJKL = TOR_POT(TIJKL,NI)*(1.0D0-ALPHA) &
                + ALPHA*TOR_POT(TIJKL,NI+1)
    END  IF
!		VTOR = VTOR + VIJKL

    if(type_label(i) .eq. 1)then
        if(type_label(j) .eq. 2)then
            VTOR_MIX = VTOR_MIX + VIJKL              
        else
            if(type_label(k) .eq. 2)then
                VTOR_MIX = VTOR_MIX + VIJKL 
            else
                if(type_label(l) .eq. 2)then
                    VTOR_MIX = VTOR_MIX + VIJKL
                else
                    VTOR = VTOR + VIJKL
                end if
            end if
        end if
    elseif(type_label(i) .eq. 2)then
        if(type_label(j) .eq. 1)then
            VTOR_MIX = VTOR_MIX + VIJKL                
        else
            if(type_label(k) .eq. 1 )then
                VTOR_MIX = VTOR_MIX + VIJKL
            else
                if(type_label(l) .eq. 1)then
                    VTOR_MIX = VTOR_MIX + VIJKL
                else
                    VTOR_CG = VTOR_CG + VIJKL
                end if
            end if
        end if
    end if 

    FX1 = -FIJKL * RXM/SINM/RIJ
    FY1 = -FIJKL * RYM/SINM/RIJ
    FZ1 = -FIJKL * RZM/SINM/RIJ

    FXI = FXI + FX1
    FYI = FYI + FY1
    FZI = FZI + FZ1
		
    FX4 = FIJKL * RXN/SINN/RKL
    FY4 = FIJKL * RYN/SINN/RKL
    FZ4 = FIJKL * RZN/SINN/RKL
		
    FXL(L) = FXL(L) + FX4
    FYL(L) = FYL(L) + FY4
    FZL(L) = FZL(L) + FZ4

    FX12 = FIJKL*(COTANM*RXM+COTANN*RXN)/RJK
    FY12 = FIJKL*(COTANM*RYM+COTANN*RYN)/RJK
    FZ12 = FIJKL*(COTANM*RZM+COTANN*RZN)/RJK

    FXL(J) = FXL(J) - FX1 + FX12
    FYL(J) = FYL(J) - FY1 + FY12
    FZL(J) = FZL(J) - FZ1 + FZ12

    FXL(K) = FXL(K) - FX4 - FX12
    FYL(K) = FYL(K) - FY4 - FY12
    FZL(K) = FZL(K) - FZ4 - FZ12


END DO 
!296	CONTINUE
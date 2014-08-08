DO  A = 1, NIJKL(I)
		
    J = JTORIJKL(I, A)
    K = KTORIJKL(I, A)
    L = LTORIJKL(I, A)

    TJ = ITYPE(J)
    TK = ITYPE(K)
    TL = ITYPE(L)
    TIJKL = ITORT(TI, TJ, TK, TL)

    RXYZ_IJ(1) = RXYZ_I(1) - SXYZ(1,J)
    RXYZ_IJ(2) = RXYZ_I(2) - SXYZ(2,J)
    RXYZ_IJ(3) = RXYZ_I(3) - SXYZ(3,J)

    RXYZ_JK(1) = SXYZ(1,J) - SXYZ(1,K)
    RXYZ_JK(2) = SXYZ(2,J) - SXYZ(2,K)
    RXYZ_JK(3) = SXYZ(3,J) - SXYZ(3,K)

    RXYZ_KL(1) = SXYZ(1,K) - SXYZ(1,L)
    RXYZ_KL(2) = SXYZ(2,K) - SXYZ(2,L)
    RXYZ_KL(3) = SXYZ(3,K) - SXYZ(3,L)

!		RXYZ_KL(1) = SXYZ(1,L) - SXYZ(1,K)
!              	RXYZ_KL(2) = SXYZ(2,L) - SXYZ(2,K)
!              	RXYZ_KL(3) = SXYZ(3,L) - SXYZ(3,K)

    RIJ = SQRT(RXYZ_IJ(1) ** 2.0 + RXYZ_IJ(2) ** 2.0 + RXYZ_IJ(3) ** 2.0)
    RJK = SQRT(RXYZ_JK(1) ** 2.0 + RXYZ_JK(2) ** 2.0 + RXYZ_JK(3) ** 2.0)
    RKL = SQRT(RXYZ_KL(1) ** 2.0 + RXYZ_KL(2) ** 2.0 + RXYZ_KL(3) ** 2.0)

!	****** vector M cross product IJxJK
    RXYZ_M(1) = RXYZ_IJ(2)*RXYZ_JK(3) - RXYZ_IJ(3)*RXYZ_JK(2)
    RXYZ_M(2) =-RXYZ_IJ(1)*RXYZ_JK(3) + RXYZ_IJ(3)*RXYZ_JK(1)
    RXYZ_M(3) = RXYZ_IJ(1)*RXYZ_JK(2) - RXYZ_IJ(2)*RXYZ_JK(1)
    RM = SQRT(RXYZ_M(1) ** 2.0 + RXYZ_M(2) ** 2.0 + RXYZ_M(3) ** 2.0)

!	if RM = 0 IJ & JK are parallel

    if (RM.GT.0) THEN
        RXYZ_M(1) = RXYZ_M(1) / RM
        RXYZ_M(2) = RXYZ_M(2) / RM
        RXYZ_M(3) = RXYZ_M(3) / RM                
        COSM = (RXYZ_IJ(1)*RXYZ_JK(1)+RXYZ_IJ(2)*RXYZ_JK(2)+RXYZ_IJ(3)*RXYZ_JK(3))/RIJ/RJK
        SINM = SQRT(1-COSM*COSM)
        COTANM = COSM/SINM
    else if (RM .EQ. 0) then                
        RXYZ_IJ(1) = RXYZ_I(1)*0.99 - SXYZ(1,J)       
        RIJ = SQRT(RXYZ_IJ(1) ** 2.0 + RXYZ_IJ(2) ** 2.0 + RXYZ_IJ(3) ** 2.0)              
        RXYZ_M(1) = RXYZ_IJ(2)*RXYZ_JK(3) - RXYZ_IJ(3)*RXYZ_JK(2)
        RXYZ_M(2) =-RXYZ_IJ(1)*RXYZ_JK(3) + RXYZ_IJ(3)*RXYZ_JK(1)
        RXYZ_M(3) = RXYZ_IJ(1)*RXYZ_JK(2) - RXYZ_IJ(2)*RXYZ_JK(1)
        RM = SQRT(RXYZ_M(1) ** 2.0 + RXYZ_M(2) ** 2.0 + RXYZ_M(3) ** 2.0)
        RXYZ_M(1) = RXYZ_M(1) / RM
        RXYZ_M(2) = RXYZ_M(2) / RM
        RXYZ_M(3) = RXYZ_M(3) / RM
        COSM = (RXYZ_IJ(1)*RXYZ_JK(1)+RXYZ_IJ(2)*RXYZ_JK(2)+RXYZ_IJ(3)*RXYZ_JK(3))/RIJ/RJK
        SINM = SQRT(1-COSM*COSM)
        COTANM = COSM/SINM         
    end if

!	****** vector N cross product JKxKL
    RXYZ_N(1) = RXYZ_JK(2)*RXYZ_KL(3) - RXYZ_JK(3)*RXYZ_KL(2)
    RXYZ_N(2) =-RXYZ_JK(1)*RXYZ_KL(3) + RXYZ_JK(3)*RXYZ_KL(1)
    RXYZ_N(3) = RXYZ_JK(1)*RXYZ_KL(2) - RXYZ_JK(2)*RXYZ_KL(1)
    RN = SQRT(RXYZ_N(1) ** 2.0 + RXYZ_N(2) ** 2.0 + RXYZ_N(3) ** 2.0)

    if (RN .GT. 0) then
        RXYZ_N(1) = RXYZ_N(1) / RN
        RXYZ_N(2) = RXYZ_N(2) / RN
        RXYZ_N(3) = RXYZ_N(3) / RN
        COSN = (RXYZ_JK(1)*RXYZ_KL(1)+RXYZ_JK(2)*RXYZ_KL(2)+RXYZ_JK(3)*RXYZ_KL(3))/RJK/RKL
        SINN = SQRT(1-COSN*COSN)
        COTANN = COSN/SINN
    elseif (RN .EQ. 0)then
        RXYZ_KL(1) = SXYZ(1,K) - SXYZ(1,L)*0.99
        RKL = SQRT(RXYZ_KL(1) ** 2.0 + RXYZ_KL(2) ** 2.0 + RXYZ_KL(3) ** 2.0)                
        RXYZ_N(1) = RXYZ_JK(2)*RXYZ_KL(3) - RXYZ_JK(3)*RXYZ_KL(2)
        RXYZ_N(2) =-RXYZ_JK(1)*RXYZ_KL(3) + RXYZ_JK(3)*RXYZ_KL(1)
        RXYZ_N(3) = RXYZ_JK(1)*RXYZ_KL(2) - RXYZ_JK(2)*RXYZ_KL(1)
        RN = SQRT(RXYZ_N(1) ** 2.0 + RXYZ_N(2) ** 2.0 + RXYZ_N(3) ** 2.0)
        RXYZ_N(1) = RXYZ_N(1) / RN
        RXYZ_N(2) = RXYZ_N(2) / RN
        RXYZ_N(3) = RXYZ_N(3) / RN
        COSN = (RXYZ_JK(1)*RXYZ_KL(1)+RXYZ_JK(2)*RXYZ_KL(2)+RXYZ_JK(3)*RXYZ_KL(3))/RJK/RKL
        SINN = SQRT(1-COSN*COSN)
        COTANN = COSN/SINN               
    endif

    COST = RXYZ_M(1)*RXYZ_N(1) + RXYZ_M(2)*RXYZ_N(2) + RXYZ_M(3)*RXYZ_N(3)
    SIGNT = -(RXYZ_M(1)*RXYZ_KL(1) + RXYZ_M(2)*RXYZ_KL(2) + RXYZ_M(3)*RXYZ_KL(3))

    if(COST .GT. 1.0) then
        COST = 1.0
    else if (COST .LT. -1.0) THEN
        COST = -1.0
    end if

    PHI = ACOS(COST)*R2D

    IF (SIGNT.LT.0) PHI=360. - PHI

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
        FIJKL = TOR_FORCE(TIJKL,NI)*(1.0-ALPHA) &
                + ALPHA*TOR_FORCE(TIJKL,NI+1)
        VIJKL = TOR_POT(TIJKL,NI)*(1.0-ALPHA) &
                + ALPHA*TOR_POT(TIJKL,NI+1)
    END  IF

    V_TORSION(TORSION_TYPE_LABEL(TIJKL)) = V_TORSION(TORSION_TYPE_LABEL(TIJKL)) + VIJKL

    FXYZ_1(1) = -FIJKL * RXYZ_M(1)/SINM/RIJ
    FXYZ_1(2) = -FIJKL * RXYZ_M(2)/SINM/RIJ
    FXYZ_1(3) = -FIJKL * RXYZ_M(3)/SINM/RIJ

    FXYZ_I(1) = FXYZ_I(1) + FXYZ_1(1)
    FXYZ_I(2) = FXYZ_I(2) + FXYZ_1(2)
    FXYZ_I(3) = FXYZ_I(3) + FXYZ_1(3)

    FXYZ_4(1) = FIJKL * RXYZ_N(1)/SINN/RKL
    FXYZ_4(2) = FIJKL * RXYZ_N(2)/SINN/RKL
    FXYZ_4(3) = FIJKL * RXYZ_N(3)/SINN/RKL
		
    FXYZ_L(1,L) = FXYZ_L(1,L) + FXYZ_4(1)
    FXYZ_L(2,L) = FXYZ_L(2,L) + FXYZ_4(2)
    FXYZ_L(3,L) = FXYZ_L(3,L) + FXYZ_4(3)

    FXYZ_12(1) = FIJKL*(COTANM*RXYZ_M(1)+COTANN*RXYZ_N(1))/RJK
    FXYZ_12(2) = FIJKL*(COTANM*RXYZ_M(2)+COTANN*RXYZ_N(2))/RJK
    FXYZ_12(3) = FIJKL*(COTANM*RXYZ_M(3)+COTANN*RXYZ_N(3))/RJK

    FXYZ_L(1,J) = FXYZ_L(1,J) - FXYZ_1(1) + FXYZ_12(1)
    FXYZ_L(2,J) = FXYZ_L(2,J) - FXYZ_1(2) + FXYZ_12(2)
    FXYZ_L(3,J) = FXYZ_L(3,J) - FXYZ_1(3) + FXYZ_12(3)

    FXYZ_L(1,K) = FXYZ_L(1,K) - FXYZ_4(1) - FXYZ_12(1)
    FXYZ_L(2,K) = FXYZ_L(2,K) - FXYZ_4(2) - FXYZ_12(2)
    FXYZ_L(3,K) = FXYZ_L(3,K) - FXYZ_4(3) - FXYZ_12(3)

END DO 

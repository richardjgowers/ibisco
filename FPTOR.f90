


SUBROUTINE FPTOR (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)
use VAR
implicit none

integer :: J, K, M, L 
integer,intent(in) :: i,ti
integer ::  TJ,TK,TL,TIJ,NI, TIJK, TIJKL
real(kind=rkind),intent(inout) :: RXI, RYI, RZI, FXI, FYI, FZI
real(kind=rkind) :: RCUTSQ, W, ALPHA, ct
real(kind=rkind) :: RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
real(kind=rkind) :: RXKJ, RYKJ, RZKJ, RKJSQ, THETA
real(kind=rkind) :: RXM, RYM, RZM, RM, RXN, RYN, RZN, RJI, RKI
real(kind=rkind) :: RXS, RYS, RZS
real(kind=rkind) :: VIJ, WIJ, FIJ, RIJ, FIJK, RKJ, FXIJK, FYIJK, FZIJK
real(kind=rkind) :: FXIJKL, FYIJKL, FZIJKL, FIJKL, PHI, COST, SIGNT, COSM, SINM
real(kind=rkind) :: RXJK, RYJK, RZJK, RXKL, RYKL, RZKL, RJK, RKL, RN
real(kind=rkind) :: COTANN, COTANM, COSN, SINN
real(kind=rkind) :: RXIK, RYIK, RZIK, RIK, VIJK, VIJKL
real(kind=rkind) :: FX1, FY1, FZ1, FX4, FY4, FZ4, FX12, FY12, FZ12
!      ******************************************************************

DO 296 M = 1, NIJKL(I)
		
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

    VTOR = VTOR + VIJKL

    FX1 = -FIJKL * RXM/SINM/RIJ
    FY1 = -FIJKL * RYM/SINM/RIJ
    FZ1 = -FIJKL * RZM/SINM/RIJ

    FXI = FXI + FX1
    FYI = FYI + FY1
    FZI = FZI + FZ1
		
    FX4 = FIJKL * RXN/SINN/RKL
    FY4 = FIJKL * RYN/SINN/RKL
    FZ4 = FIJKL * RZN/SINN/RKL
		
    FX(L) = FX(L) + FX4
    FY(L) = FY(L) + FY4
    FZ(L) = FZ(L) + FZ4

    FX12 = FIJKL*(COTANM*RXM+COTANN*RXN)/RJK
    FY12 = FIJKL*(COTANM*RYM+COTANN*RYN)/RJK
    FZ12 = FIJKL*(COTANM*RZM+COTANN*RZN)/RJK

    FX(J) = FX(J) - FX1 + FX12
    FY(J) = FY(J) - FY1 + FY12
    FZ(J) = FZ(J) - FZ1 + FZ12

    FX(K) = FX(K) - FX4 - FX12
    FY(K) = FY(K) - FY4 - FY12
    FZ(K) = FZ(K) - FZ4 - FZ12

296	CONTINUE

        RETURN
        END
!	*********************************************************************************************

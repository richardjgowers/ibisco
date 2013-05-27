  DO A = 1, NOOPIJKL(I)

     J = JOOPIJKL(I, A)
     K = KOOPIJKL(I, A)
     L = LOOPIJKL(I, A)

     TJ = ITYPE(J)
     TK = ITYPE(K)
     TL = ITYPE(L)
     TOIJKL = IOOPT(TI,TJ,TK,TL)

     RXIJ = RXI - SX(J) 
     RYIJ = RYI - SY(J) 
     RZIJ = RZI - SZ(J) 
     RXJK = SX(K) - SX(J) 
     RYJK = SY(K) - SY(J) 
     RZJK = SZ(K) - SZ(J) 
     RXKL = SX(K) - SX(L) 
     RYKL = SY(K) - SY(L) 
     RZKL = SZ(K) - SZ(L)

     RIJ = SQRT(RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0)
     RJK = SQRT(RXJK ** 2.0D0 + RYJK ** 2.0D0 + RZJK ** 2.0D0)
     RKL = SQRT(RXKL ** 2.0D0 + RYKL ** 2.0D0 + RZKL ** 2.0D0)

     !	****** vector M cross product IJxJK
     RXM = RYIJ*RZJK - RZIJ*RYJK
     RYM =-RXIJ*RZJK + RZIJ*RXJK
     RZM = RXIJ*RYJK - RYIJ*RXJK
     RM = SQRT(RXM ** 2.0D0 + RYM ** 2.0D0 + RZM ** 2.0D0)
     if (RM .lt. 1.0D-08) stop 'rm = 0'
     RXM = RXM / RM
     RYM = RYM / RM
     RZM = RZM / RM                
     COSM = (RXIJ*RXJK+RYIJ*RYJK+RZIJ*RZJK)/RIJ/RJK
     SINM = SQRT(1-COSM*COSM)
     COTANM = COSM/SINM

     !	****** vector N cross product JKxKL
     RXN = RYJK*RZKL - RZJK*RYKL
     RYN =-RXJK*RZKL + RZJK*RXKL
     RZN = RXJK*RYKL - RYJK*RXKL
     RN = SQRT(RXN ** 2.0D0 + RYN ** 2.0D0 + RZN ** 2.0D0)
     if (RN .lt. 1.0D-08) stop 'rn = 0'
     RXN = RXN / RN
     RYN = RYN / RN
     RZN = RZN / RN

     COSN = (RXJK*RXKL+RYJK*RYKL+RZJK*RZKL)/RJK/RKL
     SINN = SQRT(1-COSN*COSN)
     COTANN = COSN/SINN

     COST = RXM*RXN + RYM*RYN + RZM*RZN
     RXmN = RYm*RZn - RZm*RYn
     RYmN =-RXm*RZn + RZm*RXn
     RZmN = RXm*RYn - RYm*RXn

     SIGNT = sign(1.0D0,RXMn*RXjk + RYMn*RYjk + RZMn*RZjk)


     IF (COST .GT. 1.0D0 + eps) THEN
        write(*,*) i,j,k,l
        write(*,*)COST 
        stop 'COST .GT. 1.00'
     else if(COST .lT. -(1.00D0+eps))THEN
        write(*,*) i,j,k,l
        write(*,*)COST 
        stop 'COST .lT. 1.00'
     END IF

     IF (COST .GT. 1.0D0 ) THEN
        COST=1.0D0 
     else if(COST .lT. -1.00D0)THEN
        COST=-1.0D0 
     END IF

     PHI = DACOS(COST)
     phi = signt * phi
     !    phi = phi - pi*0.5D0 * anint(phi * ipi*2.0D0)
     phi = phi*R2D

     PHI_t= 180.0D0 + phi
     NI = anint(PHI_t / BINO(TOIJKL))

     if(NI .lt. 0)then
        write(*,*) i,j,k,l,signt
        write(*,*) NI,phi_t
        write(*,*)PHI_t / BINO(TOIJKL)
        stop
     end if

     IF(NI .GT. NDATO(TOIJKL)) THEN
        WRITE(*,*)'FATAL ERROR: Entry in out of plane table', TOIJKL,' does not exist'
        WRITE(1,*)'FATAL ERROR: Entry in out of plane table', TOIJKL,' does not exist'
        write(*,*) i,j,k,l
        write(*,*) NI,phi
        WRITE(*,*)'Simulation stopped at Time Step: ',STEP
        STOP
     END IF

     ALPHA=abs(abs(PHI)-abs(ANGLE_OOP(TOIJKL,NI)))/BINO(TOIJKL)
     if(alpha .gt. 1.0D0)then
        ni = ni+1
        ALPHA=abs(abs(PHI)-abs(ANGLE_OOP(TOIJKL,NI)))/BINO(TOIJKL)
     end if

     FIJKL = OOP_FORCE(TOIJKL,NI)*(1.0D0-ALPHA)+ ALPHA*OOP_FORCE(TOIJKL,NI+1)
     VIJKL = OOP_POT(TOIJKL,NI)*(1.0D0-ALPHA)+ ALPHA*OOP_POT(TOIJKL,NI+1)

     VOOP = VOOP + VIJKL

     FX1 = FIJKL * RXM/SINM/RIJ
     FY1 = FIJKL * RYM/SINM/RIJ
     FZ1 = FIJKL * RZM/SINM/RIJ

     FXI = FXI + FX1
     FYI = FYI + FY1
     FZI = FZI + FZ1

     FX4 = -FIJKL * RXN/SINN/RKL
     FY4 = -FIJKL * RYN/SINN/RKL
     FZ4 = -FIJKL * RZN/SINN/RKL

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

  END DO

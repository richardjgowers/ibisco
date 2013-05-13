#include "ibi-preprocess.h"

SUBROUTINE FPOUTPLANE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)
	
USE VAR

IMPLICIT NONE

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
REAL*8	      RXJL, RYJL, RZJL, RJL
REAL*8	      COTANN, COTANM, COSN, SINN
REAL*8	      RXIK, RYIK, RZIK, RIK, VIJK, VIJKL
REAL(KIND=8) :: rx1,ry1,rz1,rx2,ry2,rz2,rx3,ry3,rz3,rx4,ry4,rz4,phi_t
INTEGER       JBEG, JEND, JNAB, TI, TJ, TIJ, NI
INTEGER	      TK, TL, TIJK, TOIJKL
REAL*8	      FX1, FY1, FZ1, FX4, FY4, FZ4, FX12, FY12, FZ12
real(kind=8) :: rxmn,rymn,rzmn,rmn
real(kind=8),parameter :: eps=1.0D-4
!      ******************************************************************

DO M = 1, NOOPIJKL(I)

    J = JOOPIJKL(I, M)
    K = KOOPIJKL(I, M)
    L = LOOPIJKL(I, M)

    TJ = ITYPE(J)
    TK = ITYPE(K)
    TL = ITYPE(L)
    TOIJKL = IOOPT(TI,TJ,TK,TL)

    rx1=RXI
    ry1=RYI
    rz1=RZI

    rx2=SX(j)
    ry2=SY(j)
    rz2=SZ(j)

    rx3=SX(k)
    ry3=SY(k)
    rz3=SZ(k)

    rx4=SX(l)
    ry4=SY(l)
    rz4=SZ(l)

    RXIJ = rx1 - rx2 
    RYIJ = ry1 - ry2 
    RZIJ = rz1 - rz2 
    RXJK = rx3 - rx2 
    RYJK = ry3 - ry2 
    RZJK = rz3 - rz2 
    RXKL = rx3 - rx4 
    RYKL = ry3 - ry4 
    RZKL = rz3 - rz4

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
                  WRITE(*,*)'Simulation stopped at Time Step: ',timestepcheck
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

#ifdef DEBUG_OOP
if(vijkl .lt. 0.0D0 .or. alpha .gt. 1.0D0 .or. alpha .lt. 0.0D0)then
    write(*,*)
    write(*,*)i,j,k,l,'NEGATIVO'
    write(*,*)VIJKL,ni,alpha,signt
    write(*,*)OOP_POT(TOIJKL,NI)
    write(*,*)OOP_POT(TOIJKL,NI+1)
    write(*,*)
    stop
end if


if(toijkl .le. 4)then
voop_LD =VOOP_LD + VIJKL
    if(toijkl .eq. 1 .or. toijkl .eq. 2)then
voop_L =VOOP_L + VIJKL
    elseif(toijkl .eq. 3 .or. toijkl .eq. 4)then
voop_D =VOOP_D + VIJKL
    end if

write(111,*)i,j,k,l
write(111,*)ti,tj,tk,tl
write(111,*)phi,alpha,ni
write(111,*)vijkl
write(111,*)


else
voop_ring =VOOP_ring + VIJKL
end if

#endif

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


        RETURN

END SUBROUTINE FPOUTPLANE


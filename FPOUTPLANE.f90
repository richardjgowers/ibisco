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
	REAL*8	      RXIL, RYIL, RZIL
REAL(KIND=8) :: rx1,ry1,rz1,rx2,ry2,rz2,rx3,ry3,rz3,rx4,ry4,rz4,phi_t
	INTEGER       JBEG, JEND, JNAB, TI, TJ, TIJ, NI
	INTEGER	      TK, TL, TIJK, TOIJKL
	REAL*8	      FX1, FY1, FZ1, FX4, FY4, FZ4, FX12, FY12, FZ12
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

    RXIJ = rx2 - rx1 
    RYIJ = ry2 - ry1 
    RZIJ = rz2 - rz1 
    RXJK = rx3 - rx2 
    RYJK = ry3 - ry2 
    RZJK = rz3 - rz2 
    RXKL = rx4 - rx3 
    RYKL = ry4 - ry3 
    RZKL = rz4 - rz3

! This vector allows to calculate the direction of the torsion

    RXIL = rx4 - rx1 
    RYIL = ry4 - ry1 
    RZIL = rz4 - rz1

		RIJ = SQRT(RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0)
		RJK = SQRT(RXJK ** 2.0D0 + RYJK ** 2.0D0 + RZJK ** 2.0D0)
		RKL = SQRT(RXKL ** 2.0D0 + RYKL ** 2.0D0 + RZKL ** 2.0D0)

!	****** vector M cross product IJxJK
		RXM = RYIJ*RZJK - RZIJ*RYJK
		RYM =-RXIJ*RZJK + RZIJ*RXJK
		RZM = RXIJ*RYJK - RYIJ*RXJK
		RM = SQRT(RXM ** 2.0D0 + RYM ** 2.0D0 + RZM ** 2.0D0)
                if (RM.GT.0) THEN
		RXM = RXM / RM
		RYM = RYM / RM
		RZM = RZM / RM
                
                COSM = (RXIJ*RXJK+RYIJ*RYJK+RZIJ*RZJK)/RIJ/RJK
                SINM = SQRT(1-COSM*COSM)
                COTANM = COSM/SINM

                ELSEIF (RM.EQ.0) THEN
                 
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

                ENDIF
!		COSM = (RXIJ*RXJK+RYIJ*RYJK+RZIJ*RZJK)/RIJ/RJK
!		SINM = SQRT(1-COSM*COSM)
!		COTANM = COSM/SINM
!	****** vector N cross product JKxKL
		RXN = RYJK*RZKL - RZJK*RYKL
		RYN =-RXJK*RZKL + RZJK*RXKL
		RZN = RXJK*RYKL - RYJK*RXKL
		RN = SQRT(RXN ** 2.0D0 + RYN ** 2.0D0 + RZN ** 2.0D0)

                if (RN.GT.0) THEN
		RXN = RXN / RN
		RYN = RYN / RN
		RZN = RZN / RN

		COSN = (RXJK*RXKL+RYJK*RYKL+RZJK*RZKL)/RJK/RKL
		SINN = SQRT(1-COSN*COSN)
		COTANN = COSN/SINN
                ELSEIF (RN.EQ.0) THEN
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

                ENDIF

		COST = RXM*RXN + RYM*RYN + RZM*RZN
		SIGNT = -(RXM*RXIL + RYM*RYIL + RZM*RZIL)

		IF (COST .GT. 1.00) THEN
			COST = 1.00
		END IF

		PHI = DACOS(COST)*R2D
        phi_t = phi
        !if(toijkl .le. 4) THEN		
		    IF (SIGNT .LT. 0) then
                PHI_t= 180.0D0 - PHI_t
            else
                PHI_t= 180.0D0 + PHI_t
            end if
		!END IF

		NI = INT (PHI_t / BINO(TOIJKL))

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
	

!		IF (PHI == 360.0) THEN
!		FIJKL = OOP_FORCE(TOIJKL,NI)
!		VIJKL = OOP_POT(TOIJKL,NI)

!		ELSE
		
		ALPHA=abs(PHI-abs(ANGLE_OOP(TOIJKL,NI)))/BINO(TOIJKL)

        if(alpha .gt. 1.0D0)then
            ni = ni+1
            ALPHA=abs(PHI-abs(ANGLE_OOP(TOIJKL,NI)))/BINO(TOIJKL)
        end if
 !write(2000,*)ti,'indx',i,alpha
!   write(2000,*)  PHI,ANGLE_OOP(TOIJKL,NI),BINO(TOIJKL)
! write(2000,*)

		FIJKL = OOP_FORCE(TOIJKL,NI)*(1.0D0-ALPHA) &
		+ ALPHA*OOP_FORCE(TOIJKL,NI+1)

		VIJKL = OOP_POT(TOIJKL,NI)*(1.0D0-ALPHA) &
		+ ALPHA*OOP_POT(TOIJKL,NI+1)

!		END  IF

		VOOP = VOOP + VIJKL

!if(toijkl .le. 4)then
!fijkl = 0.0

!voop_LD =VOOP_LD + VIJKL
!    if(toijkl .eq. 1 .or. toijkl .eq. 2)then
!voop_L =VOOP_L + VIJKL
!    elseif(toijkl .eq. 3 .or. toijkl .eq. 4)then
!voop_D =VOOP_D + VIJKL
!    end if
!else
!voop_ring =VOOP_ring + VIJKL
!end if

     !   if(TOIJKL .gt. 3)then
!    write(2000,*)ti,'indx',i,j,k,l
!    write(2000,*)signt,phi,ti,tj,tk,tl
!    write(2000,*)
	!END  IF
!if(i .eq. 5354)then 
!        write(2000,*)ti,'indx',i,j,k,l,timestepcheck
!    write(2000,*)phi_t, vijkl,signt
!    write(2000,*)

!end if

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

    FX12 = -FIJKL*(COTANM*RXM+COTANN*RXN)/RJK
    FY12 = -FIJKL*(COTANM*RYM+COTANN*RYN)/RJK
    FZ12 = -FIJKL*(COTANM*RZM+COTANN*RZN)/RJK

    FX(J) = FX(J) - FX1 + FX12
    FY(J) = FY(J) - FY1 + FY12
    FZ(J) = FZ(J) - FZ1 + FZ12

    FX(K) = FX(K) - FX4 - FX12
    FY(K) = FY(K) - FY4 - FY12
    FZ(K) = FZ(K) - FZ4 - FZ12

!end if
!end if
!end if
!end if
	END DO


        RETURN

END SUBROUTINE FPOUTPLANE


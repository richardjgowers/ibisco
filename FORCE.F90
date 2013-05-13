#include "ibi-preprocess.h"

SUBROUTINE FORCE ( )

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  integer :: I, J, K, M, L 
  integer :: JBEG, JEND, JNAB 
  INTEGER :: TI, TJ, TK, TL, TIJ, TIJK, TIJKL, NI, TOIJKL
  real(kind=rkind) :: RCUTSQ, ALPHA
  real(kind=rkind) :: RXI, RYI, RZI, FXI, FYI, FZI
  real(kind=rkind) :: RXIJ, RYIJ, RZIJ, RIJSQ, RIJ, FXIJ, FYIJ, FZIJ
  real(kind=rkind) :: VIJ, FIJ
  REAL(KIND=RKIND), dimension(NATOMS) :: FXL,FYL,FZL !Temp openmp variables
!FPANGLE
  REAL(KIND=RKIND) :: W, CT
  REAL(KIND=RKIND) :: RXKJ, RYKJ, RZKJ, RKJ, THETA
  REAL(KIND=RKIND) :: RXM,RYM,RZM,RM,RXN,RYN,RZN,RN,RXS,RYS,RZS
  REAL(KIND=RKIND) :: FIJK, VIJK
  REAL(KIND=RKIND) :: FXIJK, FYIJK, FZIJK
!FPTOR
  REAL(KIND=RKIND) :: RXKL,RYKL,RZKL,RKL
  REAL(KIND=RKIND) :: RXJK,RYJK,RZJK,RJK
  REAL(KIND=RKIND) :: COSM,SINM,COTANM,COSN,SINN,COTANN,COST,SIGNT,PHI
  REAL(KIND=RKIND) :: FIJKL,VIJKL
  REAL(KIND=RKIND) :: FX1, FY1, FZ1, FX4, FY4, FZ4, FX12, FY12, FZ12
!FPOOP
!  REAL(KIND=RKIND) :: rx1,ry1,rz1,rx2,ry2,rz2,rx3,ry3,rz3,rx4,ry4,rz4,phi_t
!  real(kind=RKIND) :: rxmn,rymn,rzmn,rmn
!  real(kind=8),parameter :: eps=1.0D-4
  !      *******************************************************************

  RCUTSQ = RCUT*RCUT

  VNBOND = 0.0D0
  VBOND  = 0.0D0
  VANGLE = 0.0D0
  VTOR   = 0.0D0
  VOOP   = 0.0D0
#ifdef DEBUG_OOP
      VOOP_ring   = 0.0D0
      VOOP_LD   = 0.0D0
      VOOP_L   = 0.0D0
      VOOP_D   = 0.0D0
#endif

  PT11 = 0.0
  PT22 = 0.0
  PT33 = 0.0
  PT12 = 0.0
  PT13 = 0.0
  PT23 = 0.0

!$OMP PARALLEL DEFAULT(NONE)&
!$OMP& SHARED(NATOMS,POINT,RX,RY,RZ,ITYPE,LIST,INBONDT)&
!$OMP& SHARED(BOXXINV,BOXYINV,BOXZINV,BOXX,BOXY,BOXZ)&
!$OMP& SHARED(RCUTSQ,BINNB,RNBOND,NBOND_FORCE,NBOND_POT)&
!$OMP& SHARED(SX,SY,SZ,timestepcheck,RCUT)&
!$OMP& SHARED(NBONDS,JBOND,typeBond,IBONDT)&
!$OMP& SHARED(NDATB,BINB,RBOND,BOND_FORCE,BOND_POT)&
!$OMP& SHARED(NIJK,JANGLEIJK,KANGLEIJK,IANGT)&
!$OMP& SHARED(BINA,ANGLE,BEND_FORCE,BEND_POT,R2D)&
!$OMP& SHARED(NIJKL,JTORIJKL,KTORIJKL,LTORIJKL,ITORT)&
!$OMP& SHARED(BINT,NDATT,TOR_FORCE,TOR_POT,ANGLE_TOR)&
!$OMP& PRIVATE(I,J,K,L,TI,TJ,TK,TL,M,TIJ,TIJK,TIJKL)&
!$OMP& PRIVATE(JBEG,JEND,JNAB,RXI,RYI,RZI,FXI,FYI,FZI)&
!$OMP& PRIVATE(RXIJ,RYIJ,RZIJ,RIJSQ,RIJ)&
!$OMP& PRIVATE(NI,ALPHA,CT,THETA)&
!$OMP& PRIVATE(FIJ,VIJ,FXIJ,FYIJ,FZIJ,FXL,FYL,FZL)&
!$OMP& PRIVATE(RXKJ,RYKJ,RZKJ,RKJ)&
!$OMP& PRIVATE(RXJK,RYJK,RZJK,RJK)&
!$OMP& PRIVATE(RXKL,RYKL,RZKL,RKL)&
!$OMP& PRIVATE(RXM,RYM,RZM,RM,RXN,RYN,RZN,RN,RXS,RYS,RZS)&
!$OMP& PRIVATE(FIJK,VIJK,FXIJK,FYIJK,FZIJK)&
!$OMP& PRIVATE(COSM,SINM,COTANM,COSN,SINN,COTANN,COST,SIGNT,PHI)&
!$OMP& PRIVATE(FIJKL,VIJKL)&
!$OMP& PRIVATE(FX1,FY1,FZ1,FX4,FY4,FZ4,FX12,FY12,FZ12)&
!$OMP& REDUCTION(+:VNBOND,VBOND,VANGLE,VTOR)&
!$OMP& REDUCTION(+:PT11,PT22,PT33,PT12,PT13,PT23)&
!$OMP& REDUCTION(+:FX,FY,FZ,FXNB,FYNB,FZNB)

!Local counters used to reduce the amount of communication between threads
!Local counters for force
DO I=1,NATOMS
   FXL(I) = 0.0
   FYL(I) = 0.0
   FZL(I) = 0.0
   FXNB(I) = 0.0
   FYNB(I) = 0.0
   FZNB(I) = 0.0
END DO

  !       ** USE THE LIST TO FIND THE NEIGHBOURS **
!$OMP DO SCHEDULE(STATIC,1)
  DO I = 1, NATOMS !200

     JBEG = POINT(I)
     JEND = POINT(I+1) - 1

     !       ** CHECK THAT ATOM I HAS NEIGHBOURS **

     IF( JBEG .LE. JEND ) THEN

        RXI = RX(I)
        RYI = RY(I)
        RZI = RZ(I)
        FXI = 0.0
        FYI = 0.0
        FZI = 0.0
        TI = ITYPE(I)

        DO JNAB = JBEG, JEND !199

           !		TAKE THE INDEX OF NEIGHBOUR ATOMS
           J = LIST(JNAB)

           !		TAKE THE TYPE OF NEIGHBOUR ATOM AND NON-BONDED INTERACTIONS	
           TJ = ITYPE(J)
           TIJ = INBONDT(TI, TJ)

           IF( TIJ .NE. 0) THEN	
              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J)
              RZIJ = RZI - RZ(J)

              RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
              RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
              RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ

              RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0

              IF ( RIJSQ < RCUTSQ ) THEN

                 RIJ = DSQRT(RIJSQ)
                 NI = INT (RIJ / BINNB(TIJ))

                 !		LINEAR INTEPOLATION
                 ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)

                 FIJ = NBOND_FORCE(TIJ,NI)*(1.0-ALPHA) &
                      + ALPHA*NBOND_FORCE(TIJ,NI+1) 

                 VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
                      + ALPHA*NBOND_POT(TIJ,NI+1)

                 VNBOND = VNBOND + VIJ

                 FXIJ  = FIJ * RXIJ
                 FYIJ  = FIJ * RYIJ
                 FZIJ  = FIJ * RZIJ

!                 IF (DPDINPUT.EQ.1.and.RIJ.LT.RCUTDPD) THEN
!                    IF (WDPDT.EQ.2) THEN
!                       call FDR(I,J,FXIJ,FYIJ,FZIJ,RXIJ,RYIJ,RZIJ,RIJ)
!                    ELSEIF (WDPDT.EQ.1) THEN
!                       call FDR_STEP(I,J,FXIJ,FYIJ,FZIJ,RXIJ,RYIJ,RZIJ,RIJ)
!                    ENDIF
!                 ENDIF

                 FXI   = FXI + FXIJ
                 FYI   = FYI + FYIJ
                 FZI   = FZI + FZIJ

                 FXL(J) = FXL(J) - FXIJ
                 FYL(J) = FYL(J) - FYIJ
                 FZL(J) = FZL(J) - FZIJ

                 !		ADD THE NON-BONDED PART OF PRESSURE
                 PT11 = PT11 + FXIJ * RXIJ
                 PT22 = PT22 + FYIJ * RYIJ
                 PT33 = PT33 + FZIJ * RZIJ

                 PT12 = PT12 + FYIJ * RXIJ
                 PT13 = PT13 + FZIJ * RXIJ
                 PT23 = PT23 + FZIJ * RYIJ

                 !***************LOWE ANDERSON***************************************
!                 IF (LAINPUT.EQ.1.AND.RIJ.LT.RCUTDPD) THEN
!                    call LOWEAND (RIJSQ,RXIJ,RYIJ,RZIJ,I,J)
!                 ENDIF
                 !***************LOWE ANDERSON***************************************
              END IF   ! endif  IF ( RIJSQ < RCUTSQ )
           END IF     ! endif   IF( TIJ .NE. 0) THEN
        END DO !199        CONTINUE

!        IF (DPD_BONDED.EQ.1) THEN
!           IF (LAINPUT.EQ.1) THEN
!              call  LA_ADDITION (I,RXI,RYI,RZI)
!           ELSEIF (DPDINPUT.EQ.1) THEN
!              IF(WDPDT.EQ.2) THEN
!                call FDR_ADDITION (I,RXI,RYI,RZI,FXI,FYI,FZI )
!              ELSEIF(WDPDT.EQ.1) THEN
!                 call FDR_STEP_ADDITION (I,RXI,RYI,RZI,FXI,FYI,FZI )
!              ENDIF
!           ENDIF
!        ENDIF

        FXL(I) = FXL(I) + FXI
        FYL(I) = FYL(I) + FYI
        FZL(I) = FZL(I) + FZI

     ENDIF
  END DO !200     CONTINUE
!$OMP END DO

  !Collate forces from all threads
DO I=1,NATOMS
   FX(I) = FX(I) + FXL(I)
   FY(I) = FY(I) + FYL(I) 
   FZ(I) = FZ(I) + FZL(I)
END DO

!$OMP END PARALLEL

!SAVE THE NON-BONDED PART OF FORCE
DO I=1, NATOMS
   FXNB(I) = FX(I)
   FYNB(I) = FY(I)
   FZNB(I) = FZ(I)
END DO

!$OMP PARALLEL DEFAULT(NONE)&
!$OMP& SHARED(NATOMS,POINT,RX,RY,RZ,ITYPE,LIST,INBONDT)&
!$OMP& SHARED(BOXXINV,BOXYINV,BOXZINV,BOXX,BOXY,BOXZ)&
!$OMP& SHARED(RCUTSQ,BINNB,RNBOND,NBOND_FORCE,NBOND_POT)&
!$OMP& SHARED(SX,SY,SZ,timestepcheck,RCUT)&
!$OMP& SHARED(NBONDS,JBOND,typeBond,IBONDT)&
!$OMP& SHARED(NDATB,BINB,RBOND,BOND_FORCE,BOND_POT)&
!$OMP& SHARED(NIJK,JANGLEIJK,KANGLEIJK,IANGT)&
!$OMP& SHARED(BINA,ANGLE,BEND_FORCE,BEND_POT,R2D)&
!$OMP& SHARED(NIJKL,JTORIJKL,KTORIJKL,LTORIJKL,ITORT)&
!$OMP& SHARED(BINT,NDATT,TOR_FORCE,TOR_POT,ANGLE_TOR)&
!$OMP& PRIVATE(I,J,K,L,TI,TJ,TK,TL,M,TIJ,TIJK,TIJKL)&
!$OMP& PRIVATE(JBEG,JEND,JNAB,RXI,RYI,RZI,FXI,FYI,FZI)&
!$OMP& PRIVATE(RXIJ,RYIJ,RZIJ,RIJSQ,RIJ)&
!$OMP& PRIVATE(NI,ALPHA,CT,THETA)&
!$OMP& PRIVATE(FIJ,VIJ,FXIJ,FYIJ,FZIJ,FXL,FYL,FZL)&
!$OMP& PRIVATE(RXKJ,RYKJ,RZKJ,RKJ)&
!$OMP& PRIVATE(RXJK,RYJK,RZJK,RJK)&
!$OMP& PRIVATE(RXKL,RYKL,RZKL,RKL)&
!$OMP& PRIVATE(RXM,RYM,RZM,RM,RXN,RYN,RZN,RN,RXS,RYS,RZS)&
!$OMP& PRIVATE(FIJK,VIJK,FXIJK,FYIJK,FZIJK)&
!$OMP& PRIVATE(COSM,SINM,COTANM,COSN,SINN,COTANN,COST,SIGNT,PHI)&
!$OMP& PRIVATE(FIJKL,VIJKL)&
!$OMP& PRIVATE(FX1,FY1,FZ1,FX4,FY4,FZ4,FX12,FY12,FZ12)&
!$OMP& REDUCTION(+:VNBOND,VBOND,VANGLE,VTOR)&
!$OMP& REDUCTION(+:PT11,PT22,PT33,PT12,PT13,PT23)&
!$OMP& REDUCTION(+:FX,FY,FZ,FXNB,FYNB,FZNB)

!Reset local force counters
DO I=1,NATOMS
   FXL(I) = 0.0
   FYL(I) = 0.0
   FZL(I) = 0.0
END DO

!$OMP DO SCHEDULE(STATIC,1)
  DO I = 1, NATOMS
     RXI = SX(I)
     RYI = SY(I)
     RZI = SZ(I)
     FXI = 0.0
     FYI = 0.0
     FZI = 0.0
     TI = ITYPE(I)

     !	*******************************************************************************************
     !	*********************CALCULATE THE BONDED FORCE AND POTENTIAL******************************
!     if(.not. nobond) CALL FPBOND (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)
     INCLUDE 'FPBOND.inc'
     !	*******************************************************************************************
     !	*************************CALCULATE THE ANGLE FORCE AND POTENTIAL***************************
!     CALL FPANGLE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)
     INCLUDE 'FPANGLE.inc'
     !	*******************************************************************************************
     !	**********************CALCULATE THE TORSION FORCE AND POTENTIAL****************************
!     CALL FPTOR (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)
     INCLUDE 'FPTOR.inc'
     !	*******************************************************************************************
     !	**********************CALCULATE THE IMPROPER TORSION FORCE AND POTENTIAL********************
     CALL FPOUTPLANE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)
!     INCLUDE 'FPOUTPLANE.inc'

     FXL(I) = FXL(I) + FXI
     FYL(I) = FYL(I) + FYI
     FZL(I) = FZL(I) + FZI

  END DO !300        CONTINUE
!$OMP END DO

!Collate bonded forces
DO I=1,NATOMS
   FX(I) = FX(I) + FXL(I) 
   FY(I) = FY(I) + FYL(I)
   FZ(I) = FZ(I) + FZL(I)
END DO

!$OMP END PARALLEL

  DO I = 1, NATOMS
     PT11 = PT11 + (FX(I) - FXNB(I))*SX(I)
     PT22 = PT22 + (FY(I) - FYNB(I))*SY(I)
     PT33 = PT33 + (FZ(I) - FZNB(I))*SZ(I)
     PT12 = PT12 + (FY(I) - FYNB(I))*SX(I)
     PT13 = PT13 + (FZ(I) - FZNB(I))*SX(I)
     PT23 = PT23 + (FZ(I) - FZNB(I))*SY(I)
  END DO

  RETURN
END SUBROUTINE FORCE

        !	*********************************************************************************************

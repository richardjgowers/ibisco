SUBROUTINE FORCE_HYBR ( )

USE VAR
USE OMP_LIB
IMPLICIT NONE

INTEGER :: I, J, K, M, L, hh, h,A,me,kt,vsite
INTEGER :: JBEG, JEND, JNAB, TI, TJ, TL, TIJ, TIJKL,  NI, TK, TIJK
REAL(kind=rkind) :: RCUTSQ, RCUTIBRSQ, ALPHA, FCUT
REAL(kind=rkind) :: RXI, RYI, RZI, FXI, FYI, FZI, FXIJa, FYIJa, FZIJa
REAL(kind=rkind) :: RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
REAL(kind=rkind) :: VIJ, FIJ, RIJ
REAL(kind=rkind), DIMENSION(NATOMS) :: FXL, FYL, FZL
!HFPANGLE
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

 REAL(KIND=RKIND), DIMENSION(NVIRTA) :: FXVL, FYVL, FZVL

!      REAL,PARAMETER :: hrij = 0.4
DO I = 1, NATOMS
           FX(I) = 0.0D0
           FY(I) = 0.0D0
           FZ(I) = 0.0D0
END DO

VNBOND     = 0.0D0
VBOND      = 0.0D0
VANGLE     = 0.0D0
VTOR       = 0.0D0
VOOP       = 0.0D0
VFXNB      = 0.0D0
VFYNB      = 0.0D0
VFZNB      = 0.0D0
VBOND_MIX  = 0.0D0
VANGLE_MIX = 0.0D0
VTOR_MIX   = 0.0D0
VOOP_MIX   = 0.0D0
VNBOND_MIX = 0.0D0
VBOND_CG   = 0.0D0
VANGLE_CG  = 0.0D0
VTOR_CG    = 0.0D0
VOOP_CG    = 0.0D0
VNBOND_CG  = 0.0D0
FXI = 0.0
FYI = 0.0
FZI = 0.0

!$OMP PARALLEL DEFAULT(NONE)&
!$OMP& SHARED(A,num_bead,num_vs,indx_atm,init_numbcomp,vitype,virtual_center,INDEX_AB,POINT,RX,RY,RZ)&
!$OMP& SHARED(masscoeff,NATOMS,NVIRTA,ITYPE,LIST,TYPE_LABEL,INBONDT,BOXXINV,BOXYINV,BOXZINV,BOXX,BOXY,BOXZ)&
!$OMP& SHARED(FCUTA,FCUTB,BINNB,NDATNB,MASS,INVTOTBMASS,RNBOND,NBOND_FORCE,NBOND_POT,timestepcheck)&
!$OMP& SHARED(SX,SY,SZ,NBONDS,JBOND,IBONDT,typeBond,BINB,BINA,NDATB,RCUT,BOND_FORCE,BOND_POT,RBOND)&
!$OMP& SHARED(NIJK,JANGLEIJK,IANGT,R2D,BEND_FORCE,BEND_POT,ANGLE,NIJKL)&
!$OMP& SHARED(KANGLEIJK,KTORIJKL,JTORIJKL,LTORIJKL,ITORT,BINT,NDATT,TOR_POT,TOR_FORCE,ANGLE_TOR)&
!$OMP& PRIVATE(FXI,FYI,FZI,FXL,FYL,FZL,FXVL,FYVL,FZVL)&
!$OMP& PRIVATE(me,TI,RXJK,RYJK,RZJK,RXKL,RYKL,RZKL,RJK,RKL)&
!$OMP& PRIVATE(hh,I,JBEG,JEND,RXI,RYI,RZI,RXM,RYM,RZM,RM,RXN,RYN,RZN,RXS,RYS,RZS)&
!$OMP& PRIVATE(JNAB,J,TJ,TK,TIJ,TIJK,RXIJ,RYIJ,RZIJ,RIJSQ,TL,TIJKL)&
!$OMP& PRIVATE(RIJ,NI,ALPHA,FIJ,VIJ,VIJK,FIJK,FXIJ,FYIJ,FZIJ,FXIJK,FYIJK,FZIJK)&
!$OMP& PRIVATE(COSM,COTANM,RN,COSN,COTANN,SINM,SINN,COST,SIGNT,PHI,VIJKL)&
!$OMP& PRIVATE(FX1,FY1,FZ1,FX4,FY4,FZ4,FX12,FY12,FZ12)&
!$OMP& PRIVATE(H,vsite,K,FXIJa,FYIJa,FZIJa,FIJKL)&
!$OMP& PRIVATE(FXNB,FYNB,FZNB,RXKJ,RYKJ,RZKJ,RKJ,CT,THETA)&
!$OMP& REDUCTION(+:VNBOND_CG,VNBOND_MIX,VNBOND)&
!$OMP& REDUCTION(+:VBOND,VBOND_MIX,VBOND_CG)&
!$OMP& REDUCTION(+:VANGLE,VANGLE_MIX,VANGLE_CG)&
!$OMP& REDUCTION(+:VTOR,VTOR_MIX,VTOR_CG)&
!$OMP& REDUCTION(+:PT11,PT22,PT33,PT12,PT13,PT23)&
!$OMP& REDUCTION(+:FX,FY,FZ)



!me = omp_get_thread_num()

DO A=1,NATOMS
   FXL(A) = 0.0
   FYL(A) = 0.0
   FZL(A) = 0.0
END DO

DO A=1,NVIRTA
   FXVL(A) = 0.0
   FYVL(A) = 0.0
   FZVL(A) = 0.0
END DO

!###############################################################################
! Loop through BEADS
!$OMP DO SCHEDULE(STATIC,1)
DO hh=1,num_bead !DO 200
   I = INDEX_AB(hh)
   JBEG = POINT(I)
   JEND = POINT(I+1) - 1

   !       ** CHECK THAT BEAD I HAS NEIGHBOURS **
   IF( JBEG .LE. JEND ) THEN
      RXI = RX(I)
      RYI = RY(I)
      RZI = RZ(I)
      !    FXI = FX(I)
      !    FYI = FY(I)
      !    FZI = FZ(I)
      FXI=0.0
      FYI=0.0
      FZI=0.0
      TI = ITYPE(I)

      DO JNAB = JBEG, JEND !Do 199
         !TAKE THE INDEX OF NEIGHBOUR ATOMS
         J = LIST(JNAB)
         if(type_label(j) .eq. 2)then
            !TAKE THE TYPE OF NEIGHBOUR ATOM AND NON-BONDED INTERACTIONS	
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
               IF ( RIJSQ < FCUTB ) THEN
                  RIJ = DSQRT(RIJSQ)
                  NI = INT (RIJ / BINNB(TIJ))
                  IF(NI .GT. NDATNB(TIJ)) THEN
                     WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                     WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ
                     WRITE(1,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                     !WRITE(*,*)'Simulation stopped at time step: ', timestepcheck
                     STOP
                  END IF

                  !		LINEAR INTEPOLATION
                  ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)
                  FIJ = NBOND_FORCE(TIJ,NI)*(1.0-ALPHA) &
                       + ALPHA*NBOND_FORCE(TIJ,NI+1) 
                  VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
                       + ALPHA*NBOND_POT(TIJ,NI+1)

                  VNBOND_CG = VNBOND_CG + VIJ

                  FXIJ  = FIJ * RXIJ
                  FYIJ  = FIJ * RYIJ
                  FZIJ  = FIJ * RZIJ

                  !     QUI HO TOLTO LA PARTE DOVE C'ERA IL DPD, TANTO PER ORA NON SERVE

                  FXI   = FXI + FXIJ
                  FYI   = FYI + FYIJ
                  FZI   = FZI + FZIJ
                  !FX(J) = FX(J) - FXIJ
                  !FY(J) = FY(J) - FYIJ
                  !FZ(J) = FZ(J) - FZIJ
                  FXL(J) = FXL(J) - FXIJ
                  FYL(J) = FYL(J) - FYIJ
                  FZL(J) = FZL(J) - FZIJ

                  !ADD THE NON-BONDED PART OF PRESSURE
                  PT11 = PT11 + FXIJ * RXIJ
                  PT22 = PT22 + FYIJ * RYIJ
                  PT33 = PT33 + FZIJ * RZIJ
                  PT12 = PT12 + FYIJ * RXIJ
                  PT13 = PT13 + FZIJ * RXIJ
                  PT23 = PT23 + FZIJ * RYIJ

               END IF   ! endif  IF ( RIJSQ < RCUTSQ )
            END IF     ! endif   IF( TIJ .NE. 0) THEN
         else !Else if typelabel
            vsite = virtual_center(J)
            TIJ = INBONDT(TI, vitype(vsite))
            IF( TIJ .NE. 0) THEN	
               RXIJ = RXI - RX(J)
               RYIJ = RYI - RY(J)
               RZIJ = RZI - RZ(J)
               RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
               RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
               RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ
               RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0
               IF ( RIJSQ < FCUTB ) THEN
                  RIJ = DSQRT(RIJSQ)
                  NI = INT (RIJ / BINNB(TIJ))
                  IF(NI .GT. NDATNB(TIJ)) THEN
                     WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                     WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ
                     WRITE(1,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                     WRITE(*,*)'Simulation stopped at time step: ', timestepcheck
                     STOP
                  END IF
                  !		LINEAR INTEPOLATION
                  ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)
                  FIJ = NBOND_FORCE(TIJ,NI)*(1.0-ALPHA) &
                       + ALPHA*NBOND_FORCE(TIJ,NI+1) 
                  VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
                       + ALPHA*NBOND_POT(TIJ,NI+1)

                  VNBOND_MIX = VNBOND_MIX + VIJ

                  FXIJ  = FIJ * RXIJ
                  FYIJ  = FIJ * RYIJ
                  FZIJ  = FIJ * RZIJ

                  FXI   = FXI + FXIJ
                  FYI   = FYI + FYIJ
                  FZI   = FZI + FZIJ
                  
!                  do H=1,init_numbcomp(vsite)
!                     K = indx_atm(vsite,H)
!                     FXL(K) = FXL(K) - FXIJ*masscoeff(vsite,H)
!                     FYL(K) = FYL(K) - FYIJ*masscoeff(vsite,H)
!                     FZL(K) = FZL(K) - FZIJ*masscoeff(vsite,H)                    
!                  end do
                  FXVL(vsite) = FXVL(vsite) - FXIJ
                  FYVL(vsite) = FYVL(vsite) - FYIJ
                  FZVL(vsite) = FZVL(vsite) - FZIJ

                  !		ADD THE NON-BONDED PART OF PRESSURE
                  PT11 = PT11 + FXIJ * RXIJ
                  PT22 = PT22 + FYIJ * RYIJ
                  PT33 = PT33 + FZIJ * RZIJ
                  PT12 = PT12 + FYIJ * RXIJ
                  PT13 = PT13 + FZIJ * RXIJ
                  PT23 = PT23 + FZIJ * RYIJ

               END IF   ! endif  IF ( RIJSQ < RCUTB )
            END IF     ! endif   IF( TIJ .NE. 0) THEN
         end if ! if(type_label(j) .eq. 1)
      END DO
      !FX(I) = FX(I) + FXI
      !FY(I) = FY(I) + FYI
      !FZ(I) = FZ(I) + FZI
      FXL(I) = FXL(I) + FXI
      FYL(I) = FYL(I) + FYI
      FZL(I) = FZL(I) + FZI

   ENDIF!If JEND > JBEG
END DO !DO 200
!$OMP END DO

!###############################################################################
!###############################################################################
! Loop through the VS

!$OMP DO SCHEDULE(STATIC,1)
do hh=num_bead+1,num_vs !do 220
   I = INDEX_AB(hh)
   JBEG = POINT(I)
   JEND = POINT(I+1) - 1

   !       ** CHECK THAT BEAD I HAS NEIGHBOURS **

   IF( JBEG .LE. JEND ) THEN

      RXI = RX(I)
      RYI = RY(I)
      RZI = RZ(I)
      !    FXI = FX(I)
      !    FYI = FY(I)
      !    FZI = FZ(I)
      FXI = 0.0
      FYI = 0.0
      FZI = 0.0
      TI = ITYPE(I)

      DO JNAB = JBEG, JEND !Do 219

         !		TAKE THE INDEX OF NEIGHBOUR ATOMS or BEADS
         J = LIST(JNAB)
         if(type_label(j) .eq. 1)then !If J is an atom

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
               IF ( RIJSQ < FCUTA ) THEN
                  RIJ = DSQRT(RIJSQ)
                  NI = INT (RIJ / BINNB(TIJ))
                  IF(NI .GT. NDATNB(TIJ)) THEN
                     WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                     WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ
                     WRITE(1,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                     !WRITE(*,*)'Simulation stopped at time step: ', timestepcheck
                     STOP
                  END IF

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

                  !     QUI HO TOLTO LA PARTE DOVE C'ERA IL DPD, TANTO PER ORA NON SERVE
                  FXI   = FXI + FXIJ
                  FYI   = FYI + FYIJ
                  FZI   = FZI + FZIJ

                  !FX(J) = FX(J) - FXIJ
                  !FY(J) = FY(J) - FYIJ
                  !FZ(J) = FZ(J) - FZIJ
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

               END IF   ! endif  IF ( RIJSQ < RCUTSQ )
            END IF     ! endif   IF( TIJ .NE. 0) THEN

         else ! Here we have the interaction with beads

            TJ = ITYPE(J)
            vsite=virtual_center(I)
            TIJ = INBONDT(vitype(vsite), TJ)

            IF( TIJ .NE. 0) THEN	
               RXIJ = RXI - RX(J)
               RYIJ = RYI - RY(J)
               RZIJ = RZI - RZ(J)
               RXIJ = RXIJ - ANINT ( RXIJ * BOXXINV ) * BOXX
               RYIJ = RYIJ - ANINT ( RYIJ * BOXYINV ) * BOXY
               RZIJ = RZIJ - ANINT ( RZIJ * BOXZINV ) * BOXZ
               RIJSQ = RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0

               IF ( RIJSQ < FCUTB ) THEN
                  RIJ = DSQRT(RIJSQ)
                  NI = INT (RIJ / BINNB(TIJ))
                  IF(NI .GT. NDATNB(TIJ)) THEN
                     WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                     WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ
                     WRITE(1,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                     WRITE(*,*)'Simulation stopped at time step: ', timestepcheck
                     STOP
                  END IF
                  !		LINEAR INTEPOLATION
                  ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)
                  FIJ = NBOND_FORCE(TIJ,NI)*(1.0-ALPHA) &
                       + ALPHA*NBOND_FORCE(TIJ,NI+1) 
                  VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
                       + ALPHA*NBOND_POT(TIJ,NI+1)

                  VNBOND_MIX = VNBOND_MIX + VIJ

                  FXIJ  = FIJ * RXIJ
                  FYIJ  = FIJ * RYIJ
                  FZIJ  = FIJ * RZIJ
 
                  FXVL(vsite) = FXVL(vsite) + FXIJ
                  FYVL(vsite) = FYVL(vsite) + FYIJ
                  FZVL(vsite) = FZVL(vsite) + FZIJ
               
!                 do H=1,init_numbcomp(vsite)
!                    K = indx_atm(vsite,H)
!                    FXL(K) = FXL(K) + FXIJ*masscoeff(vsite,H)
!                    FYL(K) = FYL(K) + FYIJ*masscoeff(vsite,H)
!                    FZL(K) = FZL(K) + FZIJ*masscoeff(vsite,H)                 
!                 end do

                  !FX(J) = FX(J) - FXIJ
                  !FY(J) = FY(J) - FYIJ
                  !FZ(J) = FZ(J) - FZIJ
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

               END IF   ! endif  IF ( RIJSQ < RCUTSQ )
            END IF     ! endif   IF( TIJ .NE. 0) THEN
         end if ! if(type_label(j) .eq. 1)
      END DO !do 219

      !FX(I) = FX(I) + FXI
      !FY(I) = FY(I) + FYI
      !FZ(I) = FZ(I) + FZI
      FXL(I) = FXL(I) + FXI
      FYL(I) = FYL(I) + FYI
      FZL(I) = FZL(I) + FZI
   ENDIF
end do !do 220
!$OMP END DO

DO A=1,NVIRTA
   DO H=1,init_numbcomp(A)
      K = indx_atm(A,H)
      FXL(K) = FXL(K) + FXVL(A)*masscoeff(A,H)
      FYL(K) = FYL(K) + FYVL(A)*masscoeff(A,H)
      FZL(K) = FZL(K) + FZVL(A)*masscoeff(A,H)
   END DO
END DO

!###############################################################################
!###############################################################################
! Loop through atoms

!$OMP DO SCHEDULE(STATIC,1)
do hh=num_vs+1,natoms !Do 210
   i = INDEX_AB(hh)
   JBEG = POINT(I)
   JEND = POINT(I+1) - 1

   !       ** CHECK THAT BEAD I HAS NEIGHBOURS **
   IF( JBEG .LE. JEND ) THEN
      RXI = RX(I)
      RYI = RY(I)
      RZI = RZ(I)
      !    FXI = FX(I)
      !    FYI = FY(I)
      !    FZI = FZ(I)
      FXI = 0.0
      FYI = 0.0
      FZI = 0.0
      TI = ITYPE(I)

      DO JNAB = JBEG, JEND !DO 209

         !TAKE THE INDEX OF NEIGHBOUR ATOMS
         J = LIST(JNAB)
         !TAKE THE TYPE OF NEIGHBOUR ATOM AND NON-BONDED INTERACTIONS	
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
            IF ( RIJSQ < FCUTA ) THEN
               RIJ = DSQRT(RIJSQ)
               NI = INT (RIJ / BINNB(TIJ))
               IF(NI .GT. NDATNB(TIJ)) THEN
                  WRITE(*,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                  WRITE(*,*)'Atom/Bead ',I,J,'RIJ ',RIJ
                  WRITE(1,*)'FATAL ERROR: Entry in non bonded table', TIJ,' does not exist'
                  !WRITE(*,*)'Simulation stopped at time step: ', timestepcheck
                  STOP
               END IF

               !LINEAR INTEPOLATION

               ALPHA=(RIJ-RNBOND(TIJ,NI))/BINNB(TIJ)
               FIJ = NBOND_FORCE(TIJ,NI)*(1.0-ALPHA) &
                    + ALPHA*NBOND_FORCE(TIJ,NI+1) 
               VIJ = NBOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
                    + ALPHA*NBOND_POT(TIJ,NI+1)
               VNBOND = VNBOND + VIJ
               FXIJ  = FIJ * RXIJ
               FYIJ  = FIJ * RYIJ
               FZIJ  = FIJ * RZIJ

               !QUI HO TOLTO LA PARTE DOVE C'ERA IL DPD, TANTO PER ORA NON SERVE

               FXI   = FXI + FXIJ
               FYI   = FYI + FYIJ
               FZI   = FZI + FZIJ
               !FX(J) = FX(J) - FXIJ
               !FY(J) = FY(J) - FYIJ
               !FZ(J) = FZ(J) - FZIJ
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

            END IF   ! endif  IF ( RIJSQ < RCUTSQ )
         END IF     ! endif   IF( TIJ .NE. 0) THEN
      END DO !DO 209

      !FX(I) = FX(I) + FXI
      !FY(I) = FY(I) + FYI
      !FZ(I) = FZ(I) + FZI
      FXL(I) = FXL(I) + FXI
      FYL(I) = FYL(I) + FYI
      FZL(I) = FZL(I) + FZI

   ENDIF

END DO !DO 210
!$OMP END DO

do A=1,NATOMS
   FX(A) = FX(A) + FXL(A)
   FY(A) = FY(A) + FYL(A)
   FZ(A) = FZ(A) + FZL(A)
END DO


!$OMP END PARALLEL

!###############################################################################
!###############################################################################

!	SAVE THE NON-BONDED PART OF FORCE
	DO I=1, NATOMS
	
	FXNB(I) = FX(I) + VFXNB(I)
	FYNB(I) = FY(I) + VFYNB(I)
	FZNB(I) = FZ(I) + VFZNB(I)

	END DO	

!$OMP PARALLEL DEFAULT(NONE)&
!$OMP& SHARED(A,num_bead,num_vs,indx_atm,init_numbcomp,vitype,virtual_center,INDEX_AB,POINT,RX,RY,RZ)&
!$OMP& SHARED(masscoeff,NATOMS,ITYPE,LIST,TYPE_LABEL,INBONDT,BOXXINV,BOXYINV,BOXZINV,BOXX,BOXY,BOXZ)&
!$OMP& SHARED(FCUTA,FCUTB,BINNB,NDATNB,MASS,INVTOTBMASS,RNBOND,NBOND_FORCE,NBOND_POT,timestepcheck)&
!$OMP& SHARED(SX,SY,SZ,NBONDS,JBOND,IBONDT,typeBond,BINB,BINA,NDATB,RCUT,BOND_FORCE,BOND_POT,RBOND)&
!$OMP& SHARED(NIJK,JANGLEIJK,IANGT,R2D,BEND_FORCE,BEND_POT,ANGLE,NIJKL)&
!$OMP& SHARED(KANGLEIJK,KTORIJKL,JTORIJKL,LTORIJKL,ITORT,BINT,NDATT,TOR_POT,TOR_FORCE,ANGLE_TOR)&
!$OMP& PRIVATE(FXI,FYI,FZI,FXL,FYL,FZL,me,TI,RXJK,RYJK,RZJK,RXKL,RYKL,RZKL,RJK,RKL)&
!$OMP& PRIVATE(hh,I,JBEG,JEND,RXI,RYI,RZI,RXM,RYM,RZM,RM,RXN,RYN,RZN,RXS,RYS,RZS)&
!$OMP& PRIVATE(JNAB,J,TJ,TK,TIJ,TIJK,RXIJ,RYIJ,RZIJ,RIJSQ,TL,TIJKL)&
!$OMP& PRIVATE(RIJ,NI,ALPHA,FIJ,VIJ,VIJK,FIJK,FXIJ,FYIJ,FZIJ,FXIJK,FYIJK,FZIJK)&
!$OMP& PRIVATE(COSM,COTANM,RN,COSN,COTANN,SINM,SINN,COST,SIGNT,PHI,VIJKL)&
!$OMP& PRIVATE(FX1,FY1,FZ1,FX4,FY4,FZ4,FX12,FY12,FZ12)&
!$OMP& PRIVATE(H,vsite,K,FXIJa,FYIJa,FZIJa,FIJKL)&
!$OMP& PRIVATE(FXNB,FYNB,FZNB,RXKJ,RYKJ,RZKJ,RKJ,CT,THETA)&
!$OMP& REDUCTION(+:VNBOND_CG,VNBOND_MIX,VNBOND)&
!$OMP& REDUCTION(+:VBOND,VBOND_MIX,VBOND_CG)&
!$OMP& REDUCTION(+:VANGLE,VANGLE_MIX,VANGLE_CG)&
!$OMP& REDUCTION(+:VTOR,VTOR_MIX,VTOR_CG)&
!$OMP& REDUCTION(+:PT11,PT22,PT33,PT12,PT13,PT23)&
!$OMP& REDUCTION(+:FX,FY,FZ)

DO I=1,NATOMS
   FXL(I) = 0.0
   FYL(I) = 0.0
   FZL(I) = 0.0
END DO


!$OMP DO SCHEDULE(STATIC,1)
      DO  I = 1, NATOMS

      RXI = SX(I)
      RYI = SY(I)
      RZI = SZ(I)
      !FXI = FX(I)
      !FYI = FY(I)
      !FZI = FZ(I)
       FXI = 0.0 
       FYI = 0.0
       FZI = 0.0 
	
      TI = ITYPE(I)

!	*******************************************************************************************
!	*********************CALCULATE THE BONDED FORCE AND POTENTIAL******************************
!      CALL HFPBOND (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)
      INCLUDE 'HFPBOND.inc' 

!	*******************************************************************************************
!	*************************CALCULATE THE ANGLE FORCE AND POTENTIAL***************************
!      CALL HFPANGLE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)
      INCLUDE 'HFPANGLE.inc' 

!	*******************************************************************************************
!	**********************CALCULATE THE TORSION FORCE AND POTENTIAL****************************
!      CALL HFPTOR (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)
       INCLUDE 'HFPTOR.inc' 

!	*******************************************************************************************
!	**********************CALCULATE THE IMPROPER TORSION FORCE AND POTENTIAL********************
!     CALL FPOUTPLANE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

      FXL(I) = FXL(I) + FXI
      FYL(I) = FYL(I) + FYI
      FZL(I) = FZL(I) + FZI
End do 

!$OMP END DO 

!Collate forces again
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
      END SUBROUTINE FORCE_HYBR

!	*********************************************************************************************

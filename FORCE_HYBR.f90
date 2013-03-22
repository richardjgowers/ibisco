
SUBROUTINE FORCE_HYBR ( )

USE VAR
USE OMP_LIB
IMPLICIT NONE

INTEGER :: I, J, hh, h,A,me,k,kt,vsite
INTEGER :: JBEG, JEND, JNAB, TI, TJ, TIJ, NI
REAL(kind=rkind) :: RCUTSQ, RCUTIBRSQ, ALPHA, FCUT
REAL(kind=rkind) :: RXI, RYI, RZI, FXI, FYI, FZI, FXIJa, FYIJa, FZIJa
REAL(kind=rkind) :: RXIJ, RYIJ, RZIJ, RIJSQ, FXIJ, FYIJ, FZIJ
REAL(kind=rkind) :: VIJ, FIJ, RIJ
REAL(kind=rkind), DIMENSION(NATOMS) :: FXL, FYL, FZL

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
!$OMP& SHARED(masscoeff,NATOMS,ITYPE,LIST,TYPE_LABEL,INBONDT,BOXXINV,BOXYINV,BOXZINV,BOXX,BOXY,BOXZ)&
!$OMP& SHARED(FCUTA,FCUTB,BINNB,NDATNB,MASS,INVTOTBMASS,RNBOND,NBOND_FORCE,NBOND_POT,timestepcheck)&
!$OMP& PRIVATE(FXI,FYI,FZI,FXL,FYL,FZL,me)&
!$OMP& PRIVATE(hh,I,JBEG,JEND,RXI,RYI,RZI,TI)&
!$OMP& PRIVATE(JNAB,J,TJ,TIJ,RXIJ,RYIJ,RZIJ,RIJSQ)&
!$OMP& PRIVATE(RIJ,NI,ALPHA,FIJ,VIJ,FXIJ,FYIJ,FZIJ)&
!$OMP& PRIVATE(H,vsite,k,FXIJa,FYIJa,FZIJa)&
!$OMP& REDUCTION(+:VNBOND_CG,VNBOND_MIX,VNBOND)&
!$OMP& REDUCTION(+:PT11,PT22,PT33,PT12,PT13,PT23)&
!$OMP& REDUCTION(+:FX,FY,FZ)

!me = omp_get_thread_num()

DO A=1,NATOMS
   FXL(A) = 0.0
   FYL(A) = 0.0
   FZL(A) = 0.0
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
                  
                  do H=1,init_numbcomp(vsite)
                     K = indx_atm(vsite,H)
                     FXL(K) = FXL(K) - FXIJ*masscoeff(vsite,H)
                     FYL(K) = FYL(K) - FYIJ*masscoeff(vsite,H)
                     FZL(K) = FZL(K) - FZIJ*masscoeff(vsite,H)                    
                  end do
                  
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
         if(type_label(j) .eq. 1)then !If J is a bead

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
                
                  do H=1,init_numbcomp(vsite)
                     K = indx_atm(vsite,H)
                     FXL(K) = FXL(K) + FXIJ*masscoeff(vsite,H)
                     FYL(K) = FYL(K) + FYIJ*masscoeff(vsite,H)
                     FZL(K) = FZL(K) + FZIJ*masscoeff(vsite,H)                 
                  end do

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
!   !$OMP ATOMIC
   FX(A) = FX(A) + FXL(A)
!   !$OMP ATOMIC
   FY(A) = FY(A) + FYL(A)
!   !$OMP ATOMIC
   FZ(A) = FZ(A) + FZL(A)
END DO

!$OMP END PARALLEL

!###############################################################################
!###############################################################################

!	SAVE THE NON-BONDED PART OF FORCE
	DO I=1, NATOMS
	
	FXNB(I) = FX(I) !+ VFXNB(I)
	FYNB(I) = FY(I) !+ VFYNB(I)
	FZNB(I) = FZ(I) !+ VFZNB(I)

	END DO	

      DO 300 I = 1, NATOMS

      RXI = SX(I)
      RYI = SY(I)
      RZI = SZ(I)
      FXI = FX(I)
      FYI = FY(I)
      FZI = FZ(I)
	
      TI = ITYPE(I)

!	*******************************************************************************************
!	*********************CALCULATE THE BONDED FORCE AND POTENTIAL******************************
      CALL HFPBOND (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

!	*******************************************************************************************
!	*************************CALCULATE THE ANGLE FORCE AND POTENTIAL***************************
      CALL HFPANGLE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

!	*******************************************************************************************
!	**********************CALCULATE THE TORSION FORCE AND POTENTIAL****************************
      CALL HFPTOR (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

!	*******************************************************************************************
!	**********************CALCULATE THE IMPROPER TORSION FORCE AND POTENTIAL********************
     CALL FPOUTPLANE (I, RXI, RYI, RZI, FXI, FYI, FZI, TI)

      FX(I) = FXI
      FY(I) = FYI
      FZ(I) = FZI

300	CONTINUE

      DO I = 1, NATOMS

            PT11 = PT11 + (FX(I) - FXNB(I))*SX(I)
            PT22 = PT22 + (FY(I) - FYNB(I))*SY(I)
            PT33 = PT33 + (FZ(I) - FZNB(I))*SZ(I)

            PT12 = PT12 + (FY(I) - FYNB(I))*SX(I)
            PT13 = PT13 + (FZ(I) - FZNB(I))*SX(I)
            PT23 = PT23 + (FZ(I) - FZNB(I))*SY(I)

      END DO


      RETURN
      END

!	*********************************************************************************************

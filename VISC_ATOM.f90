!***************************************************************************************************
!            SUBROUTINE VISC_ATOM (TIMESTEP, BOXZ, NUMA, VX, ZN, TRANSFER)
             
             SUBROUTINE VISC_ATOM (TIMESTEP)
  
              USE  RNEMD
              USE VAR
              IMPLICIT NONE
              INTEGER  TIMESTEP
              INTEGER*4  NFLAG, I, J, K, IFIRST 
              REAL*8 POSINF1, POSSUP1,    &
          	     POSINF2, POSSUP2, VX_SLAB1, VX_SLAB_MIDL, &
           	     OLD_VX_SLABMIDL,  OLD_VX_SLAB1,  &
           	     CURRENT_V , AMJ, AMK
  

!          -VX        +VX
!          _____________________
!          |S| | | | |S| | | | |
!          |L| | | | |L| | | | |
!          |A| | | | |A| | | | |
!          |B| | | | |B| | | | |
!          | | | | | | | | | | |
!          |1| | | | |M| | | | |
!
!      -boxz/2________0_________+boxz/2 ------------> Z direction
!          ^ ^       ^ ^
!         /   \     /   \
!        /     |   |     \
!       |      |   |      |
!    POSINF1   |   |   POSSUP2
!         POSSUP1  |
!                 POSINF2
!     DEFINITION OF THE THICKNESS OF SLABS

      SLAB_THICKNESS = BOXZ/NUMSLAB

      POSINF1 =  -BOXZ/2.0D0 
      POSSUP1 =   SLAB_THICKNESS - BOXZ/2.0D0

      POSINF2 =  (NUMSLAB/2.0D00 ) * SLAB_THICKNESS - BOXZ/2.0D00      ! this values should equal to zero
      POSSUP2 =  ( NUMSLAB/2.0D00 + 1.0D00 ) * SLAB_THICKNESS - BOXZ/2.0D0
      
!     LOOP, WHICH FINDS THE TWO ATOMS WHICH WILL HAVE VELOCITIES EXCHANGED
!     LOOP FOR ATOM IN SLAB 1
     
      NFLAG = 0
        DO I = 1, NATOMS
          
          IF ((RZ(I).GE.POSINF1).AND.((RZ(I).LT.POSSUP1))) THEN
              IF(NFLAG.EQ.0) THEN
               VX_SLAB1 = VX(I)
               J=I
               NFLAG=1
              ELSE
                CURRENT_V =  VX(I)
                IF (CURRENT_V .LT. VX_SLAB1) THEN
                  VX_SLAB1 = CURRENT_V
                  J = I
                ENDIF
              ENDIF
           ENDIF
         ENDDO
         
!     SECOND ATOM, MIDDLE SLAB -  1+N/2

      NFLAG = 0
      K = 0                           
      DO I = 1, NATOMS

         IF ((RZ(I).GE.POSINF2).AND.((RZ(I).LT.POSSUP2))) THEN
            IF (AM(I) .EQ. AM(J)) THEN
               IF(NFLAG.EQ.0) THEN
                  VX_SLAB_MIDL = VX(I)
                  K=I
                  NFLAG = 1
               ELSE
                 CURRENT_V = VX(I)
                 IF (CURRENT_V .GT. VX_SLAB_MIDL) THEN
                     VX_SLAB_MIDL = CURRENT_V
                     K=I
                 ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
         
 
  !   MOMENTUM TRANSFERRED
      IF (K.GT.0) THEN                     !hjqian: ONLY SWAP WHEN DO HAVE ATOM IN MIDDLE SLAB
      AMJ = AM(J) * MASSSCALE
      AMK = AM(K) * MASSSCALE
      TRANSFER = - AMJ * VX_SLAB1*VSCALE + AMK * VX_SLAB_MIDL*VSCALE   ! UNIT: KG * M/S
      OLD_VX_SLAB1 = VX(J)
      OLD_VX_SLABMIDL = VX(K)
      VX(J) = OLD_VX_SLABMIDL
      VX(K) = OLD_VX_SLAB1
!     WRITE(555,*) J,K,VX(J),VX(K),VX_SLAB1,VX_SLAB_MIDL
      ENDIF
      END SUBROUTINE

!***************************************************************************

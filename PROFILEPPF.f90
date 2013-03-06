       subroutine PROFILEPPF ()

       use VAR 

       IMPLICIT NONE
       
       integer IZ,J

       REAL*8 CMASS,RZI,PVX

       
      DO IZ = 1,SLIDE

       NUM(IZ) = 0
       DENSITY(IZ) = 0.d0
       TEEMP(IZ) = 0.d0

       VXPRO(IZ) = 0.d0
      ENDDO

      DEN1 = 0.D0
      DEN2 = 0.D0
      DO J = 1, NATOMS

         CMASS = BEADMASS(J)
         RZI = SZ(J)
         RZI = RZI - BOXZ*ANINT(RZI*BOXZINV)
         if (RZI.GE.0) then
         DEN2 = DEN2 + CMASS
         elseif(RZI.LT.0) then
         DEN1 = DEN1 + CMASS
         endif

         IZ = int((RZI*BOXZINV + 0.5)*real(SLIDE))
         IZ = 1+MOD((IZ+SLIDE),SLIDE)
         Z_POSITION(J) = IZ

         NUM(IZ) = NUM(IZ) + 1

         VXPRO(IZ) = VXPRO(IZ) + VTX(J)
         DENSITY(IZ) = DENSITY(IZ) + CMASS
      ENDDO
      DO J = 1,SLIDE
         VXPRO(J) = VXPRO(J)/NUM(J)
      ENDDO
      DO  J = 1,NATOMS
         CMASS = BEADMASS(J)
!        RZI = SZ(J)
!        RZI = RZI - BOXZ*ANINT(RZI*BOXZINV)
!        IZ = int((RZI*BOXZINV + 0.5)*real(SLIDE))
!        IZ = 1+MOD((IZ+SLIDE),SLIDE)
         IZ = Z_POSITION(J)
         VPX(J) = VTX(J) - VXPRO(IZ)
         PVX = VPX(J)
         TEEMP(IZ) = TEEMP(IZ) +              &
                     CMASS*(PVX**2 +          &
                         VTY(J)**2 +          &
                         VTZ(J)**2 )
      ENDDO

      RETURN
      END


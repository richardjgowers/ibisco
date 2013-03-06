       subroutine PROFILERNEMD ()

       use VAR 
       use RNEMD
       IMPLICIT NONE
       
       integer IZ,J

       REAL*8 CMASS,RZI

       
      DO IZ = 1,NUMSLAB
       VX_SLAB(IZ) = 0.d0
       NUMA_SLAB(IZ) = 0
      ENDDO

      DO J = 1, NATOMS

         CMASS = BEADMASS(J)
         RZI = SZ(J)
         RZI = RZI - BOXZ*ANINT(RZI*BOXZINV)

         IZ = int((RZI*BOXZINV + 0.5)*real(NUMSLAB))
         IZ = 1+MOD((IZ+NUMSLAB),NUMSLAB)
         Z_POSITION(J) = IZ

         NUMA_SLAB(IZ) = NUMA_SLAB(IZ) + 1

         VX_SLAB(IZ) = VX_SLAB(IZ) + VTX(J)
      ENDDO
      DO J = 1,NUMSLAB
         VX_SLAB(J) = VX_SLAB(J)/NUMA_SLAB(J)
      ENDDO
      DO  J = 1,NATOMS
         CMASS = BEADMASS(J)
!        RZI = SZ(J)
!        RZI = RZI - BOXZ*ANINT(RZI*BOXZINV)
!        IZ = int((RZI*BOXZINV + 0.5)*real(NUMSLAB))
!        IZ = 1+MOD((IZ+NUMSLAB),NUMSLAB)
         IZ = Z_POSITION(J)
         VPX(J) = VTX(J) - VX_SLAB(IZ)
      ENDDO

      RETURN
      END


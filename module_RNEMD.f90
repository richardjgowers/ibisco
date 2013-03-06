
      MODULE RNEMD

      INTEGER    VISCINPUT, NUMPRO,NSWP
      INTEGER    NEXCH,NEMDPROF,NEMDTRAJ,NUMSLAB
!       INTEGER, POINTER :: Z_POSITION(:)
      REAL*8   SLAB_THICKNESS , TRANSFER,TTRANSF
      REAL*8, POINTER :: AM(:)    ! Mass array for particle
      REAL*8, POINTER :: VXMEAN_SLAB(:),VX_SLAB(:),NUMA_SLAB(:), &
                              SLAB_TEMP(:),SLAB_DENS(:),&
                              SLAB_MTEMP(:),SLAB_MDENS(:)
          
      END MODULE RNEMD


!***************************************************************************

       SUBROUTINE VISC_PROFILE (TIMESTEP)

       USE RNEMD
       USE VAR
         
      IMPLICIT NONE
      INTEGER  TIMESTEP 
      INTEGER  I, J, IFIRST,IZ
      REAL*8   MM            ! COORDINATE OF CENTRE OF  GIVEN SLAB
      REAL*8  RZI
      REAL*8  SLAB_VOL,SLAB_TEM,TIMEINV
      REAL*8  BOUND_LO,BOUND_HI
      REAL*8  CUR_GRADIENT1, CUR_GRADIENT2, CUR_GRADIENT
      REAL*8  GRADIENT1,GRADIENT2,GRADIENT
      REAL*8  CUR_FLUX,  CUR_VISCOSITY, TIME,RTIME, AREA,FLUX, MVISCOSITY

      REAL*8  DUMMY
      REAL*8         VOLSCALE,VXP,MKTEMPT
      PARAMETER  ( VOLSCALE = 1.0D-27)
      REAL*8 SLAB_MASS(NUMSLAB)
      REAL*8 ZSLAB(NUMSLAB)
        
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SLAB_VOL =   (BOXX * BOXY * SLAB_THICKNESS)* VOLSCALE 

        DO I = 1, NUMSLAB
           SLAB_MASS(I) = 0.d0
           SLAB_TEMP(I) = 0.D0
        ENDDO
 
        DO I = 1, NATOMS
           IZ = Z_POSITION(I) 
           SLAB_MASS(IZ) = SLAB_MASS(IZ) + AM(I)
        ENDDO
        
        DO I = 1, NUMSLAB
           ZSLAB(I)     = SLAB_THICKNESS*(REAL(I)-0.5) - BOXZ/2.0D0          ! UNIT:NANOMETER
           SLAB_MASS(I) = SLAB_MASS(I)*MASSSCALE
           SLAB_DENS(I) = SLAB_MASS(I)/SLAB_VOL
        ENDDO
        
        DO I = 1, NATOMS
           IZ = Z_POSITION(I)
           VXP = VPX(I)
           SLAB_TEMP(IZ) = SLAB_TEMP(IZ)+0.5d0*AM(I)* &
                           (VXP**2+VY(I)**2+VZ(I)**2)
        ENDDO
        
!cccccccccccAVERAGE OVER ALL THE TIMESTEPScccccccccccccccccccccccccc
        DO IZ = 1, NUMSLAB
           VXMEAN_SLAB(IZ) = VXMEAN_SLAB(IZ) + VX_SLAB(IZ)
           SLAB_MTEMP(IZ) = SLAB_MTEMP(IZ) + SLAB_TEMP(IZ)
           SLAB_MDENS(IZ) = SLAB_MDENS(IZ) + SLAB_DENS(IZ)
        ENDDO
!cccccccccccccccccccc OUTPUT THE AVERAGED PROFILES IN SLAB ccccccccccccccccccccc
        IF (MOD(TIMESTEP,NEMDPROF).EQ.0.and.MOD(TIMESTEP,NEXCH).NE.0) THEN
         TIMEINV = 1.D0/DBLE(NUMPRO)
!cccccccccccAVERAGE OVER ALL THE TIMESTEPScccccccccccccccccccccccccc
        WRITE(202,*) 'STEP: ', TIMESTEP
        DO IZ = 1, NUMSLAB
           MM = ZSLAB(IZ)
           MKTEMPT    = 2.D0/DBLE(3.D0*NUMA_SLAB(IZ)-3.D0)
           SLAB_TEM    = SLAB_MTEMP(IZ)*TEMPSCALE*MKTEMPT
           WRITE(202,'(I10,4(E16.8, 2X))') IZ,MM,VXMEAN_SLAB(IZ)*VSCALE*TIMEINV,SLAB_TEM*TIMEINV,SLAB_MDENS(IZ)*TIMEINV
        ENDDO
        ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       do the linear regression to calcualte the velocity gradient
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF (MOD(TIMESTEP,NEMDTRAJ).EQ.0) THEN
       TIMEINV = 1.D0/DBLE(NUMPRO)
      CALL LINEAR (NUMSLAB/2-1, ZSLAB(2), VX_SLAB(2),   &
                  CUR_GRADIENT1, DUMMY, DUMMY)

      CALL LINEAR (NUMSLAB/2-1, ZSLAB(NUMSLAB/2+2),            &
                  VX_SLAB(NUMSLAB/2+2),                &
                 CUR_GRADIENT2, DUMMY, DUMMY)
      CUR_GRADIENT = 0.5D00 * (CUR_GRADIENT1 - CUR_GRADIENT2)/TIMESCALE      ! UNIT: S**(-1)


      TIME         = TIMESTEP * DT* TIMESCALE       ! TIME*TIMESCALE   =>  CONVERT REDUCED UNIT TO SECOND

      RTIME = NSWP*DT*TIMESCALE*NEXCH

      AREA         =  BOXX*BOXY * ( RSCALE**2)        ! AREA =>  CONVERT REDUCED UNIT TO METER**2

      CUR_FLUX      = TRANSFER / (2.0D00 * AREA * NEXCH * DT* TIMESCALE)   ! DT* TIMESCALE CONVERT REDUCED UNIT TO SECOND
      CUR_VISCOSITY = - 1.0D+3 * CUR_FLUX / CUR_GRADIENT  !UNIT: MPA*S OR CP

      CALL LINEAR (NUMSLAB/2-1, ZSLAB(2), VXMEAN_SLAB(2),   &
                 GRADIENT1, DUMMY, DUMMY)

      CALL LINEAR (NUMSLAB/2-1, ZSLAB(NUMSLAB/2+2),            &
                  VXMEAN_SLAB(NUMSLAB/2+2),                &
                 GRADIENT2, DUMMY, DUMMY)

      GRADIENT = 0.5D00 * (GRADIENT1 - GRADIENT2)*TIMEINV/TIMESCALE
      FLUX         = TTRANSF / (2.0D00 * AREA * RTIME)  ! UNIT: KG/(M* S**2)
      MVISCOSITY     = - 1.0D+3 * FLUX / GRADIENT    !UNIT: MPA*S OR CP

!     ** VISCOSITY IN CP **

        WRITE(201,'(I10,7(1x,E12.5))') TIMESTEP, TIME,  &
          CUR_FLUX,FLUX, CUR_GRADIENT, GRADIENT,   &
          CUR_VISCOSITY,MVISCOSITY
       ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
        RETURN
        END SUBROUTINE  

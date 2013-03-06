       SUBROUTINE VISC_NEMD ( TIMESTEP )
   
       USE RNEMD
       USE VAR
       IMPLICIT NONE
       INTEGER  I, J, TIMESTEP
       REAL*8   ZSLAB(NUMSLAB)
      
      IF (ENSEMBLE == 2) THEN
        SLAB_THICKNESS  =  BOXZ/ REAL( NUMSLAB)
      ENDIF
	
      IF(MOD(TIMESTEP, NEXCH) .EQ. 0) THEN
         
         NSWP = NSWP + 1        
         call VISC_ATOM (TIMESTEP)
         TTRANSF = TTRANSF + TRANSFER
      END IF

!     IF (MOD(TIMESTEP, NEMDPROF) .EQ. 0) THEN
!      IF (MOD(TIMESTEP, NEXCH) .NE. 0) THEN
         NUMPRO = NUMPRO + 1
         CALL VISC_PROFILE (TIMESTEP)
!      ENDIF
!     ENDIF
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
!     IF (TIMESTEP.LE.TSET) THEN
!        TTRANSF = 0.D0
!        DO I = 1,  NUMSLAB
!        VXMEAN_SLAB (I) =  0.0D00
!      ENDDO
!     ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      RETURN
      END SUBROUTINE



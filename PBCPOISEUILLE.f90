      subroutine PBCPOISEUILLE(IT)
      
      USE VAR
      IMPLICIT NONE
      INTEGER  I, J,IZ,IT,tt
      REAL*8   RXI, RYI, RZI,VXI, VYI, VZI,PVX
      REAL*8   V_AVE1,V_AVE2,AVTINV
      REAL*8   VISCOSITY1,VISCOSITY2,RVISCOSITY1,RVISCOSITY2
      REAL*8   VGR(500),RVGR(500)
      REAL*8   TREAL,slidbin,slidvol,slidvolINV,dslidbinINV
      REAL*8   MKTEM
       
!     AVTINV = 1.D0/DBLE(AVT)
      SLIDE2 = SLIDE/2
      slidbin = BOXZ/DBLE(SLIDE)
      dslidbinINV = 0.5D0/slidbin
      slidvol = BOXX*BOXY*slidbin
      slidvolINV = 1.D0/slidvol
      TREAL = IT * DT * TIMESCALE * 1.0D+12 + INITIME !ps
     
      DO IZ = 1,SLIDE
       PZX(IZ) = PZX(IZ)*slidvolINV
      ENDDO

      DO J = 1,SLIDE
         DENSITY(J) = DENSITY(J)*slidvolINV
         MKTEM = 1.D0/(3.D0*(NUM(J)-1))
         TEEMP(J) = TEEMP(J)*MKTEM*TEMPSCALE
      ENDDO       
      
      VGR(1) = VXPRO(2) - VXPRO(1)
      VGR(1) = VGR(1)*dslidbinINV*2
      VVISCOSITY(1) = -1.D0*PZX(1)/VGR(1)      !local viscosity
      do J = 2,SLIDE-1     
         VGR(J) = VXPRO(J+1)-VXPRO(J-1)  !calculate velocity gradient
         VGR(J) = VGR(J)*dslidbinINV
         VVISCOSITY(J) = -1.D0*PZX(J)/VGR(J)
      ENDDO
      VGR(SLIDE) = VXPRO(SLIDE)-VXPRO(SLIDE-1)
      VGR(SLIDE) = VGR(SLIDE)*dslidbinINV*2
      VVISCOSITY(SLIDE)=-1.D0*PZX(SLIDE)/VGR(SLIDE)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!AVERAGEG PROFILES!!!!!!!!!!!!!!!!!
      tt = iabs(mod(IT,AVT))
      IF (IT.LE.TSET) THEN
         AVTINV = 1.D0/DBLE(IT)
      ELSE
         AVTINV = 1.D0/DBLE(IT-TSET)
      ENDIF

         DO J = 1,SLIDE
           RVXPRO(J) = RVXPRO(J)+VXPRO(J)
           RDENSITY(J) = RDENSITY(J) + DENSITY(J)
           RTEEMP(J) = RTEEMP(J) + TEEMP(J)
           RPZX(J) = RPZX(J) + PZX(J)
         ENDDO

      RVGR(1) = RVXPRO(2) - RVXPRO(1)
      RVGR(1) = RVGR(1)*dslidbinINV*2
      RVVISCOSITY(1) = -1.D0*RPZX(1)/RVGR(1)      !local viscosity
      do J = 2,SLIDE-1
         RVGR(J) = RVXPRO(J+1)-RVXPRO(J-1)  !calculate velocity gradient
         RVGR(J) = RVGR(J)*dslidbinINV
         RVVISCOSITY(J) = -1.D0*RPZX(J)/RVGR(J)
      ENDDO
      RVGR(SLIDE) = RVXPRO(SLIDE)-RVXPRO(SLIDE-1)
      RVGR(SLIDE) = RVGR(SLIDE)*dslidbinINV*2
      RVVISCOSITY(SLIDE)=-1.D0*RPZX(SLIDE)/RVGR(SLIDE)

!output the profiles    
  
      IF (tt.eq.0) then
       WRITE(193,*) 'TIMESTEP: ',IT
       write(139,*) 'TIMESTEP: ',IT
       write(183,*) 'TIMESTEP: ',IT
   DO J = 1,SLIDE
    VVISCOSITY(J) = VVISCOSITY(J)*MASSSCALE/RSCALE/TIMESCALE*1.d3 !in cP(mPa.s)
    RVVISCOSITY(J) = RVVISCOSITY(J)*MASSSCALE/RSCALE/TIMESCALE*1.d3
    VGR(J) = VGR(J)/TIMESCALE
    RVGR(J) = RVGR(J)*AVTINV/TIMESCALE
 write(193,110) J,(J-0.5d0)*slidbin,RVXPRO(J)*AVTINV*VSCALE,RDENSITY(J)*AVTINV*DSCALE,RTEEMP(J)*AVTINV,RPZX(J)*AVTINV*(-1.d0)*PSCALE
 write(139,110) J,(J-0.5d0)*slidbin,VXPRO(J)*VSCALE,DENSITY(J)*DSCALE,TEEMP(J),PZX(J)*PSCALE*(-1.D0)
 write(183,220) J,(J-0.5D0)*slidbin, VGR(J),RVGR(J),PZX(J)*(-1.d0)*PSCALE,RPZX(J)*(-1.d0)*PSCALE*AVTINV,VVISCOSITY(J),RVVISCOSITY(J)
   ENDDO
      endif 
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       CALCULATE THE VISCOSITY   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      V_AVE1 = 0.d0
      V_AVE2 = 0.d0
      DO J = 1,SLIDE2
         V_AVE1 = V_AVE1 + VXPRO(J)
      ENDDO    
      DO J = SLIDE2+1,SLIDE
         V_AVE2 = V_AVE2 + VXPRO(J)
      ENDDO
     
      DEN1 = DEN1/(BOXX*BOXY*BOXZ2)
      DEN2 = DEN2/(BOXX*BOXY*BOXZ2)
      RDEN1 = RDEN1 + DEN1
      RDEN2 = RDEN2 + DEN2
      V_AVE1 = V_AVE1/DBLE(SLIDE2)
      V_AVE2 = V_AVE2/DBLE(SLIDE2)
      V_AVE1 = abs(V_AVE1)
      V_AVE2 = abs(V_AVE2)
      RV_AVE1 = RV_AVE1 + V_AVE1
      RV_AVE2 = RV_AVE2 + V_AVE2
      
      VISCOSITY1 = DEN1*EXTER_FX*BOXZ2**2/(12.D0*V_AVE1)    
      VISCOSITY2 = DEN2*EXTER_FX*BOXZ2**2/(12.D0*V_AVE2)
      VISCOSITY = (VISCOSITY1+VISCOSITY2)*0.5D0

      RVISCOSITY1 = RDEN1*EXTER_FX*BOXZ2**2/(12.D0*RV_AVE1)
      RVISCOSITY2 = RDEN2*EXTER_FX*BOXZ2**2/(12.D0*RV_AVE2)
      RVISCOSITY = (RVISCOSITY1+RVISCOSITY2)*0.5D0
      
      VISCOSITY = VISCOSITY*MASSSCALE/RSCALE/TIMESCALE*1.d3         !in cP(mPa.s)
      RVISCOSITY = RVISCOSITY*MASSSCALE/RSCALE/TIMESCALE*1.d3       !in cP(mPa.s)
      
      IF (IT.EQ.TSET) THEN
         RVXPRO = 0.D0
         RTEEMP = 0.D0
         RPZX   = 0.D0
         RDENSITY =0.D0
         RDEN1 = 0.D0
         RDEN2 = 0.D0
         RV_AVE1 = 0.d0
         RV_AVE2 = 0.d0
      ENDIF

 110  FORMAT(I3,5(1X,1e12.5))
 220  FORMAT(I3,1x,1E12.5,6(1X,1e12.5))
      return
      end

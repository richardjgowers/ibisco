#include "ibi-preprocess.h"

      SUBROUTINE OUTPUT (I)
      USE VAR

      IMPLICIT NONE
      INTEGER :: I
      REAL*8 :: TREAL
      REAL*8 DENSP, RDENSP
      REAL*8 VOLP, RVOLP, PT23P, RPT23P, PT13P, RPT13P, PT12P, RPT12P
      REAL*8 PT33P, RPT33P, PT22P, RPT22P, PT11P, RPT11P, PRESSP, RPP
      REAL*8 TP, RTEMPP, VTORP, RVTORP, VANGLEP, RVANGLEP, VBONDP, RVBONDP
      REAL*8 VNBONDP, RVNBONDP, EKP, REKP, TOTPOTP, RVP, TOTEP, REP
      REAL*8 RVOOPP, VOOPP 
    real(kind=rkind) :: VBONDP_MIX, VANGLEP_MIX, VTORP_MIX, VOOPP_MIX, VNBONDP_MIX
    real(kind=rkind) :: RVBONDP_At, RVANGLEP_At, RVTORP_At, RVOOPP_At, RVNBONDP_At
    real(kind=rkind) :: VNBONDP_CG, VBONDP_CG, VANGLEP_CG, VTORP_CG, VOOPP_CG
    real(kind=rkind) :: RVBONDP_MIX, RVANGLEP_MIX, RVTORP_MIX, RVOOPP_MIX, RVNBONDP_MIX
    real(kind=rkind) :: RVNBONDP_CG, RVBONDP_CG, RVANGLEP_CG, RVTORP_CG, RVOOPP_CG
    real(kind=rkind) :: totVNBOND, totVBOND, totVANGLE, totVTORP, totVOOPP, totVTORS
#ifdef DEBUG_OOP
    real(kind=rkind) :: VOOPP_L,VOOPP_LD,VOOPP_D,VOOPP_ring
#endif
!      REAL*8 DENSP, RDENSP, BOXZP, RBOXZP,BOXYP, RBOXYP, BOXXP, RBOXXP
!      REAL*8 :: TREAL, TOTE
!       *******************************************************************

      OPEN (22, FILE = 'timestep')
      TREAL = I * DT * TIMESCALE * 1.0D+12 + INITIME !ps
      TP = EK * MKTEMP * TEMPSCALE
!      TOTPOT = 0.0
            
!      TOTPOT = VNBOND + VBOND + VANGLE + VTOR
!      WRITE (100,'(4e24.10)') TREAL, EK*CONV, TOTPOT*CONV, EK*CONV+TOTPOT*CONV
!      WRITE (101,'(5e24.10)') TREAL, VNBOND*CONV, VBOND*CONV, VANGLE*CONV, VTOR*CONV
!      WRITE (110,*) TREAL , TP
      WRITE(*,*)I, TP ,(PT11+PT22+PT33)*PSCALE/3.0d0
!      write(*,*)I, (PT11+PT22+PT33)*PSCALE/1000.0D0/3.0d0
!      write(111,*)TREAL ,(PT11+PT22+PT33)/3.0d0

      EKP = EK *CONV
      REKP = REK *CONV
!Total

!      VBONDP = VBOND *CONV
      RVBONDP = RVBOND *CONV
      
!      VANGLEP = VANGLE *CONV
      RVANGLEP = RVANGLE *CONV

!      VTORP = VTOR *CONV
      RVTORP = RVTOR *CONV

!      VOOPP = VOOP *CONV
      RVOOPP = RVOOP *CONV

!      VNBONDP = VNBOND *CONV
      RVNBONDP = RVNBOND *CONV

!  Mixed
      VBONDP_MIX = VBOND_MIX *CONV
      RVBONDP_MIX = RVBOND_MIX *CONV
      
      VANGLEP_MIX = VANGLE_MIX *CONV
      RVANGLEP_MIX = RVANGLE_MIX *CONV

      VTORP_MIX = VTOR_MIX *CONV
      RVTORP_MIX = RVTOR_MIX *CONV

      VOOPP_MIX = VOOP_MIX *CONV
      RVOOPP_MIX = RVOOP_MIX *CONV

      VNBONDP_MIX = VNBOND_MIX *CONV
      RVNBONDP_MIX = RVNBOND_MIX *CONV

! Beads

      VBONDP_CG = VBOND_CG *CONV
      RVBONDP_CG = RVBOND_CG *CONV
      
      VANGLEP_CG = VANGLE_CG *CONV
      RVANGLEP_CG = RVANGLE_CG *CONV

      VTORP_CG = VTOR_CG *CONV
      RVTORP_CG = RVTOR_CG *CONV

      VOOPP_CG = VOOP_CG *CONV
      RVOOPP_CG = RVOOP_CG *CONV

      VNBONDP_CG = VNBOND_CG *CONV
      RVNBONDP_CG = RVNBOND_CG *CONV

! Atoms

      VBONDP = VBOND *CONV
      RVBONDP_At = RVBOND_At *CONV
      
      VANGLEP = VANGLE *CONV
      RVANGLEP_At = RVANGLE_At *CONV

      VTORP = VTOR *CONV
      RVTORP_At = RVTOR_At *CONV

      VOOPP = VOOP *CONV
      RVOOPP_At = RVOOP_At *CONV

      VNBONDP = VNBOND *CONV
      RVNBONDP_At = RVNBOND_At *CONV

      TOTPOTP = VNBONDP + VBONDP + VANGLEP + VTORP + VOOP 
      TOTPOTP = TOTPOTP + VNBONDP_CG + VBONDP_CG + VANGLEP_CG + VTORP_CG + VOOP_CG
      TOTPOTP = TOTPOTP + VNBONDP_MIX + VBONDP_MIX + VANGLEP_MIX + VTORP_MIX + VOOP_MIX
      RVP = RV *CONV

#ifdef DEBUG_OOP
      VOOPP_ring = VOOP_ring *CONV
      VOOPP_LD = VOOP_LD *CONV
      VOOPP_L = VOOP_L *CONV
      VOOPP_D = VOOP_D *CONV
#endif

      TOTEP = TOTPOTP + EKP
      REP = RE *CONV

      RTEMPP = RTEMP * TEMPSCALE

      PT11P = PT11 *PSCALE
      RPT11P = RPT11 *PSCALE

      PT22P = PT22 *PSCALE
      RPT22P = RPT22 *PSCALE

      PT33P = PT33 *PSCALE
      RPT33P = RPT33 *PSCALE

      PT12P = PT12 *PSCALE
      RPT12P = RPT12 *PSCALE

      PT13P = PT13 *PSCALE
      RPT13P = RPT13 *PSCALE

      PT23P = PT23 *PSCALE
      RPT23P = RPT23 *PSCALE

      PRESSP = (PT11P+PT22P+PT33P)/3.0d0
      RPP = RP *PSCALE

      VOLP = BOXX * BOXY * BOXZ 
      RVOLP = RVOL 

      DENSP = TOTMASS  * DSCALE / VOLP
      RDENSP = RDENS * DSCALE

    totVBOND  = VBONDP+VBONDP_MIX+VBONDP_CG
    totVANGLE = VANGLEP+VANGLEP_MIX+VANGLEP_CG
    TotVNBOND = VNBONDP+VNBONDP_MIX+VNBONDP_CG
    TotVTORS  = VTORP+VTORP_MIX+VTORP_CG


      WRITE (115, *)'Step:                     ', I
      WRITE (115, 100)'Simulated_time:           ', TREAL
      WRITE (115, 100)'Total_energy:             ', TOTEP, REP
      WRITE (115, 100)'Potential_energy:         ', TOTPOTP, RVP
      WRITE (115, 100)'Kinetic_energy:           ', EKP, REKP
      WRITE (115, 100)'Tot._Nonbonded_energy:         ', totVNBOND, RVNBONDP
      WRITE (115, 100)'      Nonbonded_Atom_energy:  ', VNBONDP, RVNBONDP_At
      WRITE (115, 100)'      Nonbonded_Beads_energy: ', VNBONDP_CG, RVNBONDP_CG
      WRITE (115, 100)'      Nonbonded_mix_energy:   ', VNBONDP_MIX, RVNBONDP_MIX
      WRITE (115, 100)'Tot._Bond_energy:        ', totVBOND , RVBONDP
      WRITE (115, 100)'      Bond_Atom_energy:  ', VBONDP, RVBONDP_At
      WRITE (115, 100)'      Bond_Beads_energy: ', VBONDP_CG, RVBONDP_CG
      WRITE (115, 100)'      Bond_mix_energy:   ', VBONDP_MIX, RVBONDP_MIX
      WRITE (115, 100)'Tot._Angle_energy:         ', totVANGLE, RVANGLEP
      WRITE (115, 100)'      Angle_Atom_energy:  ', VAngleP, RVAngleP_At
      WRITE (115, 100)'      Angle_Beads_energy: ', VAngleP_CG, RVAngleP_CG
      WRITE (115, 100)'      Angle_mix_energy:   ', VAngleP_MIX, RVAngleP_MIX
      WRITE (115, 100)'Tot._Torsion energy:       ', totVTORS, RVTORP
      WRITE (115, 100)'     Torsion_Atom_energy: ', VTORP, RVTORP_At
      WRITE (115, 100)'     Torsion_Beads_energy:', VTORP_CG, RVTORP_CG
      WRITE (115, 100)'     Torsion_mix_energy:  ', VTORP_MIX, RVTORP_MIX

      WRITE (115, 100)'Improper_torsion_energy:  ', VOOPP, RVOOPP
#ifdef DEBUG_OOP
      WRITE (115, *)'           Ring:          ', VOOPP_ring
      WRITE (115, *)'             LD:          ', VOOPP_LD
      WRITE (115, *)'             L:          ', VOOPP_L
      WRITE (115, *)'             D:          ', VOOPP_D
#endif
      WRITE (115, 100)'Temperature:              ', TP, RTEMPP
      WRITE (115, 100)'Pressure:                 ', PRESSP, RPP
      WRITE (115, 100)'Pressure(x):              ', PT11P, RPT11P
      WRITE (115, 100)'Pressure(y):              ', PT22P, RPT22P
      WRITE (115, 100)'Pressure(z):              ', PT33P, RPT33P
      WRITE (115, 100)'Pressure(xy):             ', PT12P, RPT12P
      WRITE (115, 100)'Pressure(xz):             ', PT13P, RPT13P
      WRITE (115, 100)'Pressure(yz):             ', PT23P, RPT23P
      WRITE (115, 100)'Box_volume:               ', VOLP, RVOLP
      WRITE (115, 100)'Box_length(x):            ', BOXX, RBOXX
      WRITE (115, 100)'Box_length(y):            ', BOXY, RBOXY
      WRITE (115, 100)'Box_length(z):            ', BOXZ, RBOXZ
      WRITE (115, 100)'Mass_density:             ', DENSP, RDENSP
      WRITE (115, 100)
      WRITE (115, 100)

      WRITE (22, *)'Step:                      ', I
      WRITE (22, *)'Simulated_time:            ', TREAL
      WRITE (22, *)'Total_energy:              ', TOTEP, REP
      WRITE (22, *)'Potential_energy:          ', TOTPOTP, RVP
      WRITE (22, *)'Kinetic_energy:            ', EKP, REKP
      WRITE (22, *)'Nonbonded_energy:          ', VNBONDP, RVNBONDP
      WRITE (22, *)'Bond_energy:              ', VBONDP, RVBONDP
      WRITE (22, *)'Angle_energy:              ', VANGLEP, RVANGLEP
      WRITE (22, *)'Torsion_energy:            ', VTORP, RVTORP
      WRITE (22, *)'Improper_torsion_energy:   ', VOOPP, RVOOPP
      WRITE (22, *)'Temperature:               ', TP, RTEMPP
      WRITE (22, *)'Pressure:                  ', PRESSP, RPP
      WRITE (22, *)'Pressure(x):               ', PT11P, RPT11P
      WRITE (22, *)'Pressure(y):               ', PT22P, RPT22P
      WRITE (22, *)'Pressure(z):               ', PT33P, RPT33P
      WRITE (22, *)'Pressure(xy):              ', PT12P, RPT12P
      WRITE (22, *)'Pressure(xz):              ', PT13P, RPT13P
      WRITE (22, *)'Pressure(yz):              ', PT23P, RPT23P
      WRITE (22, *)'Box_volume:                ', VOLP, RVOLP
      WRITE (22, *)'Box_length(x):             ', BOXX, RBOXX
      WRITE (22, *)'Box_length(y):             ', BOXY, RBOXY
      WRITE (22, *)'Box_length(z):             ', BOXZ, RBOXZ
      WRITE (22, *)'Mass_density:              ', DENSP, RDENSP

        CLOSE (22)      

100 format (A,2(F16.5,2X))

      RETURN
      END

!      *********************************************************************************************

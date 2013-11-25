#include "ibi-preprocess.h"

SUBROUTINE OUTPUT (I)
  USE VAR

  IMPLICIT NONE
  INTEGER :: I
  REAL*8 :: TREAL
  REAL*8 :: AV_CONV, AV_PSCALE, INV_AV
  REAL*8 DENSP, RDENSP
  REAL*8 VOLP, RVOLP, PT23P, RPT23P, PT13P, RPT13P, PT12P, RPT12P
  REAL*8 PT33P, RPT33P, PT22P, RPT22P, PT11P, RPT11P, PRESSP, RPP
  REAL*8 TP, RTEMPP, VTORP, RVTORP, VANGLEP, RVANGLEP, VBONDP, RVBONDP
  REAL*8 VNBONDP, RVNBONDP, TOTPOTP, RVP, TOTEP, REP
  REAL*8 RVOOPP, VOOPP 
#ifdef DEBUG_OOP
  real(kind=rkind) :: VOOPP_L,VOOPP_LD,VOOPP_D,VOOPP_ring
#endif

  INV_AV = 1.0 / NSAMPLING
  AV_CONV = INV_AV * CONV
  AV_PSCALE = INV_AV * PSCALE

  TREAL = I * DT * TIMESCALE * 1.0D+12 + INITIME !ps
  TP = EK * MKTEMP * TEMPSCALE

  WRITE(*,*) I, TEMP*TEMPSCALE ,PRES(1)*PSCALE

  !Total

  !      VBONDP = VBOND *CONV
  RVBONDP = RVBOND *CONV

  !      VANGLEP = VANGLE *CONV
  RVANGLEP = RVANGLE *CONV

  !      VTORP = VTOR *CONV
  RVTORP = RVTOR *CONV

  !      VOOPP = VOOP *CONV
  RVOOPP = RVOOP *CONV

  !      VNBONDP = (VNBOND_BEAD + VNBOND_ATOM) *CONV
  RVNBONDP = RVNBOND *CONV

  ! Beads
!  VBONDP_CG = V_BOND(2)*CONV
!  RVBONDP_CG = RVBOND_CG *CONV

 ! VANGLEP_CG = V_ANGLE(2)*CONV
 ! RVANGLEP_CG = RVANGLE_CG *CONV

  !VTORP_CG = V_TORSION(2) *CONV
  !RVTORP_CG = RVTOR_CG *CONV

!  VOOPP_CG = V_OOP(2)*CONV
!  RVOOPP_CG = RVOOP_CG *CONV

!  VNBONDP_BEAD = V_NB(2) *CONV
!  RVNBONDP_CG = RVNBOND_CG *CONV

  ! Atoms
  VBONDP = V_BOND(1)*CONV

  VANGLEP = V_ANGLE(1)*CONV

  VTORP = V_TORSION(1)*CONV

  VOOPP = V_OOP(1)*CONV

  VNBONDP = V_NB(1)*CONV

  TOTPOTP = VNBONDP + VBONDP + VANGLEP + VTORP + VOOP 

  RVP = RV *CONV


!  VANGLEP_MIX = 0.0D0

#ifdef DEBUG_OOP
  VOOPP_ring = VOOP_ring *CONV
  VOOPP_LD = VOOP_LD *CONV
  VOOPP_L = VOOP_L *CONV
  VOOPP_D = VOOP_D *CONV
#endif

!  TOTEP = TOTPOTP + EKP
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

  WRITE (115, *)'Step:                     ', I
  WRITE (115, 100)'Simulated_time:           ', TREAL
  WRITE (115, 100)'Total_energy:             ', TOT_E*CONV, AV_TOT_E*AV_CONV
  WRITE (115, 100)'Potential_energy:         ', POT_E*CONV, AV_POT_E*AV_CONV
  WRITE (115, 100)'Kinetic_energy:           ', KIN_E*CONV, AV_KIN_E*AV_CONV
  WRITE (115, 100)'Tot._Nonbonded_energy:         ', SUM(V_NB)*CONV, SUM(AV_V_NB)*AV_CONV 
  WRITE (115, 100)'      Nonbonded_Atom_energy:  ', V_NB(1)*CONV, AV_V_NB(1)*AV_CONV
  WRITE (115, 100)'      Nonbonded_Beads_energy: ', V_NB(2)*CONV, AV_V_NB(2)*AV_CONV
  WRITE (115, 100)'      Nonbonded_mix_energy:   ', V_NB(3)*CONV, AV_V_NB(3)*AV_CONV
  WRITE (115, 100)'Tot._Bond_energy:         ', SUM(V_BOND)*CONV, SUM(AV_V_BOND)*AV_CONV
  WRITE (115, 100)'      Bond_Atom_energy:  ', V_BOND(1)*CONV, AV_V_BOND(1)*AV_CONV
  WRITE (115, 100)'      Bond_Beads_energy: ', V_BOND(2)*CONV, AV_V_BOND(2)*AV_CONV
  WRITE (115, 100)'      Bond_mix_energy:   ', V_BOND(3)*CONV, AV_V_BOND(3)*AV_CONV
  WRITE (115, 100)'Tot._Angle_energy:         ', SUM(V_ANGLE)*CONV, SUM(AV_V_ANGLE)*AV_CONV
  WRITE (115, 100)'      Angle_Atom_energy:  ', V_ANGLE(1)*CONV, AV_V_ANGLE(1)*AV_CONV
  WRITE (115, 100)'      Angle_Beads_energy: ', V_ANGLE(2)*CONV, AV_V_ANGLE(2)*AV_CONV
  WRITE (115, 100)'      Angle_mix_energy:   ', V_ANGLE(3)*CONV, AV_V_ANGLE(3)*AV_CONV
  WRITE (115, 100)'Tot._Torsion energy:       ', SUM(V_TORSION)*CONV, SUM(AV_V_TORSION)*AV_CONV
  WRITE (115, 100)'     Torsion_Atom_energy: ', V_TORSION(1)*CONV, AV_V_TORSION(1)*AV_CONV
  WRITE (115, 100)'     Torsion_Beads_energy:', V_TORSION(2)*CONV, AV_V_TORSION(2)*AV_CONV
  WRITE (115, 100)'     Torsion_mix_energy:  ', V_TORSION(3)*CONV, AV_V_TORSION(3)*AV_CONV
  WRITE (115, 100)'Improper_torsion_energy:   ', SUM(V_OOP)*CONV, SUM(AV_V_OOP)*AV_CONV
#ifdef DEBUG_OOP
  WRITE (115, *)'           Ring:          ', VOOPP_ring
  WRITE (115, *)'             LD:          ', VOOPP_LD
  WRITE (115, *)'             L:          ', VOOPP_L
  WRITE (115, *)'             D:          ', VOOPP_D
#endif
  WRITE (115, 100)'Temperature:              ', TEMP*TEMPSCALE, AV_TEMP*INV_AV*TEMPSCALE
  WRITE (115, 100)'Pressure:                 ', PRES(1)*PSCALE, AV_PRES(1)*AV_PSCALE
  WRITE (115, 100)'Pressure(x):              ', PRES(2)*PSCALE, AV_PRES(2)*AV_PSCALE
  WRITE (115, 100)'Pressure(y):              ', PRES(3)*PSCALE, AV_PRES(3)*AV_PSCALE
  WRITE (115, 100)'Pressure(z):              ', PRES(4)*PSCALE, AV_PRES(4)*AV_PSCALE
  WRITE (115, 100)'Pressure(xy):             ', PRES(5)*PSCALE, AV_PRES(5)*AV_PSCALE
  WRITE (115, 100)'Pressure(xz):             ', PRES(6)*PSCALE, AV_PRES(6)*AV_PSCALE
  WRITE (115, 100)'Pressure(yz):             ', PRES(7)*PSCALE, AV_PRES(7)*AV_PSCALE
  WRITE (115, 100)'Box_volume:               ', BOX(1), AV_BOX(1)*INV_AV
  WRITE (115, 100)'Box_length(x):            ', BOX(2), AV_BOX(2)*INV_AV
  WRITE (115, 100)'Box_length(y):            ', BOX(3), AV_BOX(3)*INV_AV
  WRITE (115, 100)'Box_length(z):            ', BOX(4), AV_BOX(4)*INV_AV
  WRITE (115, 100)'Mass_density:             ', DENS * DSCALE, AV_DENS*INV_AV*DSCALE
  WRITE (115, 100)
  WRITE (115, 100)

100 format (A,2(F16.5,2X))

  OPEN (22, FILE = 'timestep')

  WRITE (22, *)'Step:                      ', I
  WRITE (22, *)'Simulated_time:            ', TREAL
  WRITE (22, *)'Total_energy:              ', TOT_E*CONV, AV_TOT_E*AV_CONV
  WRITE (22, *)'Potential_energy:          ', POT_E*CONV, AV_POT_E*AV_CONV
  WRITE (22, *)'Kinetic_energy:            ', KIN_E*CONV, AV_KIN_E*AV_CONV
  WRITE (22, *)'Nonbonded_energy:          ', SUM(V_NB)*CONV, SUM(AV_V_NB)*AV_CONV
  WRITE (22, *)'Bond_energy:               ', SUM(V_BOND)*CONV, SUM(AV_V_BOND)*AV_CONV
  WRITE (22, *)'Angle_energy:              ', SUM(V_ANGLE)*CONV, SUM(AV_V_ANGLE)*AV_CONV
  WRITE (22, *)'Torsion_energy:            ', SUM(V_TORSION)*CONV, SUM(AV_V_TORSION)*AV_CONV
  WRITE (22, *)'Improper_torsion_energy:   ', SUM(V_OOP)*CONV, SUM(AV_V_OOP)*AV_CONV
  WRITE (22, *)'Temperature:               ', TEMP*TEMPSCALE, AV_TEMP*INV_AV*TEMPSCALE
  WRITE (22, *)'Pressure:                  ', PRES(1)*PSCALE, AV_PRES(1)*AV_PSCALE
  WRITE (22, *)'Pressure(x):               ', PRES(2)*PSCALE, AV_PRES(2)*AV_PSCALE
  WRITE (22, *)'Pressure(y):               ', PRES(3)*PSCALE, AV_PRES(3)*AV_PSCALE
  WRITE (22, *)'Pressure(z):               ', PRES(4)*PSCALE, AV_PRES(4)*AV_PSCALE
  WRITE (22, *)'Pressure(xy):              ', PRES(5)*PSCALE, AV_PRES(5)*AV_PSCALE
  WRITE (22, *)'Pressure(xz):              ', PRES(6)*PSCALE, AV_PRES(6)*AV_PSCALE
  WRITE (22, *)'Pressure(yz):              ', PRES(7)*PSCALE, AV_PRES(7)*AV_PSCALE
  WRITE (22, *)'Box_volume:                ', BOX(1), AV_BOX(1)*INV_AV
  WRITE (22, *)'Box_length(x):             ', BOX(2), AV_BOX(2)*INV_AV
  WRITE (22, *)'Box_length(y):             ', BOX(3), AV_BOX(3)*INV_AV
  WRITE (22, *)'Box_length(z):             ', BOX(4), AV_BOX(4)*INV_AV
  WRITE (22, *)'Mass_density:              ', DENS * DSCALE, AV_DENS*INV_AV*DSCALE

  CLOSE (22)      

  !Reset average totals to 0 to begin new averaging window
  AV_TOT_E = 0.0
  AV_POT_E = 0.0
  AV_KIN_E = 0.0
  AV_V_NB = 0.0
  AV_V_BOND = 0.0
  AV_V_ANGLE = 0.0
  AV_V_TORSION = 0.0
  AV_V_OOP = 0.0
  AV_TEMP = 0.0
  AV_PRES = 0.0
  AV_BOX = 0.0
  AV_DENS = 0.0

  RETURN
END SUBROUTINE OUTPUT

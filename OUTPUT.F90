#ifndef DOXY_SKIP
#include "ibi-preprocess.h"
#endif
!> @file
!> @brief Writes the s-md.out file

SUBROUTINE OUTPUT (I)

  USE VAR

  IMPLICIT NONE

  INTEGER :: I
  REAL*8 :: TREAL
  REAL*8 :: AV_CONV, AV_PSCALE, INV_AV
#ifdef DEBUG_OOP
  real*4 :: VOOPP_L,VOOPP_LD,VOOPP_D,VOOPP_ring
#endif

  INV_AV = 1.0 / NTRJ
  AV_CONV = INV_AV * CONV
  AV_PSCALE = INV_AV * PSCALE

  TREAL = I * DT * TIMESCALE * 1.0D+12 + INITIME !ps

  WRITE(*,*) I, TEMP*TEMPSCALE, PRES(1)*PSCALE

#ifdef DEBUG_OOP
  VOOPP_ring = VOOP_ring *CONV
  VOOPP_LD = VOOP_LD *CONV
  VOOPP_L = VOOP_L *CONV
  VOOPP_D = VOOP_D *CONV
#endif

  WRITE (115, *)'Step:                       ', I
  WRITE (115, 100)'Simulated_time:           ', TREAL
  WRITE (115, 100)'Total_energy:             ', TOT_E*CONV, AV_TOT_E*AV_CONV
  WRITE (115, 100)'Potential_energy:         ', POT_E*CONV, AV_POT_E*AV_CONV
  WRITE (115, 100)'Kinetic_energy:           ', SUM(KIN_E)*CONV, SUM(AV_KIN_E)*AV_CONV
  WRITE (115, 100)'      Kinetic_energy_atom:', KIN_E(1)*CONV, AV_KIN_E(1)*CONV
  WRITE (115, 100)'      Kinetic_energy_bead:', KIN_E(2)*CONV, AV_KIN_E(2)*CONV
  WRITE (115, 100)'Tot._Nonbonded_energy:         ', SUM(V_NB)*CONV, SUM(AV_V_NB)*AV_CONV 
  WRITE (115, 100)'      Nonbonded_Atom_energy:   ', V_NB(1)*CONV, AV_V_NB(1)*AV_CONV
  WRITE (115, 100)'      Nonbonded_Beads_energy:  ', V_NB(2)*CONV, AV_V_NB(2)*AV_CONV
  WRITE (115, 100)'      Nonbonded_mix_energy:    ', V_NB(3)*CONV, AV_V_NB(3)*AV_CONV
  WRITE (115, 100)'Tot._Bond_energy:         ', SUM(V_BOND)*CONV, SUM(AV_V_BOND)*AV_CONV
  WRITE (115, 100)'      Bond_Atom_energy:   ', V_BOND(1)*CONV, AV_V_BOND(1)*AV_CONV
  WRITE (115, 100)'      Bond_Beads_energy:  ', V_BOND(2)*CONV, AV_V_BOND(2)*AV_CONV
  WRITE (115, 100)'      Bond_mix_energy:    ', V_BOND(3)*CONV, AV_V_BOND(3)*AV_CONV
  WRITE (115, 100)'Tot._Angle_energy:         ', SUM(V_ANGLE)*CONV, SUM(AV_V_ANGLE)*AV_CONV
  WRITE (115, 100)'      Angle_Atom_energy:   ', V_ANGLE(1)*CONV, AV_V_ANGLE(1)*AV_CONV
  WRITE (115, 100)'      Angle_Beads_energy:  ', V_ANGLE(2)*CONV, AV_V_ANGLE(2)*AV_CONV
  WRITE (115, 100)'      Angle_mix_energy:    ', V_ANGLE(3)*CONV, AV_V_ANGLE(3)*AV_CONV
  WRITE (115, 100)'Tot._Torsion energy:       ', SUM(V_TORSION)*CONV, SUM(AV_V_TORSION)*AV_CONV
  WRITE (115, 100)'     Torsion_Atom_energy:  ', V_TORSION(1)*CONV, AV_V_TORSION(1)*AV_CONV
  WRITE (115, 100)'     Torsion_Beads_energy: ', V_TORSION(2)*CONV, AV_V_TORSION(2)*AV_CONV
  WRITE (115, 100)'     Torsion_mix_energy:   ', V_TORSION(3)*CONV, AV_V_TORSION(3)*AV_CONV
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

  FLUSH(115)

100 format (A,2(F16.5,2X))

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

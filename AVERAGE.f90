!> @file 
!> @brief Keeps track of the average values for various things
SUBROUTINE AVERAGE()

  USE VAR

  IMPLICIT NONE

  ! Calculate a few totals, why not
  POT_E = SUM(V_NB) + SUM(V_BOND) + SUM(V_ANGLE) + SUM(V_TORSION) + SUM(V_OOP)
  KIN_E(1) = EK(1)
  KIN_E(2) = EK(2)
  TOT_E = SUM(KIN_E) + POT_E
  
  PRES(1) = (PT11 + PT22 + PT33)/3.0
  PRES(2) = PT11
  PRES(3) = PT22
  PRES(4) = PT33
  PRES(5) = PT12
  PRES(6) = PT13
  PRES(7) = PT23

  BOXSIZE(1) = BOX(1) * BOX(2) * BOX(3)
  BOXSIZE(2) = BOX(1)
  BOXSIZE(3) = BOX(2)
  BOXSIZE(4) = BOX(3)

  DENS = TOTMASS / BOX(1)

! Add to totals which later form averages
! These totals are reset to 0 in OUTPUT when they are used
! Energies
  AV_TOT_E = AV_TOT_E + TOT_E
  AV_POT_E = AV_POT_E + POT_E
  AV_KIN_E = AV_KIN_E + KIN_E
! Components, this is array addition
  AV_V_NB = AV_V_NB + V_NB
  AV_V_BOND = AV_V_BOND + V_BOND
  AV_V_ANGLE = AV_V_ANGLE + V_ANGLE
  AV_V_TORSION = AV_V_TORSION + V_TORSION
  AV_V_OOP = AV_V_OOP + V_OOP
! Temperature and pressures
  AV_TEMP = AV_TEMP + TEMP
  AV_PRES = AV_PRES + PRES
! Box volume & density
  AV_BOXSIZE = AV_BOXSIZE + BOXSIZE
  AV_DENS = AV_DENS + DENS

 RETURN
END SUBROUTINE AVERAGE

!> @file
!> @brief Defines some constants

SUBROUTINE UNIT ()
  USE VAR

  IMPLICIT NONE

  MMOLY = 0.012D0         ! KG / MOLE   C

  MASSSCALE = MMOLY / NA	!KG
  RSCALE = 1.0D-9		!METER
  ESCALE = KB * 300.0	!JOUL

  TEMPSCALE = ESCALE / KB 	          !TEMP*TEMPSCALE   =>  convert reduced unit to Kelvin
  !	PSCALE = ESCALE / (RSCALE)**3.0
  PSCALE = ESCALE / (RSCALE)**3.0/1000.0D0  !PRESSURE*PSCALE  =>  convert reduced unit to KPascal
  TIMESCALE = RSCALE*SQRT(MASSSCALE/ESCALE) !TIME*TIMESCALE   =>  convert reduced unit to Second
  FSCALE = ESCALE / RSCALE		  !FORCE*FSCALE     =>  convert reduced unit to Newton
  DSCALE = MASSSCALE / (RSCALE)**3.0        !RHO*DSCALE       =>  convert reduced unit to KG/M**3
  VSCALE = RSCALE / TIMESCALE		  !V * VSCALE	    =>  convert reduced unit to M/S
  CONV = NA*ESCALE/1000.0D0		  !E*CONV	    =>	convert reduced unit to kJoul/Mole

  R2D = 180.0 / PI
  D2R = PI / 180.0

  RETURN
END SUBROUTINE UNIT

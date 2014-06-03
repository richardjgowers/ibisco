!> @file
!> @brief MTS Method module
!> @details All variables and subroutines unique to MTS methodology.

MODULE MTS

  USE VAR, ONLY : RKIND

  INTEGER :: NMTS ! Number of steps to estimate the CG force
  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE :: MTS_FX, MTS_FY, MTS_FZ ! Holding arrays for force
  REAL(KIND=RKIND), DIMENSION(:), ALLOCATABLE :: MTS_PT11, MTS_PT22, MTS_PT33, MTS_PT12, MTS_PT13, MTS_PT23
  REAL(KIND=RKIND), DIMENSION(:,:), ALLOCATABLE :: MTS_V_NB

END MODULE MTS

!> @file
!> @brief MTS Method module
!> @details All variables and subroutines unique to MTS methodology.

MODULE MTS

  INTEGER :: NMTS ! Number of steps to estimate the CG force
  REAL*4, DIMENSION(:,:,:), ALLOCATABLE :: MTS_FXYZ
  REAL*4, DIMENSION(:), ALLOCATABLE :: MTS_PT11, MTS_PT22, MTS_PT33, MTS_PT12, MTS_PT13, MTS_PT23
  REAL*4, DIMENSION(:,:), ALLOCATABLE :: MTS_V_NB

END MODULE MTS

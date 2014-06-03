!> @file
!> @brief Write trajectory frame to TRZ format
!!
!> @details TRZ format can be read either by the YASP analysis tools, or the python package \n
!! MDAnalysis.

SUBROUTINE WRITETRJ (TM)

  USE VAR

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: TM !< Current time step
  INTEGER :: I, J
  REAL*8 :: PRESS, PT_11, PT_12, PT_13, PT_22, PT_23, PT_33
  REAL*8 :: PTOT, ETOT, T, TREAL
  REAL*8 :: BOX_X, BOX_Y, BOX_Z
  REAL*8 :: smallnumber
  REAL*4, POINTER :: NSX(:), NSY(:), NSZ(:), NVX(:), NVY(:), NVZ(:)
  INTEGER*4 :: nrec
  CHARACTER(80) :: YASPTITLE 

  data nrec / 10 /

  NFRAME = NFRAME + 1

  ALLOCATE(NSX(NATOMS))
  ALLOCATE(NSY(NATOMS))
  ALLOCATE(NSZ(NATOMS))
  ALLOCATE(NVX(NATOMS))
  ALLOCATE(NVY(NATOMS))
  ALLOCATE(NVZ(NATOMS))

  !	CONVERT TO 4-BYTE PRECISION (COORDINATES AND VELOCITES ONLY)
  smallnumber = 1.0D0 / 1.0D+12 / TIMESCALE

  DO J = 1, NATOMS
     NSX(J) = SX(J)
     NSY(J) = SY(J)
     NSZ(J) = SZ(J)

     NVX(J) = VX(J) *smallnumber !/ TIMESCALE / 1.0D+12      !nanometer/ps
     NVY(J) = VY(J) *smallnumber !/ TIMESCALE / 1.0D+12      !nanometer/ps
     NVZ(J) = VZ(J) *smallnumber !/ TIMESCALE / 1.0D+12      !nanometer/ps
  END DO

  !	write header record (together with first frame)
  IF (NFRAME == 1) THEN
     YASPTITLE = TITLE
     WRITE(113)YASPTITLE
     WRITE(113)nrec
  END IF
  TREAL = TM * DT * TIMESCALE * 1.0D+12 + INITIME
  !	... frame description (number of frame, the index number of the frame, number of atoms and time (ps))
  WRITE(113) NFRAME, NTRJ*NFRAME, NATOMS, TREAL

  !	... simulation cell unit vectors
  BOX_X = BOXX
  BOX_Y = BOXY
  BOX_Z = BOXZ
  WRITE(113) BOX_X, 0.0D00, 0.0D00, &
       0.0D00, BOX_Y, 0.0D00, &
       0.0D00, 0.0D00, BOX_Z

  !     ... isotropic pressure and pressure tensor
  PRESS = (PT11+PT22+PT33)*PSCALE/3.0d0
  PT_11 = PT11*PSCALE
  PT_12 = PT12*PSCALE
  PT_13 = PT13*PSCALE
  PT_22 = PT22*PSCALE
  PT_23 = PT23*PSCALE
  PT_33 = PT33*PSCALE
  WRITE(113) PRESS,      &
       PT_11,            &
       PT_12, PT_22,     &
       PT_13, PT_23, PT_33

  PTOT = SUM(V_BOND) + SUM(V_ANGLE) + SUM(V_TORSION) + SUM(V_OOP) + SUM(V_NB)
  ETOT = PTOT + SUM(EK)
  T = SUM(EK) * MKTEMP
  !	... energy and friends
  WRITE(113) 6, ETOT*CONV, PTOT*CONV, SUM(EK)*CONV, T, 0.0D00, 0.0D00

  !	... coordinates (atom positions)
  WRITE(113)NSX
  WRITE(113)NSY
  WRITE(113)NSZ

  !	... atom velocities
  WRITE(113)NVX
  WRITE(113)NVY
  WRITE(113)NVZ


  OPEN(UNIT=114, FILE = 'config.xyz', STATUS='replace')
  WRITE(114,*)NATOMS
  WRITE(114,*) 't'
  DO I = 1, NATOMS
     WRITE(114,9050) LABEL(ITYPE(I)), NSX(I)*10, NSY(I)*10, NSZ(I)*10
  END DO
  !9040 FORMAT (1X,A8,1X,3 (G21.14,1X))
9050 FORMAT (1X,A8,1X,3 (G21.14,1X))
  CLOSE (114)

  DEALLOCATE(NSX)
  DEALLOCATE(NSY)
  DEALLOCATE(NSZ)
  DEALLOCATE(NVX)
  DEALLOCATE(NVY)
  DEALLOCATE(NVZ)

  RETURN

END SUBROUTINE WRITETRJ

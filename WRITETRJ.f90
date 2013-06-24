SUBROUTINE WRITETRJ (TM)

  USE VAR

  IMPLICIT NONE

  !      INTEGER :: I,J,K,L,M,KK,TM
  INTEGER :: I,J,TM
  REAL*8 ::PRESS, PTOT, ETOT, T, TREAL
  REAL(KIND=RKIND) :: smallnumber
  REAL*4, POINTER :: NSX(:), NSY(:), NSZ(:), NVX(:), NVY(:), NVZ(:)
  INTEGER*4 :: nrec
  CHARACTER(len=3):: WRITETYPE
  CHARACTER(80) :: YASPTITLE

  data nrec / 10 /
  !       *******************************************************************

  NFRAME = NFRAME + 1

  ALLOCATE(NSX(NATOMS))
  ALLOCATE(NSY(NATOMS))
  ALLOCATE(NSZ(NATOMS))
  ALLOCATE(NVX(NATOMS))
  ALLOCATE(NVY(NATOMS))
  ALLOCATE(NVZ(NATOMS))



  !	CONVERT TO 4-BYTE PRECISION (COORDINATES AND VELOCITES ONLY)
  smallnumber = 1.0 / 1.0D+12 / TIMESCALE

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
  WRITE(113) BOXX, 0.0D00, 0.0D00, &
       0.0D00, BOXY, 0.0D00, &
       0.0D00, 0.0D00, BOXZ

  !     ... isotropic pressure and pressure tensor
  PRESS = (PT11+PT22+PT33)*PSCALE/3.0d0
  WRITE(113) PRESS,                   &
       PT11*PSCALE,                  &
       PT12*PSCALE, PT22*PSCALE,     &
       PT13*PSCALE, PT23*PSCALE, PT33*PSCALE

  PTOT = VNBOND_ATOM + VNBOND_BEAD + VBOND + VANGLE + VTOR
  ETOT = PTOT + EK
  T = EK * MKTEMP
  !	... energy and friends
  WRITE(113) 6, ETOT*CONV, PTOT*CONV, EK*CONV, T, 0.0D00, 0.0D00

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
     IF(IBRDESCR .EQ. 0) THEN
        IF(TYPE_LABEL(I) .EQ. 1) THEN
           WRITETYPE = 'C'
           WRITE(114,9040)WRITETYPE, NSX(I)*10, NSY(I)*10, NSZ(I)*10
        ELSE
           WRITETYPE = 'O'
           WRITE(114,9040)WRITETYPE, NSX(I)*10, NSY(I)*10, NSZ(I)*10
        END IF
     ELSE
        WRITE(114,9050)LABEL(ITYPE(I)), NSX(I), NSY(I), NSZ(I)
     END IF
  END DO
9040 FORMAT (1X,A8,1X,3 (G21.14,1X))
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

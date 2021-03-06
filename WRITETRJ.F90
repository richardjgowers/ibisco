#ifndef DOXY_SKIP
#include "ibi-preprocess.h"
#endif
!> @file
!> @brief Write trajectory frame to TRZ format.  Saves the restart file.
!!
!> @details TRZ format can be read either by the YASP analysis tools, or the python package \n
!! MDAnalysis.

SUBROUTINE WRITETRJ(TM)

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: TM !< Current time step
  INTEGER :: I, J
  REAL*8 :: PRESS, PT_11, PT_12, PT_13, PT_22, PT_23, PT_33
  REAL*8 :: PTOT, ETOT, T, TREAL
  REAL*8 :: BOX_X, BOX_Y, BOX_Z
  REAL*8 :: smallnumber
  INTEGER*4 :: nrec
  CHARACTER(80) :: YASPTITLE 

  data nrec / 10 /

  NFRAME = NFRAME + 1

  !	CONVERT TO 4-BYTE PRECISION (COORDINATES AND VELOCITES ONLY)
  smallnumber = 1.0D0 / 1.0D+12 / TIMESCALE

  ! Convert all values to be saved into correct precision
  !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC, 4)&
  !$OMP& SHARED(NSX, NSY, NSZ, SXYZ)&
  !$OMP& SHARED(NVX, NVY, NVZ, VXYZ)&
  !$OMP& SHARED(smallnumber, NATOMS)&
  !$OMP& PRIVATE(J)
  DO J = 1, NATOMS
     NSX(J) = SXYZ(1,J)
     NSY(J) = SXYZ(2,J)
     NSZ(J) = SXYZ(3,J)

     NVX(J) = VXYZ(1,J) *smallnumber !/ TIMESCALE / 1.0D+12      !nanometer/ps
     NVY(J) = VXYZ(2,J) *smallnumber !/ TIMESCALE / 1.0D+12      !nanometer/ps
     NVZ(J) = VXYZ(3,J) *smallnumber !/ TIMESCALE / 1.0D+12      !nanometer/ps
  END DO
  !$OMP END PARALLEL DO

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
  BOX_X = BOX(1)
  BOX_Y = BOX(2)
  BOX_Z = BOX(3)
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

  FLUSH(113)

  ! Write xyz file of last configuration
  OPEN(UNIT=114, FILE = 'config.xyz', STATUS='replace')
  WRITE(114,*)NATOMS
  WRITE(114,*) 't'
  DO I = 1, NATOMS
     WRITE(114,9050) LABEL(ITYPE(I)), NSX(I)*10, NSY(I)*10, NSZ(I)*10
  END DO
  !9040 FORMAT (1X,A8,1X,3 (G21.14,1X))
9050 FORMAT (1X,A8,1X,3 (G21.14,1X))
  CLOSE (114)


#ifndef SKIP_WRITE_RESTART
  CALL WRITE_RESTART(TM)
#endif

  RETURN

END SUBROUTINE WRITETRJ

SUBROUTINE WRITE_RESTART(TM)

  USE VAR

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: TM !< Current time step
  INTEGER :: I, J, L, KK
  REAL*4 :: TREAL

  TREAL = TM * DT * TIMESCALE * 1.0D+12 + INITIME

    OPEN ( 112 , FILE = 'restart', STATUS='replace')
    WRITE(112,*)TITLE
    WRITE(112,*)'Time	',TREAL, '	(ps)'
    WRITE(112,*)'***** Box Length(X, Y, Z in nanometres) **'
    WRITE(112,*)BOX(1), BOX(2), BOX(3)

9040 FORMAT (1X,I6,1X,I2,1X,I1,1X,3 (G21.14,1X))
9030 FORMAT (70 ('*'),/,10 ('*'),&
         ' Record for each atom is in the form:-            ',10 ('*'),/,10 ('*'), &
         ' Index Atom_Type No._of_bonds X Y Z (coords.in m) ',10 ('*'),/,10 ('*'), &
         ' Vx Vy Vz (in m/s) Indices_of_bonded_atoms        ',10 ('*'),/,70 ('*'))
    WRITE (112,9030)
    WRITE (112,*) 'num_of_molecules  ', NMOL
    L = 0
    DO I = 1,NMOL !120
       WRITE (112,*) NATM(I),'  Atoms_in_Molecule_No. ',I,name_mol(i)
       DO J = 1,NATM(I) !110
          L = L + 1
          WRITE (112,9040) L,ITYPE(L),NBONDS(L),SXYZ(1,L),SXYZ(2,L),SXYZ(3,L)
          SELECT CASE (NBONDS(L)) ! Format statement changes depending on number of bonds
          CASE(0)
             WRITE(112,9090)VXYZ(1,L)*VSCALE*1.d-3,VXYZ(2,L)*VSCALE*1.d-3,VXYZ(3,L)*VSCALE*1.d-3
9090         FORMAT (3(F16.10,1X))
          CASE(1)
             WRITE(112,9050)VXYZ(1,L)*VSCALE*1.d-3,VXYZ(2,L)*VSCALE*1.d-3,VXYZ(3,L)*VSCALE*1.d-3, & 
                  (JBOND(L,KK),KK=1,NBONDS(L))
9050         FORMAT (3(F16.10,1X), I5)
          CASE(2)
             WRITE(112,9060)VXYZ(1,L)*VSCALE*1.d-3,VXYZ(2,L)*VSCALE*1.d-3,VXYZ(3,L)*VSCALE*1.d-3, & 
                  (JBOND(L,KK),KK=1,NBONDS(L))
9060         FORMAT (3(F16.10,2X),2(I5,2X))
          CASE(3)
             WRITE(112,9070)VXYZ(1,L)*VSCALE*1.d-3,VXYZ(2,L)*VSCALE*1.d-3,VXYZ(3,L)*VSCALE*1.d-3, & 
                  (JBOND(L,KK),KK=1,NBONDS(L))
9070         FORMAT (3(F16.10,1X), 3(I5,2X))
          CASE(4)
             WRITE(112,9080)VXYZ(1,L)*VSCALE*1.d-3,VXYZ(2,L)*VSCALE*1.d-3,VXYZ(3,L)*VSCALE*1.d-3, & 
                  (JBOND(L,KK),KK=1,NBONDS(L))
9080         FORMAT (3(F16.10,1X), 4(I5,2X))
          END SELECT
          !        WRITE (112,*) VXYZ(1,L)*VSCALE*1.d-3,VXYZ(2,L)*VSCALE*1.d-3,VXYZ(3,L)*VSCALE*1.d-3, (JBOND(L,KK),KK=1,NBONDS(L))
       END DO! 110 	   	CONTINUE
    END DO!120 	CONTINUE

    CLOSE(112)

  RETURN
END SUBROUTINE WRITE_RESTART

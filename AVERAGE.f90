SUBROUTINE AVERAGE (TM)

  USE 	VAR
  IMPLICIT NONE
  INTEGER :: I, J, K, L, M, KK, TM
  REAL*8 :: TREAL, T
!       *******************************************************************
!       ********     STORING AVERAGE DATA AND RESTART FILE      ***********

 TREAL = TM * DT * TIMESCALE * 1.0D+12 + INITIME

 !	CALCULATE ROLLING AVERAGES
 IRAV = IRAV + 1
 IF (IRAV .GT. LIMRAV) IRAV = IRAV - LIMRAV
 NRAV = NRAV + 1
 IF (NRAV .GT. LIMRAV) NRAV = LIMRAV

 T = EK * MKTEMP
 STEK(IRAV) = EK
 STVBOND(IRAV) = VBOND
 STVANGLE(IRAV) = VANGLE
 STVTOR(IRAV) = VTOR
 STVNBOND(IRAV) = VNBOND
 STVOOP(IRAV) = VOOP
 STV(IRAV) = VBOND + VANGLE + VTOR + VNBOND + VOOP
 STE(IRAV) = STEK(IRAV) + STV(IRAV)
 STEMP(IRAV) = T

 STPT11(IRAV) = PT11
 STPT22(IRAV) = PT22
 STPT33(IRAV) = PT33
 STPT12(IRAV) = PT12
 STPT13(IRAV) = PT13
 STPT23(IRAV) = PT23
 STP(IRAV) = (PT11 + PT22 + PT33)/3.0D0

 SBOXX(IRAV) = BOXX
 SBOXY(IRAV) = BOXY
 SBOXZ(IRAV) = BOXZ
 SVOL(IRAV) = BOXX * BOXY * BOXZ
 SDENS(IRAV) = TOTMASS / SVOL(IRAV)

 REK = 0.0D0
 RVBOND = 0.0D0
 RVANGLE = 0.0D0
 RVTOR = 0.0D0
 RVOOP = 0.0D0
 RVNBOND = 0.0D0
 RV = 0.0D0
 RE = 0.0D0
 RTEMP = 0.0D0
 RPT11 = 0.0D0	
 RPT22 = 0.0D0
 RPT33 = 0.0D0
 RPT12 = 0.0D0
 RPT13 = 0.0D0
 RPT23 = 0.0D0
 RP = 0.0D0
 RBOXX = 0.0D0
 RBOXY = 0.0D0
 RBOXZ = 0.0D0
 RVOL = 0.0D0
 RDENS = 0.0D0

 DO I = 1, NRAV

    REK = REK + STEK(I)
    RVBOND = RVBOND + STVBOND(I)
    RVANGLE = RVANGLE + STVANGLE(I)
    RVTOR = RVTOR + STVTOR(I)
    RVOOP = RVOOP + STVOOP(I)
    RVNBOND = RVNBOND + STVNBOND(I)
    RV = RV + STV(I)
    RE = RE + STE(I)
    RTEMP = RTEMP + STEMP(I)
    RPT11 = RPT11 + STPT11(I)
    RPT22 = RPT22 + STPT22(I)
    RPT33 = RPT33 + STPT33(I)
    RPT12 = RPT12 + STPT12(I)
    RPT13 = RPT13 + STPT13(I)
    RPT23 = RPT23 + STPT23(I)
    RP = RP + STP(I)
    RBOXX = RBOXX + SBOXX(I)
    RBOXY = RBOXY + SBOXY(I)
    RBOXZ = RBOXZ + SBOXZ(I)
    RVOL = RVOL + SVOL(I)
    RDENS = RDENS + SDENS(I)
 END DO

 REK = REK / NRAV
 RVBOND = RVBOND / NRAV
 RVANGLE = RVANGLE / NRAV
 RVTOR = RVTOR / NRAV
 RVOOP = RVOOP/ NRAV
 RVNBOND = RVNBOND / NRAV
 RV = RV / NRAV
 RE = RE / NRAV
 RTEMP = RTEMP / NRAV
 RPT11 = RPT11 / NRAV
 RPT22 = RPT22 / NRAV
 RPT33 = RPT33 / NRAV
 RPT12 = RPT12 / NRAV
 RPT13 = RPT13 / NRAV
 RPT23 = RPT23 / NRAV
 RP = RP / NRAV
 RBOXX = RBOXX / NRAV
 RBOXY = RBOXY / NRAV
 RBOXZ = RBOXZ / NRAV
 RVOL = RVOL / NRAV
 RDENS = RDENS / NRAV

 !	WRITING THE RESTART FILE
	
 OPEN ( 112 , FILE = 'restart', STATUS='replace')
 WRITE(112,*)TITLE
 WRITE(112,*)'Time	',TREAL, '	(ps)'
 WRITE(112,*)'***** Box Length(X, Y, Z in nanometres) **'
 WRITE(112,*)BOXX, BOXY, BOXZ

9040 FORMAT (1X,I6,1X,I2,1X,I1,1X,3 (G21.14,1X))
9030 FORMAT (70 ('*'),/,10 ('*'),&
          ' Record for each atom is in the form:-            ', &
          10 ('*'),/,10 ('*'), &
          ' Index Atom_Type No._of_bonds X Y Z (coords.in m) ', &
          10 ('*'),/,10 ('*'), &
          ' Vx Vy Vz (in m/s) Indices_of_bonded_atoms        ', &
          10 ('*'),/,70 ('*'))
 WRITE (112,9030)
 WRITE (112,*) 'num_of_molecules  ', NMOL
 L = 0
 DO I = 1,NMOL !120
    WRITE (112,*) NATM(I),'  Atoms_in_Molecule_No. ',I,name_mol(i)
    DO J = 1,NATM(I) !110
       L = L + 1
       WRITE (112,9040) L,ITYPE(L),NBONDS(L),SX(L),SY(L),SZ(L)
       SELECT CASE (NBONDS(L))
       CASE(1)
          WRITE(112,9050)VX(L)*VSCALE*1.d-3,VY(L)*VSCALE*1.d-3,VZ(L)*VSCALE*1.d-3,(JBOND(L,KK),KK=1,NBONDS(L))
9050      FORMAT (3(F16.10,1X), I5)
       CASE(2)
          WRITE(112,9060)VX(L)*VSCALE*1.d-3,VY(L)*VSCALE*1.d-3,VZ(L)*VSCALE*1.d-3,(JBOND(L,KK),KK=1,NBONDS(L))
9060      FORMAT (3(F16.10,2X),2(I5,2X))
       CASE(3)
          WRITE(112,9070)VX(L)*VSCALE*1.d-3,VY(L)*VSCALE*1.d-3,VZ(L)*VSCALE*1.d-3,(JBOND(L,KK),KK=1,NBONDS(L))
9070      FORMAT (3(F16.10,1X), 3(I5,2X))
       CASE(4)
          WRITE(112,9080)VX(L)*VSCALE*1.d-3,VY(L)*VSCALE*1.d-3,VZ(L)*VSCALE*1.d-3,(JBOND(L,KK),KK=1,NBONDS(L))
9080      FORMAT (3(F16.10,1X), 4(I5,2X))
       END SELECT
!        WRITE (112,*) VX(L)*VSCALE*1.d-3,VY(L)*VSCALE*1.d-3,VZ(L)*VSCALE*1.d-3, (JBOND(L,KK),KK=1,NBONDS(L))
    END DO! 110 	   	CONTINUE
 END DO!120 	CONTINUE

 CLOSE (112)

 RETURN
END SUBROUTINE AVERAGE

!	*********************************************************************************************

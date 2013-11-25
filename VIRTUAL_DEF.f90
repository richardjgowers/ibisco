SUBROUTINE VIRTUAL_DEF()

USE VAR

IMPLICIT NONE

INTEGER :: I, J, K, TI, POS
REAL(KIND=RKIND) :: SUMTOTX, SUMTOTY, SUMTOTZ
REAL(KIND=RKIND) :: SPX, SPY, SPZ

!Calculate centres of virtual sites
DO I=1,NVIRTA
   POS = I + NATOMS

   !Assign type_labels to virtual sites
   TYPE_LABEL(POS) = 3

   IF (VIRT_CENTER(I) .NE. 0)THEN !If using a functional site
      J = VIRT_CENTER(I)

      RX(POS) = SX(J)
      RY(POS) = SY(J)
      RZ(POS) = SZ(J)

      SX(POS) = SX(J)
      SY(POS) = SY(J)
      SZ(POS) = SZ(J)
   ELSE !Else using a COM
      TI = VITYPE(I)
      SUMTOTX = 0.0D0
      SUMTOTY = 0.0D0
      SUMTOTZ = 0.0D0
      !Finds COM
      DO J=1,VIRT_NUMATOMS(TI)
         K = VIRT_ATM_IND(I,J)
         SUMTOTX = SUMTOTX + MASS(ITYPE(K))*SX(K)
         SUMTOTY = SUMTOTY + MASS(ITYPE(K))*SY(K)
         SUMTOTZ = SUMTOTZ + MASS(ITYPE(K))*SZ(K)
      END DO

      RX(POS) = SUMTOTX*VIRT_INVMASS(TI)
      RY(POS) = SUMTOTY*VIRT_INVMASS(TI)
      RZ(POS) = SUMTOTZ*VIRT_INVMASS(TI)

      SX(POS) = RX(POS)
      SY(POS) = RY(POS)
      SZ(POS) = RZ(POS)
   END IF

   IF ((RX(POS)> BOXX2).OR.(RX(POS) < -BOXX2) &
        .OR.(RY(POS)> BOXY2) .OR.(RY(POS)< -BOXY2) &
        .OR.(RZ(POS)> BOXZ2).OR.(RZ(POS)<-BOXZ2)) THEN

      RX(POS) = RX(POS)/BOXX2
      RY(POS) = RY(POS)/BOXY2
      RZ(POS) = RZ(POS)/BOXZ2

      SPX = RX(POS) - 2.0*INT(RX(POS)/2.0)
      SPY = RY(POS) - 2.0*INT(RY(POS)/2.0)
      SPZ = RZ(POS) - 2.0*INT(RZ(POS)/2.0)

      RX(POS) = (SPX - 2.0*INT(SPX))*BOXX2
      RY(POS) = (SPY - 2.0*INT(SPY))*BOXY2
      RZ(POS) = (SPZ - 2.0*INT(SPZ))*BOXZ2
   END IF


END DO



RETURN 
END SUBROUTINE VIRTUAL_DEF

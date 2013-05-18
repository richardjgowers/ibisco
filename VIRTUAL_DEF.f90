SUBROUTINE VIRTUAL_DEF()

USE VAR

IMPLICIT NONE

INTEGER :: I, J, K, TI
REAL(KIND=RKIND) :: SUMTOTX, SUMTOTY, SUMTOTZ
REAL(KIND=RKIND) :: RPX, RPY, RPZ, SPX, SPY, SPZ
REAL(KIND=RKIND) :: SUMTOTBMASS

!Calculate centres of virtual sites
DO I=1,NVIRTA
   IF (VIRT_CENTER(I) .NE. 0)THEN !If using a functional site
      J = VIRT_CENTER(I)
      VIRTRX(I) = RX(J)
      VIRTRY(I) = RY(J)
      VIRTRZ(I) = RZ(J)
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

      VIRTRX(I) = SUMTOTX*VIRT_INVMASS(TI)
      VIRTRY(I) = SUMTOTY*VIRT_INVMASS(TI)
      VIRTRZ(I) = SUMTOTZ*VIRT_INVMASS(TI)

      !Applies PBC

      IF ((VIRTRX(I) .gt. BOXX2).OR.(VIRTRX(I) .lt. -BOXX2) &
           .OR.(VIRTRY(I) .gt. BOXY2).OR.(VIRTRY(I) .lt. -BOXY2)&
           .OR.(VIRTRZ(I) .gt. BOXZ2).OR.(VIRTRZ(I) .lt. -BOXZ2)) THEN

         RPX = VIRTRX(I)/BOXX2
         RPY = VIRTRY(I)/BOXY2
         RPZ = VIRTRZ(I)/BOXZ2

         SPX = RPX - 2.0D0 * INT(RPX/2.0)
         SPY = RPY - 2.0D0 * INT(RPY/2.0)
         SPZ = RPZ - 2.0D0 * INT(RPZ/2.0)

         RPX = (SPX - 2.0D0 * INT(SPX))*BOXX2
         RPY = (SPY - 2.0D0 * INT(SPY))*BOXY2
         RPZ = (SPZ - 2.0D0 * INT(SPZ))*BOXZ2

         VIRTRX(I) = RPX
         VIRTRY(I) = RPY
         VIRTRZ(I) = RPZ 
      END IF
   END IF
END DO

RETURN 
END SUBROUTINE VIRTUAL_DEF

!> @file
!> @brief Defines the positions of virtual sites in the system
!!
!> @details Virtuals sites can take their position as either COM of the atoms or defined by the 
!! position of another atom
!> @author Richard J Gowers

SUBROUTINE VIRTUAL_DEF()

USE VAR

IMPLICIT NONE

INTEGER :: I, J, K, TI, POS
REAL*4, DIMENSION(3) :: SUMTOTXYZ
REAL*4, DIMENSION(3) :: SPXYZ

!Calculate centres of virtual sites
DO I=1,NVIRTA
   POS = I + NATOMS

   !Assign type_labels to virtual sites
   TYPE_LABEL(POS) = 3

   IF (VIRT_CENTER(I) .NE. 0)THEN !If using a functional site
      J = VIRT_CENTER(I)

      RXYZ(1,POS) = SXYZ(1,J)
      RXYZ(2,POS) = SXYZ(2,J)
      RXYZ(3,POS) = SXYZ(3,J)

      SXYZ(1,POS) = SXYZ(1,J)
      SXYZ(2,POS) = SXYZ(2,J)
      SXYZ(3,POS) = SXYZ(3,J)
   ELSE !Else using a COM
      TI = VITYPE(I)
      SUMTOTXYZ(1) = 0.0D0
      !Finds COM
      DO J=1,VIRT_NUMATOMS(TI)
         K = VIRT_ATM_IND(I,J)
         SUMTOTXYZ(1) = SUMTOTXYZ(1) + MASS(ITYPE(K))*SXYZ(1,K)
         SUMTOTXYZ(2) = SUMTOTXYZ(2) + MASS(ITYPE(K))*SXYZ(2,K)
         SUMTOTXYZ(3) = SUMTOTXYZ(3) + MASS(ITYPE(K))*SXYZ(3,K)
      END DO

      RXYZ(1,POS) = SUMTOTXYZ(1)*VIRT_INVMASS(TI)
      RXYZ(2,POS) = SUMTOTXYZ(2)*VIRT_INVMASS(TI)
      RXYZ(3,POS) = SUMTOTXYZ(3)*VIRT_INVMASS(TI)

      SXYZ(1,POS) = RXYZ(1,POS)
      SXYZ(2,POS) = RXYZ(2,POS)
      SXYZ(3,POS) = RXYZ(3,POS)
   END IF

   IF ((RXYZ(1,POS)> BOXX2).OR.(RXYZ(1,POS) < -BOXX2) &
        .OR.(RXYZ(2,POS)> BOXY2) .OR.(RXYZ(2,POS)< -BOXY2) &
        .OR.(RXYZ(3,POS)> BOXZ2).OR.(RXYZ(3,POS)<-BOXZ2)) THEN

      RXYZ(1,POS) = RXYZ(1,POS)/BOXX2
      RXYZ(2,POS) = RXYZ(2,POS)/BOXY2
      RXYZ(3,POS) = RXYZ(3,POS)/BOXZ2

      SPXYZ(1) = RXYZ(1,POS) - 2.0*INT(RXYZ(1,POS)/2.0)
      SPXYZ(2) = RXYZ(2,POS) - 2.0*INT(RXYZ(2,POS)/2.0)
      SPXYZ(3) = RXYZ(3,POS) - 2.0*INT(RXYZ(3,POS)/2.0)

      RXYZ(1,POS) = (SPXYZ(1) - 2.0*INT(SPXYZ(1)))*BOXX2
      RXYZ(2,POS) = (SPXYZ(2) - 2.0*INT(SPXYZ(2)))*BOXY2
      RXYZ(3,POS) = (SPXYZ(3) - 2.0*INT(SPXYZ(3)))*BOXZ2
   END IF


END DO



RETURN 
END SUBROUTINE VIRTUAL_DEF

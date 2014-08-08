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

        RXYZ(:,POS) = SXYZ(:,J)
        SXYZ(:,POS) = SXYZ(:,J)
     ELSE !Else using a COM
        TI = VITYPE(I)
        SUMTOTXYZ(1) = 0.0
        !Finds COM
        DO J=1,VIRT_NUMATOMS(TI)
           K = VIRT_ATM_IND(I,J)
           SUMTOTXYZ(:) = SUMTOTXYZ(:) + MASS(ITYPE(K)) * SXYZ(:,K)
        END DO

        RXYZ(:,POS) = SUMTOTXYZ(:) * VIRT_INVMASS(TI)
        SXYZ(:,POS) = RXYZ(:,POS)
     END IF

     IF (ANY((SXYZ(:,POS) > BOX2(:)) .OR. (SXYZ(:,POS) < -BOX2(:)))) THEN
        RXYZ(:,POS) = RXYZ(:,POS) / BOX2(:)
        SPXYZ(:) = RXYZ(:,POS) - 2.0 * INT(RXYZ(:,POS) / 2.0)
        RXYZ(:,POS) = (SPXYZ(:) - 2.0 * INT(SPXYZ(:))) * BOX2(:)
     END IF
  END DO

  RETURN 
END SUBROUTINE VIRTUAL_DEF

SUBROUTINE RDVIRTUAL()

USE VAR

IMPLICIT NONE

INTEGER :: iost, I, J, K, NUMATOM, TI, TJ, A
INTEGER, PARAMETER :: MAX_ATOMS = 15 !Temp max number of atoms within a VS
INTEGER :: ACTUAL_MAX
REAL(KIND=RKIND), PARAMETER :: MASSTOL = 0.0001
REAL*8 :: SUMTOTBMASS
INTEGER, POINTER :: INDX_ATM(:,:) !Temp array for filling in index of atoms before max VIRT_NUMATOMS is known

OPEN(UNIT=10,FILE='virtual',status='old',IOSTAT=iost)

IF(iost .ne. 0) THEN
   WRITE(*,*) '*** FATAL ERROR! File virtual site does not exist ****'
   WRITE(1,*) '*** FATAL ERROR! File virtual site does not exist ****'
   ISTOP=1
   RETURN
END IF

READ(10,*)
READ(10,*)
READ(10,*) NVIRTA
READ(10,*)
READ(10,*) 

!Variables
!NVIRTA Number of virtual sites
ALLOCATE(VIRTRX(NVIRTA),VIRTRY(NVIRTA),VIRTRZ(NVIRTA)) !Position of the virtual site, either a COM or an atom
ALLOCATE(VIRT_NUMATOMS(NTYPE)) !Number of atoms in different virtual site types
ALLOCATE(VIRT_MASS(NTYPE), VIRT_INVMASS(NTYPE)) !Mass of VS
ALLOCATE(VIRT_MASSCOEFF(NTYPE,NTYPE)) !Masscoefficient for atom in VS. Usage VIRT_MASSCOEFF(virt type, atom type)
ALLOCATE(VITYPE(NVIRTA)) !Type of each virtual site.  Types are the same as in interaction file
ALLOCATE(VIRT_CENTER(NVIRTA)) !Center is 0 for COM VS or index of atom
ALLOCATE(VIRT_VS_IND(NATOMS)) !Returns the VS that an atom belongs to
ALLOCATE(INDX_ATM(NVIRTA,MAX_ATOMS)) !Temp array for reading info
!ALLOCATE(VIRT_ATM_IND(NVIRTA,ACTUAL_MAX)) !Returns the index of atoms in a VS
!The number of atoms and mass of a virtual site is defined by the type of virtual site
!ie all virtual sites have the same properties

!The mass coefficient of atoms within a site can change as atoms can appear in a different order within VS
!Eg -(C-C-N-C)- or -(C-N-C-C)-
!The type of an atom within a virtual site will always have the same mass coefficient

ALLOCATE(VIRT_POINT_SEC(NVIRTA+1))
ALLOCATE(VCELL(NVIRTA))
ALLOCATE(VLCLIST(NVIRTA))

indx_atm = 0
VIRT_VS_IND = 0

!Read virtual site information
VIRT_NUMATOMS = 0
DO I=1,NVIRTA
   READ(10,*,iostat=iost) K, VITYPE(I), NUMATOM, VIRT_CENTER(I)
   ITYPE(NATOMS + I) = VITYPE(I)
   IF(iost .ne. 0) THEN
      WRITE(*,*) 'ERROR READING VIRTUAL, TOO SHORT'
      WRITE(1,*) 'ERROR READING VIRTUAL, TOO SHORT'
      ISTOP = 1
      RETURN
   END IF

   IF(VIRT_CENTER(I) .ne. 0) VIRT_VS_IND(VIRT_CENTER(I)) = I

   IF(NUMATOM .gt. MAX_ATOMS) THEN
      WRITE(*,*) 'Number of atoms in VS exceeds maximum, check virtual file'
      ISTOP = 1
      RETURN
   END IF

   READ(10,*) (INDX_ATM(I,J), J=1,NUMATOM)

   IF(VIRT_NUMATOMS(VITYPE(I)) .eq. 0) THEN !If number of atoms in VITYPE not assigned
      VIRT_NUMATOMS(VITYPE(I)) = NUMATOM
   ELSE IF(VIRT_NUMATOMS(VITYPE(I)) .ne. NUMATOM) THEN !If numatoms doesn't agree with old value
      WRITE(*,*) 'Disagreement in number of atoms in a VS type, check virtual file'
      ISTOP = 1
      RETURN
   END IF
END DO

!Transfer atom index info into smallest array possible
ACTUAL_MAX = 0
ACTUAL_MAX = MAXVAL(VIRT_NUMATOMS) !Find the maximum number of atoms in a VS
ALLOCATE(VIRT_ATM_IND(NVIRTA,ACTUAL_MAX))
VIRT_ATM_IND = 0

DO I=1,NVIRTA
   DO J=1,VIRT_NUMATOMS(VITYPE(I))
      K = INDX_ATM(I,J)

      VIRT_ATM_IND(I,J) = K
!      VIRT_VS_IND(K) = I
   END DO
END DO

!Calculate VS mass information
!VS mass is taken as sum of atoms within and NOT the mass of the bead it represents

VIRT_MASS = 0.0D0 !Mass of different types of virtual sites
VIRT_INVMASS = 0.0D0 !Reciprocal of mass to reduce divide operations
VIRT_MASSCOEFF = 0.0D0 !VIRT_MASSCOEFF(virtual type, atom type) returns the mass coefficient for the atom

DO I=1,NVIRTA
   TI = VITYPE(I)
   sumtotBmass = 0.0D0
   DO J=1,VIRT_NUMATOMS(TI)
      TJ = ITYPE(VIRT_ATM_IND(I,J))
      sumtotBmass = sumtotBmass + MASS(TJ)
   END DO

   IF(VIRT_MASS(TI) .eq. 0) THEN
      VIRT_MASS(TI) = sumtotBmass
      VIRT_INVMASS(TI) = 1.0D0 / sumtotBmass
   ELSE IF(ABS(VIRT_MASS(TI) - sumtotBmass) .gt. MASSTOL) THEN
      WRITE(*,*) 'Disagreement in VS mass, check virtual file'
      WRITE(*,*) sumtotBmass, VIRT_MASS(TI)
      ISTOP = 1
      RETURN
   END IF

   DO J=1,VIRT_NUMATOMS(TI)
      TJ = ITYPE(VIRT_ATM_IND(I,J))
      VIRT_MASSCOEFF(TI,TJ) = MASS(TJ)*VIRT_INVMASS(TI)
   END DO
END DO


DO I=1,NVIRTA
   TI = VITYPE(I)
   DO A=1,VIRT_NUMATOMS(TI)
      J = VIRT_ATM_IND(I,A)
      TJ = ITYPE(J)
   END DO
END DO

RETURN
END SUBROUTINE RDVIRTUAL

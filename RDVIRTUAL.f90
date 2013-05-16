SUBROUTINE RDVIRTUAL()

USE VAR

IMPLICIT NONE

INTEGER :: iost, I, J, K
INTEGER, PARAMETER :: NBEAD = 15 !Max number of atoms within a VS
REAL(KIND=RKIND) :: SUMTOTX, SUMTOTY, SUMTOTZ
REAL(KIND=RKIND) :: RPX, RPY, RPZ, SPX, SPY, SPZ


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

NB = 2

ALLOCATE(VIRTRX(NVIRTA))
ALLOCATE(VIRTRY(NVIRTA))
ALLOCATE(VIRTRZ(NVIRTA))
ALLOCATE(INDEX_VSITE(NVIRTA))
ALLOCATE(TYPE_BEAD(NB))
ALLOCATE(NUMATOM(NB))
ALLOCATE(TOTBMASS(NVIRTA))
ALLOCATE(INVTOTBMASS(NVIRTA))
ALLOCATE(VITYPE(NVIRTA))
allocate(init_numbcomp(nvirta))
ALLOCATE(masscoeff(NVIRTA,NBEAD))
ALLOCATE(VIRT_POINT_SEC(NVIRTA+1))
ALLOCATE(VCELL(NVIRTA))
ALLOCATE(VLCLIST(NVIRTA))
ALLOCATE(INDX_ATM(NVIRTA,NBEAD))
ALLOCATE(virtual_center(NATOMS))
indx_atm = 0
virtual_center = 0

TOTBMASS = 0
INVTOTBMASS = 0
masscoeff=0

!Read virtual site information

DO I=1,NVIRTA
   read(10,*) K, VITYPE(I), INIT_NUMBCOMP(I), INDEX_VSITE(I)
   read(10,*,IOSTAT=iost) (INDX_ATM(I,J), J=1,INIT_NUMBCOMP(I))

   virtual_center(INDEX_VSITE(I)) = I

   IF(iost .ne. 0) THEN
      WRITE(*,*) 'ERROR READING VIRTUAL, TOO SHORT'
      WRITE(1,*) 'ERROR READING VIRTUAL, TOO SHORT'
      ISTOP = 1
      RETURN
   END IF
END DO

DO I=1,NVIRTA
   TOTBMASS(I) = mass(vitype(I))
   INVTOTBMASS(I) = 1/mass(vitype(I))
END DO


!Calculate centre coordinates of virtual sites
DO I=1,NVIRTA
   IF (INDEX_VSITE(I) .NE. 0)THEN !If using a functional site
      J = INDEX_VSITE(I)
      VIRTRX(I) = RX(J)
      VIRTRY(I) = RY(J)
      VIRTRZ(I) = RZ(J)
   ELSE !Else using a COM
      SUMTOTX = 0.0
      SUMTOTY = 0.0
      SUMTOTZ = 0.0
      !Finds COM
      DO J=1,INIT_NUMBCOMP(I)
         K = INDX_ATM(I,J)
         SUMTOTX = SUMTOTX + MASS(ITYPE(K))*SX(K)
         SUMTOTY = SUMTOTY + MASS(ITYPE(K))*SY(K)
         SUMTOTZ = SUMTOTZ + MASS(ITYPE(K))*SZ(K)
      END DO

      VIRTRX(I) = SUMTOTX/MASS(VITYPE(I)) 
      VIRTRY(I) = SUMTOTY/MASS(VITYPE(I)) 
      VIRTRZ(I) = SUMTOTZ/MASS(VITYPE(I)) 

      !Applies PBC

      IF ((VIRTRX(I) .gt. BOXX2).OR.(VIRTRX(I) .lt. -BOXX2) &
           .OR.(VIRTRY(I) .gt. BOXY2).OR.(VIRTRY(I) .lt. -BOXY2)&
           .OR.(VIRTRZ(I) .gt. BOXZ2).OR.(VIRTRZ(I) .lt. -BOXZ2)) THEN

         RPX = VIRTRX(I)/BOXX2
         RPY = VIRTRY(I)/BOXY2
         RPZ = VIRTRZ(I)/BOXZ2

         SPX = RPX - 2.0*INT(RPX/2.0)
         SPY = RPY - 2.0*INT(RPY/2.0)
         SPZ = RPZ - 2.0*INT(RPZ/2.0)

         RPX = (SPX - 2.0*INT(SPX))*BOXX2
         RPY = (SPY - 2.0*INT(SPY))*BOXY2
         RPZ = (SPZ - 2.0*INT(SPZ))*BOXZ2

         VIRTRX(I) = RPX
         VIRTRY(I) = RPY
         VIRTRZ(I) = RPZ 
      END IF
   END IF
END DO




RETURN
END SUBROUTINE RDVIRTUAL

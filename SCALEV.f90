
	SUBROUTINE SCALEV ()
	USE VAR
	
	IMPLICIT	NONE
        REAL*8        FACTOR, T
	INTEGER		I
!       *******************************************************************

	EK = 0.0D0
        DO 5 I = 1, NATOMS

           EK = EK + MASS(ITYPE(I))*(VX(I)**2.0D0 + VY(I)**2.0D0 + VZ(I)**2.0D0)

5       CONTINUE

        EK = 0.5D0 * EK
	T = EK *  MKTEMP
	FACTOR = TIN / T

	DO 100 I = 1, NATOMS

		CM = BEADMASS(I)
		VX(I) = SQRT(FACTOR) * VX(I)
		VY(I) = SQRT(FACTOR) * VY(I)
		VZ(I) = SQRT(FACTOR) * VZ(I)
100	CONTINUE

	RETURN
	END



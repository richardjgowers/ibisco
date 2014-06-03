!> @file
!> @brief Generate random velocities
!> @details I have no idea if this works vOv
!!
!! called from main.f90

       SUBROUTINE COMVEL ()

	USE VAR
	IMPLICIT	NONE

     	REAL*8        SUMX, SUMY, SUMZ
	REAL*8        Q, TETA, RO, X, Y
      	INTEGER       I
!       *******************************************************************

	I = 1
        DO WHILE(I <= NATOMS)

		CALL RANDOM_NUMBER (X)
		RO = -LOG(X)
		CALL RANDOM_NUMBER(Y)
		TETA = Y * 8.0D0 * ATAN(1.0D0)
		Q =( 2.0D0 * RO )** 0.50D0 * COS(TETA)
		IF ( Q ** 2.0D0 <= 2.0D0  ) THEN
			VX(I) = Q
			I = I + 1
		END IF
	END DO
	I = 1
	DO WHILE(I <= NATOMS)

		CALL RANDOM_NUMBER (X)
		RO = -LOG(X)
		CALL RANDOM_NUMBER(Y)
		TETA = Y * 8.0D0 * ATAN(1.0D0)
		Q = ( 2.0D0 * RO )** 0.50D0 * COS(TETA)
		IF ( Q ** 2.0D0 <= 2.0D0  ) THEN
			VY(I) = Q
			I = I + 1
		END IF
	END DO
	I = 1
	DO WHILE(I <= NATOMS)

		CALL RANDOM_NUMBER (X)
		RO = -LOG(X)
		CALL RANDOM_NUMBER(Y)
		TETA = Y * 8.0D0 * ATAN(1.0D0)
		Q = ( 2.0D0 * RO )** 0.50D0 * COS(TETA)
		IF ( Q ** 2.0D0 <= 2.0D0  ) THEN
			VZ(I) = Q
			I = I + 1
		END IF
	END DO

        SUMX = 0.0D0
        SUMY = 0.0D0
        SUMZ = 0.0D0

        DO 200 I = 1, NATOMS

           SUMX = SUMX + MASS(ITYPE(I))*VX(I)
           SUMY = SUMY + MASS(ITYPE(I))*VY(I)
           SUMZ = SUMZ + MASS(ITYPE(I))*VZ(I)
	   
200     CONTINUE

        SUMX = SUMX / REAL ( NATOMS )
        SUMY = SUMY / REAL ( NATOMS )
        SUMZ = SUMZ / REAL ( NATOMS )

        DO 300 I = 1, NATOMS

           VX(I) = VX(I) - SUMX/MASS(ITYPE(I))
           VY(I) = VY(I) - SUMY/MASS(ITYPE(I))
           VZ(I) = VZ(I) - SUMZ/MASS(ITYPE(I))

300     CONTINUE

        RETURN
        END

!	*********************************************************************************************

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
		TETA = Y * 8.0 * ATAN(1.0)
		Q =( 2.0 * RO )** 0.50 * COS(TETA)
		IF ( Q ** 2.0 <= 2.0  ) THEN
			VXYZ(1,I) = Q
			I = I + 1
		END IF
	END DO
	I = 1
	DO WHILE(I <= NATOMS)

		CALL RANDOM_NUMBER (X)
		RO = -LOG(X)
		CALL RANDOM_NUMBER(Y)
		TETA = Y * 8.0 * ATAN(1.0)
		Q = ( 2.0 * RO )** 0.50 * COS(TETA)
		IF ( Q ** 2.0 <= 2.0  ) THEN
			VXYZ(2,I) = Q
			I = I + 1
		END IF
	END DO
	I = 1
	DO WHILE(I <= NATOMS)

		CALL RANDOM_NUMBER (X)
		RO = -LOG(X)
		CALL RANDOM_NUMBER(Y)
		TETA = Y * 8.0 * ATAN(1.0)
		Q = ( 2.0 * RO )** 0.50 * COS(TETA)
		IF ( Q ** 2.0 <= 2.0  ) THEN
			VXYZ(3,I) = Q
			I = I + 1
		END IF
	END DO

        SUMX = 0.0
        SUMY = 0.0
        SUMZ = 0.0

        DO 200 I = 1, NATOMS

           SUMX = SUMX + MASS(ITYPE(I))*VXYZ(1,I)
           SUMY = SUMY + MASS(ITYPE(I))*VXYZ(2,I)
           SUMZ = SUMZ + MASS(ITYPE(I))*VXYZ(3,I)
	   
200     CONTINUE

        SUMX = SUMX / REAL ( NATOMS )
        SUMY = SUMY / REAL ( NATOMS )
        SUMZ = SUMZ / REAL ( NATOMS )

        DO 300 I = 1, NATOMS

           VXYZ(1,I) = VXYZ(1,I) - SUMX/MASS(ITYPE(I))
           VXYZ(2,I) = VXYZ(2,I) - SUMY/MASS(ITYPE(I))
           VXYZ(3,I) = VXYZ(3,I) - SUMZ/MASS(ITYPE(I))

300     CONTINUE

        RETURN
        END

!	*********************************************************************************************

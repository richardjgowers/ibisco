!	*****************************************************************************************
      	SUBROUTINE RDGAUSSIAN()
	USE VAR
	USE MODULEPARSING
        IMPLICIT NONE

        INTEGER  :: NGAUBND, NGAUANG, TYPINIT
        INTEGER :: I, J, k, l 
!    *******************************************************************	
	OPEN (11, IOSTAT=IOS, FILE='gaussian', STATUS='OLD')

	IF (IOS.NE.0) THEN
	        WRITE (1,*)	&
		   ' **** FATAL ERROR! File gaussian does not exist ****'
		WRITE (*,*)	&
		   ' **** FATAL ERROR! File gaussian does not exist ****'
        	ISTOP=1
         	RETURN
      	END IF

	ALLOCATE(GB(MAXINPUT, NBTYPE))
	ALLOCATE(GA(MAXINPUT, NATYPE))
	ALLOCATE(LB(NBTYPE+1))
	ALLOCATE(LA(NATYPE+1))
	NGAUBND = 0
	NGAUANG = 0
	GB = 0.0
        GA = 0.0
  
	DO I = 1, NBTYPE + NATYPE
	
	    READ (11, '(A80)') LINE
	    CALL PARSE ()

	    IF (STRNGS(2) == 'gaussians_in_bond_type') THEN
		READ (STRNGS(1),*) J
		
		READ (STRNGS(3),*) K
		DO L = 1, J
		  READ (11,*) GB(NGAUBND+(L-1)*3+1,K), GB(NGAUBND+(L-1)*3+2,K), GB(NGAUBND+(L-1)*3+3,K)
		END DO
		LB(K) = NGAUBND + 1
		NGAUBND = NGAUBND + J*3
		
	    END IF

	    IF (STRNGS(2) == 'gaussians_in_angle_type') THEN
		READ (STRNGS(1),*) J
		
		READ (STRNGS(3),*) K
		DO L = 1, J
		  READ (11,*) GA(NGAUANG+(L-1)*3+1,K), GA(NGAUANG+(L-1)*3+2,K), GA(NGAUANG+(L-1)*3+3,K)
		END DO
		LA(K) = NGAUANG + 1
		NGAUANG = NGAUANG + J*3
		
	    END IF

	END DO

!**** Read gaussian file ****
!        READ (11,*)
!        READ (11,*) NGAUBND, NGAUANG
!        TYPINIT=0

!        DO 30  I=1, NGAUBND
!        READ (11,*) TMP, J
!        GB(I,J) = TMP
!        IF (J.NE.TYPINIT) THEN
!        LB(J) = I
!        TYPINIT = J
!        ENDIF
!30      CONTINUE

        LB(NBTYPE+1) = NGAUBND + 1 
!        TYPINIT=0
 
!        DO 40 I=1,NGAUANG
!        READ (11,*) TMP, J
!        GA(I,J)=TMP
!        IF (J.NE.TYPINIT) THEN
!        LA(J) = I
!        TYPINIT = J
!        ENDIF
! 40     CONTINUE

        LA(NATYPE+1) = NGAUANG + 1
        CLOSE (11)

	RETURN
        END

!	*********************************************************************************************

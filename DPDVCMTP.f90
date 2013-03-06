	SUBROUTINE DPDVCMTP ()

	USE VAR
	IMPLICIT	NONE
        INTEGER       I, J
	REAL*8	      VXI, VYI, VZI
	REAL*8	      PT11K, PT22K, PT33K
	REAL*8	      PT12K, PT23K, PT13K
!       *******************************************************************
	EK = 0.0D0

	PT11K = 0.0D0
	PT22K = 0.0D0
	PT33K = 0.0D0

	PT12K = 0.0D0
	PT23K = 0.0D0
	PT13K = 0.0D0

        DO 100 I = 1, NATOMS
         
           VXI = VX(I)     
           VYI = VY(I)
           VZI = VZ(I)  
           CM = BEADMASS(I)

	   PT11K = PT11K + CM *VXI ** 2.0
	   PT22K = PT22K + CM *VYI ** 2.0
	   PT33K = PT33K + CM *VZI ** 2.0

	   PT12K = PT12K + CM * VXI * VYI
	   PT13K = PT13K + CM * VXI * VZI
           PT23K = PT23K + CM * VYI * VZI

100     CONTINUE

	PT11 = PT11 + PT11K 
	PT22 = PT22 + PT22K
	PT33 = PT33 + PT33K

	PT12 = PT12 + PT12K 
	PT13 = PT13 + PT13K
	PT23 = PT23 + PT23K

	EK = 0.5D0 * (PT11K + PT22K + PT33K)

	TEMP = EK * MKTEMP
        RETURN
        END
!	*********************************************************************************************

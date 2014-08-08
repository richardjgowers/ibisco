  DO A = 1, NOOPIJKL(I)

     J = JOOPIJKL(I, A)
     K = KOOPIJKL(I, A)
     L = LOOPIJKL(I, A)

     TJ = ITYPE(J)
     TK = ITYPE(K)
     TL = ITYPE(L)
     TOIJKL = IOOPT(TI,TJ,TK,TL)

     RXYZ_IJ(1) = RXYZ_I(1) - SXYZ(1,J) 
     RXYZ_IJ(2) = RXYZ_I(2) - SXYZ(2,J) 
     RXYZ_IJ(3) = RXYZ_I(3) - SXYZ(3,J) 

     RXYZ_JK(1) = SXYZ(1,K) - SXYZ(1,J) 
     RXYZ_JK(2) = SXYZ(2,K) - SXYZ(2,J) 
     RXYZ_JK(3) = SXYZ(3,K) - SXYZ(3,J) 

     RXYZ_KL(1) = SXYZ(1,K) - SXYZ(1,L) 
     RXYZ_KL(2) = SXYZ(2,K) - SXYZ(2,L) 
     RXYZ_KL(3) = SXYZ(3,K) - SXYZ(3,L)

     RIJ = SQRT(RXYZ_IJ(1) ** 2.0 + RXYZ_IJ(2) ** 2.0 + RXYZ_IJ(3) ** 2.0)
     RJK = SQRT(RXYZ_JK(1) ** 2.0 + RXYZ_JK(2) ** 2.0 + RXYZ_JK(3) ** 2.0)
     RKL = SQRT(RXYZ_KL(1) ** 2.0 + RXYZ_KL(2) ** 2.0 + RXYZ_KL(3) ** 2.0)

     !	****** vector M cross product IJxJK
     RXYZ_M(1) = RXYZ_IJ(2)*RXYZ_JK(3) - RXYZ_IJ(3)*RXYZ_JK(2)
     RXYZ_M(2) =-RXYZ_IJ(1)*RXYZ_JK(3) + RXYZ_IJ(3)*RXYZ_JK(1)
     RXYZ_M(3) = RXYZ_IJ(1)*RXYZ_JK(2) - RXYZ_IJ(2)*RXYZ_JK(1)
     RM = SQRT(RXYZ_M(1) ** 2.0 + RXYZ_M(2) ** 2.0 + RXYZ_M(3) ** 2.0)

     RXYZ_M(1) = RXYZ_M(1) / RM
     RXYZ_M(2) = RXYZ_M(2) / RM
     RXYZ_M(3) = RXYZ_M(3) / RM                
     COSM = (RXYZ_IJ(1)*RXYZ_JK(1)+RXYZ_IJ(2)*RXYZ_JK(2)+RXYZ_IJ(3)*RXYZ_JK(3))/RIJ/RJK
     SINM = SQRT(1-COSM*COSM)
     COTANM = COSM/SINM

     !	****** vector N cross product JKxKL
     RXYZ_N(1) = RXYZ_JK(2)*RXYZ_KL(3) - RXYZ_JK(3)*RXYZ_KL(2)
     RXYZ_N(2) =-RXYZ_JK(1)*RXYZ_KL(3) + RXYZ_JK(3)*RXYZ_KL(1)
     RXYZ_N(3) = RXYZ_JK(1)*RXYZ_KL(2) - RXYZ_JK(2)*RXYZ_KL(1)
     RN = SQRT(RXYZ_N(1) ** 2.0 + RXYZ_N(2) ** 2.0 + RXYZ_N(3) ** 2.0)

     RXYZ_N(1) = RXYZ_N(1) / RN
     RXYZ_N(2) = RXYZ_N(2) / RN
     RXYZ_N(3) = RXYZ_N(3) / RN

     COSN = (RXYZ_JK(1)*RXYZ_KL(1)+RXYZ_JK(2)*RXYZ_KL(2)+RXYZ_JK(3)*RXYZ_KL(3))/RJK/RKL
     SINN = SQRT(1-COSN*COSN)
     COTANN = COSN/SINN

     COST = RXYZ_M(1)*RXYZ_N(1) + RXYZ_M(2)*RXYZ_N(2) + RXYZ_M(3)*RXYZ_N(3)
     RXYZ_mN(1) = RXYZ_M(2)*RXYZ_N(3) - RXYZ_M(3)*RXYZ_N(2)
     RXYZ_mN(2) =-RXYZ_M(1)*RXYZ_N(3) + RXYZ_M(3)*RXYZ_N(1)
     RXYZ_mN(3) = RXYZ_M(1)*RXYZ_N(2) - RXYZ_M(2)*RXYZ_N(1)

     SIGNT = sign(1.0, RXYZ_mN(1) * RXYZ_JK(1) + RXYZ_mN(2) * RXYZ_JK(2) + RXYZ_mN(3) * RXYZ_JK(3))

     IF (COST .GT. 1.0 + eps) THEN
        write(*,*) i,j,k,l
        write(*,*)COST 
        stop 'COST .GT. 1.00'
     else if(COST .lT. -(1.00+eps))THEN
        write(*,*) i,j,k,l
        write(*,*)COST 
        stop 'COST .lT. 1.00'
     END IF

     IF (COST .GT. 1.0 ) THEN
        COST=1.0 
     else if(COST .lT. -1.00)THEN
        COST=-1.0 
     END IF

     PHI = ACOS(COST)
     phi = signt * phi
     !    phi = phi - pi*0.5 * anint(phi * ipi*2.0)
     phi = phi*R2D

     PHI_t= 180.0 + phi
     NI = anint(PHI_t / BINO(TOIJKL))

     if(NI .lt. 0)then
        write(*,*) i,j,k,l,signt
        write(*,*) NI,phi_t
        write(*,*)PHI_t / BINO(TOIJKL)
        stop
     end if

     IF(NI .GT. NDATO(TOIJKL)) THEN
        WRITE(*,*)'FATAL ERROR: Entry in out of plane table', TOIJKL,' does not exist'
        WRITE(1,*)'FATAL ERROR: Entry in out of plane table', TOIJKL,' does not exist'
        write(*,*) i,j,k,l
        write(*,*) NI,phi
        WRITE(*,*)'Simulation stopped at Time Step: ',STEP
        STOP
     END IF

     ALPHA=abs(abs(PHI)-abs(ANGLE_OOP(TOIJKL,NI)))/BINO(TOIJKL)
     if(alpha .gt. 1.0)then
        ni = ni+1
        ALPHA=abs(abs(PHI)-abs(ANGLE_OOP(TOIJKL,NI)))/BINO(TOIJKL)
     end if

     FIJKL = OOP_FORCE(TOIJKL,NI)*(1.0-ALPHA)+ ALPHA*OOP_FORCE(TOIJKL,NI+1)
     VIJKL = OOP_POT(TOIJKL,NI)*(1.0-ALPHA)+ ALPHA*OOP_POT(TOIJKL,NI+1)

     V_OOP(OOP_TYPE_LABEL(TOIJKL)) = V_OOP(OOP_TYPE_LABEL(TOIJKL)) + VIJKL

     FXYZ_1(1) = FIJKL * RXYZ_M(1)/SINM/RIJ
     FXYZ_1(2) = FIJKL * RXYZ_M(2)/SINM/RIJ
     FXYZ_1(3) = FIJKL * RXYZ_M(3)/SINM/RIJ

     FXYZ_I(1) = FXYZ_I(1) + FXYZ_1(1)
     FXYZ_I(2) = FXYZ_I(2) + FXYZ_1(2)
     FXYZ_I(3) = FXYZ_I(3) + FXYZ_1(3)

     FXYZ_4(1) = -FIJKL * RXYZ_N(1)/SINN/RKL
     FXYZ_4(2) = -FIJKL * RXYZ_N(2)/SINN/RKL
     FXYZ_4(3) = -FIJKL * RXYZ_N(3)/SINN/RKL

     FXYZ_L(1,L) = FXYZ_L(1,L) + FXYZ_4(1)
     FXYZ_L(2,L) = FXYZ_L(2,L) + FXYZ_4(2)
     FXYZ_L(3,L) = FXYZ_L(3,L) + FXYZ_4(3)

     FXYZ_12(1) = FIJKL*(COTANM*RXYZ_M(1)+COTANN*RXYZ_N(1))/RJK
     FXYZ_12(2) = FIJKL*(COTANM*RXYZ_M(2)+COTANN*RXYZ_N(2))/RJK
     FXYZ_12(3) = FIJKL*(COTANM*RXYZ_M(3)+COTANN*RXYZ_N(3))/RJK

     FXYZ_L(1,J) = FXYZ_L(1,J) - FXYZ_1(1) + FXYZ_12(1)
     FXYZ_L(2,J) = FXYZ_L(2,J) - FXYZ_1(2) + FXYZ_12(2)
     FXYZ_L(3,J) = FXYZ_L(3,J) - FXYZ_1(3) + FXYZ_12(3)

     FXYZ_L(1,K) = FXYZ_L(1,K) - FXYZ_4(1) - FXYZ_12(1)
     FXYZ_L(2,K) = FXYZ_L(2,K) - FXYZ_4(2) - FXYZ_12(2)
     FXYZ_L(3,K) = FXYZ_L(3,K) - FXYZ_4(3) - FXYZ_12(3)

  END DO

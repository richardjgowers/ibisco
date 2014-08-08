  !Inline version of calculating FORCE of an angle    
  DO A = 1, NIJK(I) !298
     J = JANGLEIJK(I,A)
     K = KANGLEIJK(I,A)

     TJ = ITYPE(J)
     TK = ITYPE(K)
     TIJK = IANGT(TI, TJ, TK)

     RXYZ_IJ(1) = RXYZ_I(1) - SXYZ(1,J)
     RXYZ_IJ(2) = RXYZ_I(2) - SXYZ(2,J)
     RXYZ_IJ(3) = RXYZ_I(3) - SXYZ(3,J)

     RXYZ_KJ(1) = SXYZ(1,K) - SXYZ(1,J)
     RXYZ_KJ(2) = SXYZ(2,K) - SXYZ(2,J)
     RXYZ_KJ(3) = SXYZ(3,K) - SXYZ(3,J)

     RIJ = SQRT(RXYZ_IJ(1) ** 2.0 + RXYZ_IJ(2) ** 2.0 + RXYZ_IJ(3) ** 2.0)
     RKJ = SQRT(RXYZ_KJ(1) ** 2.0 + RXYZ_KJ(2) ** 2.0 + RXYZ_KJ(3) ** 2.0)
     CT = (RXYZ_IJ(1)*RXYZ_KJ(1) + RXYZ_IJ(2)*RXYZ_KJ(2) + RXYZ_IJ(3)*RXYZ_KJ(3))/RIJ/RKJ
     THETA = ACOS(CT) * R2D

     !	****** vector M cross product -IJxKJ	
     RXYZ_M(1) = - ( RXYZ_IJ(2)*RXYZ_KJ(3) - RXYZ_IJ(3)*RXYZ_KJ(2))
     RXYZ_M(2) = - (-RXYZ_IJ(1)*RXYZ_KJ(3) + RXYZ_IJ(3)*RXYZ_KJ(1))
     RXYZ_M(3) = - ( RXYZ_IJ(1)*RXYZ_KJ(2) - RXYZ_IJ(2)*RXYZ_KJ(1))
     RM = SQRT(RXYZ_M(1) ** 2.0 + RXYZ_M(2) ** 2.0 + RXYZ_M(3) ** 2.0)

     RXYZ_M(1) = RXYZ_M(1) / RM
     RXYZ_M(2) = RXYZ_M(2) / RM
     RXYZ_M(3) = RXYZ_M(3) / RM

     !	***** vector N cross product MxIJ
     RXYZ_N(1) = ( RXYZ_M(2)*RXYZ_IJ(3) - RXYZ_M(3)*RXYZ_IJ(2) ) / RIJ
     RXYZ_N(2) = (-RXYZ_M(1)*RXYZ_IJ(3) + RXYZ_M(3)*RXYZ_IJ(1) ) / RIJ
     RXYZ_N(3) = ( RXYZ_M(1)*RXYZ_IJ(2) - RXYZ_M(2)*RXYZ_IJ(1) ) / RIJ

     !	***** vector S cross product -MxKJ
     RXYZ_S(1) = -( RXYZ_M(2)*RXYZ_KJ(3) - RXYZ_M(3)*RXYZ_KJ(2) ) / RKJ
     RXYZ_S(2) = -(-RXYZ_M(1)*RXYZ_KJ(3) + RXYZ_M(3)*RXYZ_KJ(1) ) / RKJ
     RXYZ_S(3) = -( RXYZ_M(1)*RXYZ_KJ(2) - RXYZ_M(2)*RXYZ_KJ(1) ) / RKJ

     NI = INT (THETA / BINA(TIJK))
     IF (THETA == 180.0) THEN
        FIJK = BEND_FORCE(TIJK,NI)
        VIJK = BEND_POT(TIJK,NI)
     ELSE
        ALPHA=(THETA-ANGLE(TIJK,NI))/BINA(TIJK)

        FIJK = BEND_FORCE(TIJK,NI)*(1.0D0-ALPHA) &
             + ALPHA*BEND_FORCE(TIJK,NI+1) 

        VIJK = BEND_POT(TIJK,NI)*(1.0D0-ALPHA) &
             + ALPHA*BEND_POT(TIJK,NI+1) 
     END IF

     V_ANGLE(ANGLE_TYPE_LABEL(TIJK)) = V_ANGLE(ANGLE_TYPE_LABEL(TIJK)) + VIJK

     FXYZ_IJK(1)  = FIJK * RXYZ_N(1) / RIJ
     FXYZ_IJK(2)  = FIJK * RXYZ_N(2) / RIJ
     FXYZ_IJK(3)  = FIJK * RXYZ_N(3) / RIJ

     FXYZ_I(1)   = FXYZ_I(1) + FXYZ_IJK(1)
     FXYZ_I(2)   = FXYZ_I(2) + FXYZ_IJK(2)
     FXYZ_I(3)   = FXYZ_I(3) + FXYZ_IJK(3)

     FXYZ_L(1,K) = FXYZ_L(1,K) + FIJK * RXYZ_S(1) / RKJ
     FXYZ_L(2,K) = FXYZ_L(2,K) + FIJK * RXYZ_S(2) / RKJ
     FXYZ_L(3,K) = FXYZ_L(3,K) + FIJK * RXYZ_S(3) / RKJ

     FXYZ_L(1,J) = FXYZ_L(1,J) - FIJK * RXYZ_S(1) / RKJ - FXYZ_IJK(1)
     FXYZ_L(2,J) = FXYZ_L(2,J) - FIJK * RXYZ_S(2) / RKJ - FXYZ_IJK(2)
     FXYZ_L(3,J) = FXYZ_L(3,J) - FIJK * RXYZ_S(3) / RKJ - FXYZ_IJK(3)

  END DO

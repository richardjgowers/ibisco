  !Inline version of calculating FORCE of an angle    
  DO A = 1, NIJK(I) !298
     J = JANGLEIJK(I,A)
     K = KANGLEIJK(I,A)

     TJ = ITYPE(J)
     TK = ITYPE(K)
     TIJK = IANGT(TI, TJ, TK)

     RXIJ = RXI - SX(J)
     RYIJ = RYI - SY(J)
     RZIJ = RZI - SZ(J)

     RXKJ = SX(K) - SX(J)
     RYKJ = SY(K) - SY(J)
     RZKJ = SZ(K) - SZ(J)

     RIJ = SQRT(RXIJ ** 2.0D0 + RYIJ ** 2.0D0 + RZIJ ** 2.0D0)
     RKJ = SQRT(RXKJ ** 2.0D0 + RYKJ ** 2.0D0 + RZKJ ** 2.0D0)
     CT = (RXIJ*RXKJ + RYIJ*RYKJ + RZIJ*RZKJ)/RIJ/RKJ
     THETA = DACOS(CT) * R2D

     !	****** vector M cross product -IJxKJ	
     RXM = - ( RYIJ*RZKJ - RZIJ*RYKJ)
     RYM = - (-RXIJ*RZKJ + RZIJ*RXKJ)
     RZM = - ( RXIJ*RYKJ - RYIJ*RXKJ)
     RM = SQRT(RXM ** 2.0D0 + RYM ** 2.0D0 + RZM ** 2.0D0)

     RXM = RXM / RM
     RYM = RYM / RM
     RZM = RZM / RM

     !	***** vector N cross product MxIJ
     RXN = ( RYM*RZIJ - RZM*RYIJ )/RIJ
     RYN = (-RXM*RZIJ + RZM*RXIJ )/RIJ
     RZN = ( RXM*RYIJ - RYM*RXIJ )/RIJ

     !	***** vector S cross product -MxKJ
     RXS = -( RYM*RZKJ - RZM*RYKJ )/RKJ
     RYS = -(-RXM*RZKJ + RZM*RXKJ )/RKJ
     RZS = -( RXM*RYKJ - RYM*RXKJ )/RKJ

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

     VANGLE = VANGLE + VIJK

     FXIJK  = FIJK * RXN/RIJ
     FYIJK  = FIJK * RYN/RIJ
     FZIJK  = FIJK * RZN/RIJ

     FXI   = FXI + FXIJK
     FYI   = FYI + FYIJK
     FZI   = FZI + FZIJK

     FXL(K) = FXL(K) + FIJK* RXS/RKJ
     FYL(K) = FYL(K) + FIJK* RYS/RKJ
     FZL(K) = FZL(K) + FIJK* RZS/RKJ

     FXL(J) = FXL(J) - FIJK* RXS/RKJ - FXIJK
     FYL(J) = FYL(J) - FIJK* RYS/RKJ - FYIJK
     FZL(J) = FZL(J) - FIJK* RZS/RKJ - FZIJK

  END DO

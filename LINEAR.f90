!*********************************************************************************
      SUBROUTINE LINEAR (N, X, Y, A, B, RHO)

      IMPLICIT REAL*8 (A-H, O-Z)

      DIMENSION X(N), Y(N)

      IF (N .LT. 2) THEN

         WRITE (1, *) 'LINEAR: CANNOT FIT TO ', N, ' DATAPOINTS'

         STOP
      ENDIF

      SX  = 0.0D00
      SXX = 0.0D00
      SXY = 0.0D00
      SY  = 0.0D00
      SYY = 0.0D00

      DO I = 1, N
         SX  = SX  + X(I)
         SXX = SXX + X(I) * X(I)
         SXY = SXY + X(I) * Y(I)
         SY  = SY  + Y(I)
         SYY = SYY + Y(I) * Y(I)
      ENDDO

      DX =  N * SXX - SX * SX
      DY =  N * SYY - SY * SY
      DA =  N * SXY - SX * SY
      DB = SY * SXX - SX * SXY

      A   = DA / DX
      B   = DB / DX


      DXDY = DX * DY
      IF (DXDY .GT. 0.0D00) THEN
         RHO = DA / SQRT(DXDY)
      ELSE
         RHO = 0.0D00
      ENDIF

      RETURN
      END SUBROUTINE

!**********************************************************************************




















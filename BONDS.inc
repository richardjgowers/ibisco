  !Inline method of calculating bonded forces
  
  DO A=1,NBONDS(I)

     J = JBOND(I,A)
     !            IF ANY OF NEIGHBOUR OF ITH SITE HAS INDEX J < I THEN WE 
     !            SKIP THE BONDED INTERATION OF THIS GROUP AS IT HAS BEEN CALCULATED BEFORE 
     IF (J > I) THEN
        TJ = ITYPE(J)
        TIJ = IBONDT(TI, TJ)
        RXYZ_IJ(1) = RXYZ_I(1) - SXYZ(1,J)
        RXYZ_IJ(2) = RXYZ_I(2) - SXYZ(2,J)
        RXYZ_IJ(3) = RXYZ_I(3) - SXYZ(3,J)

        RIJSQ = RXYZ_IJ(1) ** 2.0 + RXYZ_IJ(2) ** 2.0 + RXYZ_IJ(3) ** 2.0

        RIJ = SQRT(RIJSQ)
        NI = INT((RIJ - RBOND(TIJ,0)) / BINB(TIJ))

        !If doesn't fall within table, use last/first entry
        IF(NI .GT. NDATB(TIJ)) THEN
           WRITE(*,*) 'Overextended bond between atoms ', I, ' and ', J
           WRITE(*,*) 'Distance ', RIJ, ' Timestep ', STEP
           WRITE(*,*) 'Continuing on, but this might be a problem!'
           WRITE(1,*) 'Overextended bond between atoms ', I, ' and ', J
           WRITE(1,*) 'Distance ', RIJ, ' Timestep ', STEP
           WRITE(1,*) 'Continuing on, but this might be a problem!'
           NI = NDATB(TIJ)
        ELSE IF(NI .LT. 0) THEN
           NI = 0
        END IF

        !            LINEAR INTEPOLATION
        ALPHA=(RIJ-RBOND(TIJ,NI))/BINB(TIJ)

        FIJ = BOND_FORCE(TIJ,NI)*(1.0D0-ALPHA) &
             + ALPHA*BOND_FORCE(TIJ,NI+1) 

        VIJ = BOND_POT(TIJ,NI)*(1.0D0-ALPHA) &
             + ALPHA*BOND_POT(TIJ,NI+1) 

        V_BOND(BOND_TYPE_LABEL(TIJ)) = V_BOND(BOND_TYPE_LABEL(TIJ)) + VIJ           

        FXYZ_IJ(1)  = FIJ * RXYZ_IJ(1)
        FXYZ_IJ(2)  = FIJ * RXYZ_IJ(2)
        FXYZ_IJ(3)  = FIJ * RXYZ_IJ(3)

        FXYZ_I(1)   = FXYZ_I(1) + FXYZ_IJ(1)
        FXYZ_I(2)   = FXYZ_I(2) + FXYZ_IJ(2)
        FXYZ_I(3)   = FXYZ_I(3) + FXYZ_IJ(3)

        FXYZ_L(1,J) = FXYZ_L(1,J) - FXYZ_IJ(1)
        FXYZ_L(2,J) = FXYZ_L(2,J) - FXYZ_IJ(2)
        FXYZ_L(3,J) = FXYZ_L(3,J) - FXYZ_IJ(3)

     END IF
  END DO

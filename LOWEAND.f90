       subroutine LOWEAND (RIJSQ,RXIJ,RYIJ,RZIJ,I,J)
       USE VAR
       IMPLICIT NONE

       integer J,I
       real R2S,RAN
       real*8 RIJSQ,RIJ,RIJINV,CM1,CM2,MW
       real*8 MFAC,VOE,VE,DV,DVXIJ,DVYIJ,DVZIJ
       real*8 EXIJ,EYIJ,EZIJ
       REAL*8 VXIJ,VYIJ,VZIJ,RXIJ,RYIJ,RZIJ
       DOUBLE PRECISION GAUSSO,SIGGMA,GAUSS

       SIGGMA = 1.D0
       GAUSSO = 0.D0

       CM1 = BEADMASS(I)
       CM2 = BEADMASS(J)
       MW = CM1*CM2/(CM1+CM2)
       MFAC = DSQRT(1.D0/MW)
       RAN = R2S()
       IF (RAN.LE.LAFREQ) THEN
         RIJ = DSQRT(RIJSQ)
         RIJINV = 1.D0/RIJ
         
         EXIJ = RXIJ * RIJINV
         EYIJ = RYIJ * RIJINV
         EZIJ = RZIJ * RIJINV
         
         VXIJ = VX(I) - VX(J) 
         VYIJ = VY(I) - VY(J)
         VZIJ = VZ(I) - VZ(J)
         
         VE = VXIJ*EXIJ + VYIJ*EYIJ + VZIJ*EZIJ
         VOE = GAUSS(SIGGMA,GAUSSO)*TFAC*MFAC
         DV = VOE-VE

         DVXIJ=DV*EXIJ*MW
         DVYIJ=DV*EYIJ*MW
         DVZIJ=DV*EZIJ*MW

         VX(I) = VX(I) + DVXIJ/CM1
         VY(I) = VY(I) + DVYIJ/CM1
         VZ(I) = VZ(I) + DVZIJ/CM1
         
         VX(J) = VX(J) - DVXIJ/CM2
         VY(J) = VY(J) - DVYIJ/CM2
         VZ(J) = VZ(J) - DVZIJ/CM2
         
       ENDIF
       
       RETURN
       END
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       FUNCTION GAUSS(SIGGMA,GAUSSO)
       IMPLICIT NONE
       DOUBLE PRECISION  r,GAUSSO,V1,V2,SIGGMA,GAUSS
       REAL R2S

       r=2.d0
       do while (r.ge.1.d0)
          V1=2.D0*R2S()-1.D0
          V2=2.D0*R2S()-1.D0
          R = V1**2+V2**2
       ENDDO
       GAUSS = V1*SQRT(-2.D0*LOG(r)/(r))
       GAUSS = GAUSSO + SIGGMA*GAUSS
       RETURN
       END
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     SUBROUTINE: 
!     
!     The pseudorandom number generator. 
!     This subroutine has been taken as it is. 
!     
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real function  R2S()

!*************************************************************(GB 12/95)
!
!   Portable long-period (ca. 2.3 * 10^18) random number generator of
!   L'ECUYER [1] with BAYS-DURHAM shuffle [2] and added safeguards as
!   published by PRESS et al. [3] as "ran2" in 1992. In this version
!   (called "R2S" for distinction) no argument is needed, and the
!   initialization can be done by an entry point "R2INIS(iseedS)" with
!   any positive integer "iseedS" as seed. The internal state corres-
!   ponding to  iseedS = 1  is used if no initialization with "R2INIS"
!   is done before the first call of "R2S".
!
!   "R2S" returns a uniform random number in the interval  ]0., 1.[
!   (NOTE: The endpoint values are excluded!) 
!
!
!   *  INITIALIZATION of "R2S":           rdummy = R2INIS(iseedS) 
!
!   *  GENERATION of a random number:     r = R2S()
!
!
!   NOTE:
!
!   *  "rdummy" is a dummy variable of type REAL.
!
!   *  No variable declaractions are necessary in the calling module.
!
!   *  Parameter RNMX=1.-EPS should approximate the largest floating
!      point value that is less than 1. (EPS for a specific machine
!      can be determined with subroutine "MACHAR" from chapt. 20.1 of
!      ref [3]). (For "IBM RISC System/6000" workstations with "AIX XL
!      Fortran/6000", subroutine MACHAR gives  EPS=1.192092896E-07 ).
!
!
!   REFERENCES:
!
!   [1]  P. L'ECUYER, Communications of the ACM, vol. 31 (1988) 742-774.
!   [2]  in D.E. KNUTH, "Seminumerical Algorithms" (2nd. ed.), vol. 2 of
!        "The Art of Computer Programming", Addison-Wesley, Reading, MA
!        (1981) paragraphs 3.2-3.3 .
!   [3]  W.H. PRESS, S.A. TEUKOLSKY, W.T. VETTERLING, and B.P. FLANNERY,
!        "Numerical Recipes in FORTRAN" (2nd ed.), Cambridge University
!        Press, Cambridge (1992), chapt. 7.1 .
!
!
!   TEST OUTPUT (first 35 numbers for iseed=1, in row-wise sequence):
!
!   0.285381  0.253358  0.093469  0.608497  0.903420  0.195873  0.462954
!   0.939021  0.127216  0.415931  0.533739  0.107446  0.399671  0.650371
!   0.027072  0.076975  0.068986  0.851946  0.637346  0.573097  0.902278
!   0.887676  0.372177  0.347516  0.207896  0.569131  0.364677  0.392418
!   0.366707  0.432149  0.153942  0.626891  0.749454  0.341041  0.830393
!
!***********************************************************************

      parameter ( IM1  = 2147483563,     &
                 IM2  = 2147483399,     &
                 AM   = 1./IM1,     &
                 IMM1 = IM1-1,     &
                 IA1  = 40014,     &
                 IA2  = 40692,     &
                IQ1  = 53668,     &
                 IQ2  = 52774,     &
                 IR1  = 12211,     &
                 IR2  = 3791,     &
                 NTAB = 32,     &
                 NDIV = 1+IMM1/NTAB,     &
                 EPS  = 1.2e-7,     &
                 RNMX = 1.-EPS )

      integer     ivS(NTAB)

      save        ivS, iyS, idum1S, idum2S

      data        idum1S/1720212868/, idum2S/1/, iyS/1720212868/

      data        ivS/ 1720212868, 1392842846, 1031324961,  718590712,     &
                        82237802, 1816996195, 1529538438, 1789446856,     &
                       156648835,   52437849, 1441478319,   36906150,     &
                      1269685686, 1644535938,  394503142,  310212663,     &
                      1596049480,    7553450,  322224693,  445508654,     &
                        28884682,  643161691,  407948861,  479214492,     &
                      2124954851,  612891482,  112933431, 1814689225,     &
                        53445315, 1904850491, 1695805043, 1860990862 /

!***********************************************************************
!
!---- Compute "idum1S" by  idum1S = mod( IA1*idum1S, IM1 )  without
!     overflow by SCHRAGE's method (see ref. [3]):

      kS = idum1S / IQ1
      idum1S = IA1*( idum1S - kS*IQ1 ) - kS*IR1
      if  ( idum1S .lt. 0 )  idum1S = idum1S + IM1

!---- Compute "idum2S" likewise:

      kS = idum2S / IQ2
      idum2S = IA2*( idum2S - kS*IQ2 ) - kS*IR2
      if  ( idum2S .lt. 0 )  idum2S = idum2S + IM2

!---- "jS" will be in the range [1 (1) NTAB] :

      jS = 1 + iyS/NDIV

!---- Here "idum1S" is shuffled, and "idum1S" and "idum2S" are combined
!     to produce output:

      iyS = ivS(jS) - idum2S
      ivS(jS) = idum1S
      if  ( iyS .lt. 1 )  iyS = iyS + IMM1

!---- Because users don't expect endpoint values:

      R2S = min( AM*iyS, RNMX )

      return

!>>>> Initialization: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      entry  R2INIS (iseedS)

!---- Be sure to prevent a negative "iseedS" or  iseedS = 0 :

      idum1S = max( iabs(iseedS), 1 )
      idum2S = idum1S

!---- Load the shuffle table (after 8 warm-ups):

      do  10  jS= NTAB+8, 1, -1
         kS = idum1S / IQ1
         idum1S = IA1*( idum1S - kS*IQ1 ) - kS*IR1
         if  ( idum1S .lt. 0 )  idum1S = idum1S + IM1
         if  ( jS .le. NTAB )  ivS(jS) = idum1S
  10  continue
      iyS = ivS(1)
      R2INIS = iyS

      return

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


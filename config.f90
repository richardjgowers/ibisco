subroutine config

use var

implicit none

      integer :: i,k=0,il
      real(kind=rkind) :: x,y,z
      character(len=30) :: flush
 real*8 :: xf,yf,zf

k = 0

      if(timestepcheck .le. 999) then
            write(flush, 99) 'config_',timestepcheck,'.xyz'   
      elseif((timestepcheck .gt. 999) .and. (timestepcheck .le. 9999)) then
            write(flush, 100) 'config_',timestepcheck,'.xyz'  
      elseif((timestepcheck .gt. 9999) .and. (timestepcheck .le. 99999)) then
            write(flush, 101) 'config_',timestepcheck,'.xyz'  
      elseif((timestepcheck .gt. 99999) .and. (timestepcheck .le. 999999)) then
            write(flush, 102) 'config_',timestepcheck,'.xyz'
      elseif((timestepcheck .gt. 999999) .and. (timestepcheck .le. 9999999)) then
            write(flush, 103) 'config_',timestepcheck,'.xyz'
      elseif((timestepcheck .gt. 9999999) .and. (timestepcheck .le. 99999999)) then
            write(flush, 104) 'config_',timestepcheck,'.xyz'
      end if

      99 format(a,I3.3,a)
      100 format(a,I4.4,a)
      101 format(a,I5.5,a)
      102 format(a,I6.6,a)
      103 format(a,I7.7,a)
      104 format(a,I8.8,a)

      OPEN(UNIT=114, FILE = flush, STATUS='replace')
      WRITE(114,*)NATOMS+1
      WRITE(114,*) 't'
      DO Il = 1,natoms
        xf = SX(Il) - boxx*nint(boxxinv*SX(Il)) + boxx*0.5
        yf = Sy(Il) - boxy*nint(boxyinv*Sy(Il)) + boxy*0.5
        zf = Sz(Il) - boxz*nint(boxzinv*Sz(Il)) + boxz*0.5
        IF(TYPE_LABEL(Il) .EQ. 1) THEN
            WRITE(114,9040)'C', xf*10, yf*10, zf*10
        ELSE
            WRITE(114,9040)'O',xf*10, yf*10, zf*10
        END IF
      END DO

      WRITE(114,*) 'S', boxx*5., boxy*5, boxz*5
      CLOSE (114)

!      OPEN(UNIT=114, FILE = flush, STATUS='replace')
!      WRITE(114,*)NATOMS+1
!      WRITE(114,*) 't'
!      DO I = 1, NATOMS
!            IF(IBRDESCR .EQ. 0) THEN
!                  IF(TYPE_LABEL(I) .EQ. 1) THEN
!                        WRITETYPE1 = 'C'
!                        WRITE(114,9040)WRITETYPE1, SX(I)*10, SY(I)*10, SZ(I)*10
!                  ELSE
!                        WRITETYPE1 = 'O'
!                        WRITE(114,9040)WRITETYPE1, SX(I)*10, SY(I)*10, SZ(I)*10
!                  END IF
!            ELSE
!                  WRITE(114,9050)'C', SX(I)*10, SY(I)*10, SZ(I)*10
!            END IF
!      END DO
!      WRITE(114,*) 'S 0.0 0.0 0.0'
!      CLOSE (114)



9040    FORMAT (1X,A8,1X,3 (G21.14,1X))
9050    FORMAT (1X,A8,1X,3 (G21.14,1X))


return
end

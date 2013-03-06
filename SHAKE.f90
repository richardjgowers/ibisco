subroutine shake ()

use var

implicit none

integer :: i,j,it,k,l,tij,a,b,ti,tj,ll
real(kind=rkind) ::  lcfac,dx,dy,dz,rpab
real(kind=rkind) ::  pt11k,pt22k,pt33k
real(kind=rkind) ::  pt12k,pt23k,pt13k
real(kind=rkind) ::  pxab,pyab,pzab,pabsq,rabsq,diffsq
real(kind=rkind) ::  rxab,ryab,rzab,rma,rmb,gab
real(kind=rkind) :: tol2
real(kind=rkind),dimension(natoms) :: pxi,pyi,pzi,vxi,vzi,vyi,ppxi,ppyi,ppzi
logical :: done
logical,dimension(natoms) :: moved,moving

      EK = 0.0D0

      PT11K = 0.0D0
      PT22K = 0.0D0
      PT33K = 0.0D0

      PT12K = 0.0D0
      PT23K = 0.0D0
      PT13K = 0.0D0

        TOL2   = 2.0 * TOL

!   ** Berendsen Thermostat **

      if ((ENSEMBLE == 1).OR.(ENSEMBLE == 2)) THEN
            LCFAC = SQRT(1.0D0+FAC*((TIN/TEMP)-1.0D0))
      else
            LCFAC = 1.0D0
      end if

l=0
a=0
ll=0
do i = 1, nmol

! ** Leap-frog Algorithm **

    do j=1,natm(i)
        l=l+1
        CM = BEADMASS(l)
        VXI(l) = VX(l)
        VYI(l) = VY(l)
        VZI(l) = VZ(l)
        VX(l) = (VXI(l) + DT * FX(l) / CM)
        VY(l) = (VYI(l) + DT * FY(l) / CM)
        VZ(l) = (VZI(l) + DT * FZ(l) / CM)
        SX(l) = SX(l) + DT * VX(l)
        SY(l) = SY(l) + DT * VY(l)
        SZ(l) = SZ(l) + DT * VZ(l)
        !VTX(l) = 0.5D0*(VX(l) + VXI(l))
        !VTY(l) = 0.5D0*(VY(l) + VYI(l))
        !VTZ(l) = 0.5D0*(VZ(l) + VZI(l))

        PXI(l) = SX(l)
        PYI(l) = SY(l)
        PZI(l) = SZ(l)
        PPXI(l) = SX(l)
        PPYI(l) = SY(l)
        PPZI(l) = SZ(l)
    
        MOVING(l) = .FALSE.
        MOVED(l)  = .TRUE.

        FX(l) = 0.0D0
        FY(l) = 0.0D0
        FZ(l) = 0.0D0

    end do

    done = .false.

!    ** SHAKE algorithm **

    do it=1,maxit
        if ( .not. DONE ) then
            do j = 1, natm(i)
                a=a+1
                ti=itype(a)
                do k = 1,nbonds(a)
                    b = jbond(a,k)
                    tj=itype(b)
                    tij = ibondt(TI, TJ)
    write(100,*)it,a,b,tij,constr(tij)
                    if(typeBond(tij) .and. (b > a))then
                        done = .true.
                        if(moved(a) .or. moved(b))then
                            PXAB = PXI(A) - PXI(B)
                            PXAB = PXAB - ANINT ( PXAB * boxxinv ) * boxx
                            PYAB = PYI(A) - PYI(B)
                            PYAB = PYAB - ANINT ( PYAB * boxyinv ) * boxy
                            PZAB = PZI(A) - PZI(B)
                            PZAB = PZAB - ANINT ( PZAB * boxzinv ) * boxz
    
                            PABSQ  = PXAB ** 2 + PYAB ** 2 + PZAB ** 2
                            RABSQ  = constr(tij)**2
                            DIFFSQ = RABSQ - PABSQ ! RABSQ square of the constraint, PABSQ actual value of the distance between A and B
    
                           if ( ABS(DIFFSQ) .GT. ( RABSQ * TOL2 ) ) THEN
    
!                       Unconstrained positions
                               RXAB = PPXI(A) - PPXI(B)
                               RXAB = RXAB - ANINT ( RXAB * boxxinv ) * boxx
                               RYAB = PPYI(A) - PPYI(B)
                               RYAB = RYAB - ANINT ( RXAB * boxxinv ) * boxx
                               RZAB = PPZI(A) - PPZI(B)
                               RZAB = RZAB - ANINT ( RXAB * boxxinv ) * boxx
                               DIFFSQ = (RXAB**2 + RXAB**2 + RXAB**2) - RABSQ
                               RPAB = RXAB * PXAB + RYAB * PYAB + RZAB * PZAB
    
                               IF ( RPAB .LT. ( RABSQ * RPTOL ) ) THEN
                                  write(*,*)RPAB,RABSQ, RPTOL,RABSQ * RPTOL
                                  write(*,*)it,a,b
                                  STOP 'CONSTRAINT FAILURE '
                               ENDIF
        
                               RMA = 1.0 / BEADMASS(a)
                               RMB = 1.0 / BEADMASS(b)
                               GAB = DIFFSQ / ( 2.0 * ( RMA + RMB ) * RPAB )
                               DX  = RXAB * GAB
                               DY  = RYAB * GAB
                               DZ  = RZAB * GAB
        
                               PXI(A) = PPXI(A) + RMA * DX
                               PYI(A) = PPYI(A) + RMA * DY
                               PZI(A) = PPZI(A) + RMA * DZ
                               PXI(B) = PPXI(B) - RMB * DX
                               PYI(B) = PPYI(B) - RMB * DY
                               PZI(B) = PPZI(B) - RMB * DZ
        
                               PPXI(l) = PXI(l)
                               PPYI(l) = PYI(l)
                               PPZI(l) = PZI(l)

                               MOVING(A) = .TRUE.
                               MOVING(B) = .TRUE.
                               DONE = .FALSE.
                            end if !if ( ABS(DIFFSQ) .GT. ( RABSQ * TOL2 ) ) 
                        end if !if(moved(a) .or. moved(b))
                    end if !if(typeBond(tij) .and. (b > a))
                end do
            end do
       end if  !if ( ( .not. DONE ) .and. ( IT .le. MAXIT ) )

        do a = 1, natm(i)
            MOVED(A) = MOVING(A)
            MOVING(A) = .FALSE.
        end do

    end do

    IF ( .NOT. DONE ) THEN
        WRITE(*,'('' TOO MANY CONSTRAINT ITERATIONS '')')
        WRITE(*,'('' MOLECULE '',I5)') I
        STOP
    ENDIF

    do j = 1, natm(i)
        ll=ll+1
        CM = BEADMASS(ll)
        VX(ll) = (VX(ll) + (SX(ll) - PXI(ll)) / DT)*LCFAC
        VY(ll) = (VY(ll) + (SY(ll) - PYI(ll)) / DT)*LCFAC
        VZ(ll) = (VZ(ll) + (SZ(ll) - PZI(ll)) / DT)*LCFAC

        SX(ll) = PXI(ll)
        SY(ll) = PYI(ll)
        SZ(ll) = PZI(ll)

        VTX(ll) = 0.5D0*(VX(ll) + VXI(ll))
        VTY(ll) = 0.5D0*(VY(ll) + VYI(ll))
        VTZ(ll) = 0.5D0*(VZ(ll) + VZI(ll))

        PT11K = PT11K + CM * VTX(ll) ** 2.0
        PT22K = PT22K + CM * VTY(ll) ** 2.0
        PT33K = PT33K + CM * VTZ(ll) ** 2.0

        PT12K = PT12K + CM * VTY(ll) * VTY(ll)
        PT13K = PT13K + CM * VTY(ll) * VTZ(ll)
        PT23K = PT23K + CM * VTY(ll) * VTZ(ll)

    end do


end do


      PT11 = PT11 + PT11K
      PT22 = PT22 + PT22K 
      PT33 = PT33 + PT33K 

      PT12 = PT12 + PT12K
      PT13 = PT13 + PT13K
      PT23 = PT23 + PT23K

      EK = 0.5D0 * (PT11K + PT22K + PT33K)

      TEMP = EK * MKTEMP 

return

end subroutine

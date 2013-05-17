! Debugging
#define DEBUG_OOP .TRUE.! Select to include the debugging stuff of the OOP



!    *******************************************************************

      PROGRAM IBISCO

      USE MODULEPARSING

      USE VAR
      USE RNEMD

      IMPLICIT NONE


!      real etime          ! Declare the type of etime()
!      real elapsed(2),tarray(2),result     ! For receiving user and system time
!      real total          ! For receiving total time
!      real t1,t2
!              character(8)  :: date
!              character(10) :: time
!              character(5)  :: zone
!              integer,dimension(8) :: values

      INTEGER :: I,J,H
      integer :: temp_step

      REAL :: R2INIS,dummy,sumtotBmass

! character(len=2) ::writetype1 


  
      WRITE(*,'(//  '' COARSE GRAINING SIMULATION  '')')
      OPEN ( 115 , FILE = 's-md.out')
      OPEN ( 116 , FILE = 's-md.tp')
      OPEN ( 113 , FILE = 's-md.trj', form='UNFORMATTED', access='SEQUENTIAL')
      WRITE( 115, *)'IBIsCO Revision 22:'
      OPEN (1, FILE='ERROR')
      WRITE(*,*)
      ISTOP = 0

      !DEFINE REDUSED UNITS
      CALL UNIT ()

      !READ DATA.NEW FILE
      CALL RDCONTROL ()
      RLISTSQ = RLIST**2
      RLISTIBRSQ = RLISTIBR**2
      IF (ISTOP == 1) STOP

      IF (PPF_INPUT.EQ.1) THEN
         RDEN1 = 0.D0
         RDEN2 = 0.D0
         RV_AVE1 = 0.D0
         RV_AVE2 = 0.D0

         RVXPRO = 0.d0
         RDENSITY = 0.d0
         RTEEMP = 0.d0
         RPZX = 0.d0

         OPEN ( 173 , FILE = 'Visco_PPF.dat')
         OPEN ( 183 , FILE = 'VGr_Vis.dat' )
         OPEN ( 193 , FILE = 'R_PPF.pro')
         OPEN ( 139 , FILE = 'C_PPF.pro') 
         WRITE(173,*) 'Tsteps, Time(s),C_viscosity(cP),R_viscosity(cP)'
         write(139,*) '1_Binindex, 2_Bin(nm),3_C_velocity(m/s),', &
              '4_C_Density(kg/m^3),5_C_Temp(K),6_C_stressZX(kPa)'
         write(193,*) '1_Binindex, 2_Bin(nm),3_R_velocity(m/s),', &
              '4_R_Density(kg/m^3),5_R_Temp(K),6_R_stressZX(kPa)'
         write(183,*) '1_Binindex, 2_bin(nm), 3_CV_Grad(s^(-1)),',&
              '4_RV_Grad(s^(-1)), 5_CPxz(kPa), 6_RPxz(kPa), 7_C_Visc(cP), 8_R_Visc(cP)'

         !       write(139,*) '1_Binindex, 2_Bin,3_C_velocity,4_C_Density,5_C_Temp,6_C_stressZX'
         !       write(193,*) '1_Binindex, 2_Bin,3_R_velocity,4_R_Density,5_R_Temp,6_R_stressZX'
         !       write(183,*) '1_Binindex, 2_bin(nm), 3_CV_Grad(s^(-1)), 4_RV_Grad(s^(-1)), 5_CPxz(Pa), 6_RPxz(Pa), 7_C_Visc(cP), 8_R_Visc(cP)'
      ENDIF

      VVCONST = 0.5d0
      VVCONST0 = 0.65d0
      !       if (LAINPUT.EQ.1) THEN
      dummy = R2INIS(iseed)
      !       ENDIF

      rfac  = dsqrt(3.d0)     

!     ALLOCATE SOME VARIABLES
      CALL ALLOCATEVAR ()

      !READ PARAM FILE
      CALL RDINTERACT  ()
      IF (ISTOP == 1) STOP

!      if(IBRDESCR .EQ. 0) THEN
!        name_file_virt='virtual'
!      end if





!      IF (IBRDESCR .EQ. 0) THEN   
!        if(virtsite .eq. 1)then         
!            CALL VIRTUAL_SITE()
!        end if
!        IF(ISTOP .EQ. 1) STOP
!      END IF
!     STORE THE POSITION OF THE VIRTUAL SITE

      NVIRTA = 0 !By default NVIRTA is 0

      IF(IBRDESCR .EQ. 0) THEN
         CALL RDVIRTUAL()
      END IF

!      READ CONFIG.NEW FILE (INITIAL VELOCITIES AND CONFIGURATIONS)
      CALL  RDCOOR  ()
      IF (ISTOP == 1) STOP

!     SHIFT THE ATOMS INSIDE THE BOX
      CALL SHIFT ()

!Calculate centers of VS
      CALL VIRTUAL_DEF()

!Calculate mass coeffs
DO I=1,NVIRTA
   sumtotBmass = 0.0D0
   DO J=1,init_numbcomp(I)
      sumtotBmass = sumtotBmass + mass(itype(indx_atm(I,J)))
   END DO
   DO J=1,init_numbcomp(I)
      masscoeff(I,J) = MASS(ITYPE(INDX_ATM(I,J))) / sumtotBmass
   END DO
END DO

!      IF (IBRDESCR .EQ. 0) THEN   
!        if(virtsite .eq. 0)then         
!            CALL VIRTUAL_SITE_COM()
!            IF(ISTOP .EQ. 1) STOP

!        else

!        end if
!      END IF

      TIN = TEMP
      TFAC = DSQRT(TIN)

      gamma = 0.5D0*sigma**2/TIN

      IF ((BOXX <= 3.0D0*RLIST).OR.(BOXY <= 3.0D0*RLIST).OR. &
      (BOXZ <= 3.0D0*RLIST)) THEN
      NEIGHBORLIST = 'NEIGHBOR_NOLIST'
      ELSE
      NEIGHBORLIST = 'NEIGHBOR_WITHLIST'
      END IF

      DO J = 1, NATOMS
            BEADMASS(J) = MASS(ITYPE(J))
      ENDDO

      !MAKE TABLE FORCES
      CALL FTABLE ()

      !     GIVE THE VELOCITIES IF THE INITIAL TIME IS ZERO
      IF (RESTART == 1) THEN
            CALL COMVEL ()
            CALL SCALEV ()
      END IF

      EK = 0.0D0
      VNBOND = 0.0D0
      VBOND  = 0.0D0
      VANGLE = 0.0D0
      VTOR   = 0.0D0
      VOOP   = 0.0D0
      TOTMASS = 0.0D0

!     CALCULATE THE INITIAL TEMPERATURE AND KINETIC ENERGY
      DO 100 I = 1, NATOMS
            EK = EK + BEADMASS(i)*(VX(I)**2.0D0 + VY(I)**2.0D0 + VZ(I)**2.0D0)
            TOTMASS = TOTMASS + MASS(ITYPE(I))
      100 CONTINUE

      EK = 0.5D0 * EK
      TEMP = EK * MKTEMP

!     CREATE LISTS OF BONDS, ANGLES AND TORSIONS
      CALL SETLIS
      IF (ISTOP == 1) STOP
      IF (NOANGLE .EQ. 1) THEN
            IF (NOTORS .EQ. 1) THEN
                  WRITE(*,*) '     ********************* WARNING ********************'
                  WRITE(*,*) '     *  Some angle(s) and torsion(s) are not defined  *'
                  WRITE(*,*) '     *    Simulation will run anyway                  *'
                  WRITE(*,*) '     **************************************************'
                  WRITE(*,*)
            ELSE
                  WRITE(*,*) '     **************** WARNING **************'
                  WRITE(*,*) '     *    Some angles are not defined      *'
                  WRITE(*,*) '     *     Simulation will run anyway      *'
                  WRITE(*,*) '     ***************************************'
                  WRITE(*,*)
            END IF
      ELSE 
            IF (NOTORS .EQ. 1) THEN
                  WRITE(*,*) '     ***************** WARNING ***************'
                  WRITE(*,*) '     *    Some torsions are not defined      *'
                  WRITE(*,*) '     *     Simulation will run anyway        *'
                  WRITE(*,*) '     *****************************************'
                  WRITE(*,*)
           END IF
      END IF

!     CREAT A LIST FOR NON-DONDED INTERACTIONS
!     CALL NON_BOND_ARRAY()
!     IF (ISTOP == 1) STOP


!     IF YOU WANT TO USE GUSSIAN FUNCTION FOR BOND AND BEND INTERACTIONS,
!     READ GUSSIAN FILE AND MAKE TABLE FOR POTENTIAL AND FORCE
      IF (INTERACT == 0) THEN
            CALL RDGAUSSIAN()
            IF (ISTOP == 1) STOP

            CALL BONDTABLE ()
            CALL ANGLETABLE ()
      END IF

!     SETS UP A LIST OF THE TWENTY SIX NEIGHBOURING
!     CELLS OF EACH OF THE SMALL CELLS IN THE CENTRAL BOX.

      NCELLX = INT(BOXX/RLIST)
      NCELLY = INT(BOXY/RLIST)
      NCELLZ = INT(BOXZ/RLIST)

      NUMCELL = NCELLX * NCELLY * NCELLZ 
      MAPSIZ = 13 * NUMCELL 

      MAXNUMCELL = NUMCELL * 10
      MAXMAPSIZ = MAPSIZ * 10

      ALLOCATE (MAP(MAXMAPSIZ))
      ALLOCATE (CELL(NATOMS))
      ALLOCATE (HEAD(MAXNUMCELL))
      ALLOCATE (LCLIST(NATOMS))
      ALLOCATE (NCELL(MAXNUMCELL))
      ALLOCATE (INDX(MAXNUMCELL, NATOMS))
      ALLOCATE (VHEAD(MAXNUMCELL))

      CALL MAPS () 




!       MAKEING THE NEIGHBOUR LIST
      IF (NEIGHBORLIST == 'NEIGHBOR_NOLIST') THEN
!     IF(VIRTSITE .EQ. 0 .OR. VIRTSITE .EQ. 1) THEN
            IF(IBRDESCR .EQ. 0) THEN
!                 CALL VIRT_NEIGHBOR_NOLIST ()
            ELSE
                  CALL NEIGHBOR_NOLIST ()
            END IF
      ELSE
            IF(IBRDESCR .EQ. 0) THEN
                  IF(VIRTSITE .EQ. 0) THEN
                        CALL VIRT_NEIGHBOR_WITHLIST_COM ()
                  ELSE IF (VIRTSITE .EQ. 1) THEN
                        CALL VIRT_NEIGHBOR_WITHLIST ()
                  END IF
            ELSE
                  CALL NEIGHBOR_WITHLIST ()
            END IF
      END IF

      IF (VISCINPUT == 1) Then
            NUMPRO = 0
            NSWP = 0
            TTRANSF = 0.d0
            VXMEAN_SLAB =  0.D0
            SLAB_MDENS = 0.D0
            SLAB_MTEMP = 0.D0  
            OPEN ( 201, FILE = 'md.ntr')
            OPEN ( 202, FILE = 'md.prf')
            write(202,*) '1Slab, 2Rz (nm), 3Vx_prof (m/s), 4Temp (K), 5Density (kg/m^3)' 
            write(201,*) '1Tsteps,2Time(s),3C_Flux(kg/(m.s^2)),', &
                  '4Flux(kg/(m.s^2)),5C_Vgrad(s^(-1)),6Vgrad(s^(-1)),7C_Visc(cP/mPa.s),8Visc(cP)'
            SLAB_THICKNESS = BOXZ/ DBLE(REAL(NUMSLAB))

            DO  J= 1, NATOMS
                  AM(J) = BEADMASS(J)
            END DO
      END IF

!     WRITE THE TOPOLOGY FILE 
      CALL WRITETP ()

      call analysis ()

      !!*************************************************
      !*********** START OF MAIN LOOP*****************
      !*************************************************
       
!        call cpu_time ( t1 )
!              call date_and_time(date,time,zone,values)
!              call date_and_time(DATE=date,ZONE=zone)
!              call date_and_time(TIME=time)
!              call date_and_time(VALUES=values)
!              print '(a,2x,a,2x,a)', date, time, zone
!              print '(8i5))', values


IF (DPDINPUT.EQ.1) THEN    !LOOP FOR DPD: DPD-GW algorithm
    CALL FORCE()
    DO I = 1, NSTEP
        CALL LOOPDPD(I)
    END DO
ELSE ! LOOP FOR NVE/NVT/NPT/LA: leap-frog algorithm
    IF (IBRDESCR .EQ. 0) THEN
       fcutB = RCUTIBR*RCUTIBR
       fcutA = RCUT*RCUT      
        IF(MTS_CHECK .EQ. 0) THEN
            idt=1/DT
            idt2=1/(DT**2)
            idt3=1/(DT**3)
            idt4=1/(DT**4)
            idt5=1/(DT**5)
            mtsdt2 = (DT**2)*0.5
            dt3 = (DT**3)/6.
            dt4 = (DT**4)/24.
            dt5 = (DT**5)/120.
            I0=1
            I1=2
            I2=3
            I3=4
            if(virtsite .eq. 0)then
                if(type_mts .eq. 3)then
                    I4=0
                    I5=0
                    TSMTS = 4
                    Kmts = -TSMTS
                    MTSPARAM2=5
                    CALL LOOP_HYBR_MTS3_COM()
                elseif(type_mts .eq. 4)then
                    I4=5
                    I5=0
                    TSMTS = 5
                    Kmts = -TSMTS
                    MTSPARAM2=6
                    CALL LOOP_HYBR_MTS4_COM()
                elseif(type_mts .eq. 5)then
                    I4=5
                    I5=6
                    TSMTS = 6
                    Kmts = -TSMTS
                    MTSPARAM2=7
                    CALL LOOP_HYBR_MTS5_COM()
                end if
            else
                if(type_mts .eq. 3)then
                    I4=0
                    I5=0
                    TSMTS = 4
                    Kmts = -TSMTS
                    MTSPARAM2=5
                    CALL LOOP_HYBR_MTS3()
                elseif(type_mts .eq. 4)then
                    I4=5
                    I5=0
                    TSMTS = 5
                    Kmts = -TSMTS
                    MTSPARAM2=6
                    CALL LOOP_HYBR_MTS4()
                elseif(type_mts .eq. 5)then
                    I4=5
                    I5=6
                    TSMTS = 6
                    Kmts = -TSMTS
                    MTSPARAM2=7
                    CALL LOOP_HYBR_MTS5()
                end if
            end if
        ELSE IF (MTS_CHECK .EQ. 1) THEN
            if(VIRTSITE .EQ. 0)then
                CALL LOOP_HYBR_COM()
            else
                CALL LOOP_HYBR()
            end if                    
        END IF
    ELSE IF (IBRDESCR .EQ. 1) THEN
        DO I = 1, NSTEP
            CALL LOOP(I)
        END DO
    END IF
END IF

!      total = etime(elapsed)
!      write (*,*) 'End: total=', total, ' user=', elapsed(1),' system=', elapsed(2)

!              call date_and_time(date,time,zone,values)
!              call date_and_time(DATE=date,ZONE=zone)
!              call date_and_time(TIME=time)
!              call date_and_time(VALUES=values)
!              print '(a,2x,a,2x,a)', date, time, zone
!              print '(8i5))', values


!        call cpu_time ( t2 )
!        write ( *, * ) 'Elapsed CPU time = ', t2 - t1

      !*************************************************
      !***********  END OF MAIN LOOP********************
      !*************************************************
      
      CLOSE (201)
      CLOSE (202)
  
      Write(*,*)'**************************************'
      Write(*,*)'********* END OF DYNAMICS ************'
      Write(*,*)'¤.¸.¤.¸.¤°´¯`´¯`°¤.¸.¤.¸.¤°´¯`´¯`°¤.¸'
      Write(*,*)'  (|_(|(|_(|       (|_(|     '
      Write(*,*)'  (="',':','")       (="',':','")      (="',':','")    '
      Write(*,*)'  (,('')('')     (,('')('')    (,('')('')  '
      Write(*,*)'¤.¸.¤.¸.¤°´¯`´¯`°¤.¸.¤.¸.¤°´¯`´¯`°¤.¸'
      Write(*,*)'**************************************'

      END

!*********************************************************************************************

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
!    integer :: temp_step

!    REAL :: R2INIS,dummy

! character(len=2) ::writetype1 

    WRITE(*,'(//  '' IBIsCO TIME  '')')
    OPEN ( 115 , FILE = 's-md.out')
    OPEN ( 116 , FILE = 's-md.tp')
    OPEN ( 113 , FILE = 's-md.trj', form='UNFORMATTED', access='SEQUENTIAL')
    WRITE( 115, *)'IBIsCO Revision 22:'
    OPEN (1, FILE='ERROR')
    WRITE(*,*)
    ISTOP = 0

    !DEFINE REDUCED UNITS
    CALL UNIT ()

    !Read control parameters
    NVIRTA = 0 !By default no virtual sites
    CALL RDCONTROL ()
    IF (ISTOP == 1) STOP 'Failed in RDCONTROL'

!    IF (PPF_INPUT.EQ.1) THEN
!       RDEN1 = 0.D0
!       RDEN2 = 0.D0
!       RV_AVE1 = 0.D0
!       RV_AVE2 = 0.D0
!
!       RVXPRO = 0.d0
!       RDENSITY = 0.d0
!       RTEEMP = 0.d0
!       RPZX = 0.d0
!
!       OPEN ( 173 , FILE = 'Visco_PPF.dat')
!       OPEN ( 183 , FILE = 'VGr_Vis.dat' )
!       OPEN ( 193 , FILE = 'R_PPF.pro')
!       OPEN ( 139 , FILE = 'C_PPF.pro') 
!       WRITE(173,*) 'Tsteps, Time(s),C_viscosity(cP),R_viscosity(cP)'
!       write(139,*) '1_Binindex, 2_Bin(nm),3_C_velocity(m/s),', &
!            '4_C_Density(kg/m^3),5_C_Temp(K),6_C_stressZX(kPa)'
!       write(193,*) '1_Binindex, 2_Bin(nm),3_R_velocity(m/s),', &
!            '4_R_Density(kg/m^3),5_R_Temp(K),6_R_stressZX(kPa)'
!       write(183,*) '1_Binindex, 2_bin(nm), 3_CV_Grad(s^(-1)),',&
!            '4_RV_Grad(s^(-1)), 5_CPxz(kPa), 6_RPxz(kPa), 7_C_Visc(cP), 8_R_Visc(cP)'

       !       write(139,*) '1_Binindex, 2_Bin,3_C_velocity,4_C_Density,5_C_Temp,6_C_stressZX'
       !       write(193,*) '1_Binindex, 2_Bin,3_R_velocity,4_R_Density,5_R_Temp,6_R_stressZX'
       !       write(183,*) '1_Binindex, 2_bin(nm), 3_CV_Grad(s^(-1)), 4_RV_Grad(s^(-1)), 5_CPxz(Pa), 6_RPxz(Pa), 7_C_Visc(cP), 8_R_Visc(cP)'
!    ENDIF

!    VVCONST = 0.5d0
!    VVCONST0 = 0.65d0
    !       if (LAINPUT.EQ.1) THEN
!    dummy = R2INIS(iseed)
    !       ENDIF

!    rfac  = dsqrt(3.0d0)     

    !     ALLOCATE SOME VARIABLES
    CALL ALLOCATEVAR ()

    !READ PARAM FILE
    CALL RDINTERACT  ()
    IF (ISTOP == 1) STOP 'Failed in RDINTERACT'

    !Read coordinate file
    CALL  RDCOOR  ()
    IF (ISTOP == 1) STOP 'Failed in RDCOOR'

    CALL ALLOCATEVAR2()

    NCOARSE = 0
    IF(IBRDESCR .EQ. 0) THEN	
       !Read virtual file			
       CALL RDVIRTUAL()

       IF(ISTOP .eq. 1) STOP 'Failed at RDVIRTUAL'
       !Make linked lists of atoms or bead/VS
       CALL MAKE_LISTS()

       IF(ISTOP .eq. 1) STOP 'Failed at MAKE_LISTS'
    END IF

    ! Calculate centers of VS, either COM or adopt an atom's coords
    CALL VIRTUAL_DEF()

    !     SHIFT THE ATOMS INSIDE THE BOX
    CALL SHIFT ()

    ALLOCATE(POINT_ATOM(NUMATOMS+1), POINT_BEAD(NCOARSE+1))

    TIN = TEMP
    TFAC = DSQRT(TIN)
    gamma = 0.5D0*sigma**2/TIN

    !Decide which type of neighbour list to use
!    IF ((BOXX <= 3.0D0*RLIST).OR.(BOXY <= 3.0D0*RLIST).OR. &
!         (BOXZ <= 3.0D0*RLIST)) THEN
!       NEIGHBORLIST = 'NEIGHBOR_NOLIST'
!    ELSE
!       NEIGHBORLIST = 'NEIGHBOR_WITHLIST'
!    END IF

    !MAKE TABLE FORCES
    CALL FTABLE ()

    !     GIVE THE VELOCITIES IF THE INITIAL TIME IS ZERO
    IF (RESTART == 1) THEN
       CALL COMVEL ()
       CALL SCALEV ()
    END IF

    !     CALCULATE THE INITIAL TEMPERATURE AND KINETIC ENERGY
    DO I = 1, NATOMS
       EK = EK + MASS(ITYPE(I))*(VX(I)**2.0D0 + VY(I)**2.0D0 + VZ(I)**2.0D0)
       TOTMASS = TOTMASS + MASS(ITYPE(I))
    END DO

    EK = 0.5D0 * EK	
    TEMP = EK * MKTEMP	

!     CREATE LISTS OF BONDS, ANGLES AND TORSIONS
    CALL SETLIS
    IF (ISTOP == 1) STOP 'Failed in SETLIS'
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

CONNECTIONS = 0
CONNECTED_TO = 0
MAX_CONTACT = 0 !Is the largest gap between 2 connected things

    CALL BUILD_CONNECTIVITY(NUMATOMS,ATOM,NONBEXC_ATOM)
    IF(ISTOP .eq. 1) STOP 'Failed in BUILD_CONNECTIVITY'

    IF(IBRDESCR .eq. 0) THEN
       CALL BUILD_CONNECTIVITY(NCOARSE,BEAD,NONBEXC_BEAD)
       IF(ISTOP .eq. 1) STOP 'Failed in BUILD_CONNECTIVITY beads'
    END IF

    DO I=1,NBONDS(11)
       WRITE(8000,*) I, JBOND(11,I)
    END DO
    DO I=1,NIJK(11)
       WRITE(8000,*) I, JANGLEIJK(11,I), KANGLEIJK(11,I)
    END DO
    DO I=1,NOANGLEIJK(11)
       WRITE(8000,*) I, NOJANGLEIJK(11,I), NOKANGLEIJK(11,I)
    END DO

DO I=1,NCOARSE
    write(3000,*) BEAD(I), CONNECTIONS(BEAD(I)), (CONNECTED_TO(BEAD(I),J), J=1,CONNECTIONS(BEAD(I))) 
END DO

DO I=1,NUMATOMS
   write(4000,*) ATOM(I), CONNECTIONS(ATOM(I)), (CONNECTED_TO(ATOM(I),J), J=1,CONNECTIONS(ATOM(I)))
END DO
!     IF YOU WANT TO USE GUSSIAN FUNCTION FOR BOND AND BEND INTERACTIONS,
!     READ GUSSIAN FILE AND MAKE TABLE FOR POTENTIAL AND FORCE
!    IF (INTERACT == 0) THEN
!       CALL RDGAUSSIAN()
!       IF (ISTOP == 1) STOP 'Failed in RDGAUSSIAN'

!       CALL BONDTABLE ()
!       CALL ANGLETABLE ()
!    END IF



      CALL MAPS (MAP_ATOM,MAPSIZE_ATOM &
           , NCELLX_ATOM, NCELLY_ATOM, NCELLZ_ATOM) 
      CALL LINKS (HEAD_ATOM,MAXNUMCELL_ATOM,ATOM,NUMATOMS,CELL_ATOM &
           , NCELLX_ATOM, NCELLY_ATOM, NCELLZ_ATOM, LCLIST_ATOM)

      IF(IBRDESCR .eq. 0) THEN
         CALL MAPS (MAP_BEAD,MAPSIZE_BEAD &
              , NCELLX_BEAD, NCELLY_BEAD, NCELLZ_BEAD) 
         CALL LINKS (HEAD_BEAD,MAXNUMCELL_BEAD,BEAD,NUMBEADS,CELL_BEAD &
              , NCELLX_BEAD, NCELLY_BEAD, NCELLZ_BEAD, LCLIST_BEAD)
      END IF

      CALL UPDATE_NEIGHBOURLIST()

!      IF (VISCINPUT == 1) Then
!            NUMPRO = 0
!            NSWP = 0
!            TTRANSF = 0.d0
!            VXMEAN_SLAB =  0.D0
!            SLAB_MDENS = 0.D0
!            SLAB_MTEMP = 0.D0  
!            OPEN ( 201, FILE = 'md.ntr')
!            OPEN ( 202, FILE = 'md.prf')
!            write(202,*) '1Slab, 2Rz (nm), 3Vx_prof (m/s), 4Temp (K), 5Density (kg/m^3)' 
!            write(201,*) '1Tsteps,2Time(s),3C_Flux(kg/(m.s^2)),', &
!                  '4Flux(kg/(m.s^2)),5C_Vgrad(s^(-1)),6Vgrad(s^(-1)),7C_Visc(cP/mPa.s),8Visc(cP)'
!            SLAB_THICKNESS = BOXZ/ DBLE(REAL(NUMSLAB))
!
!            DO  J= 1, NATOMS
!                  AM(J) = BEADMASS(J)
!            END DO
!      END IF

!     WRITE THE TOPOLOGY FILE 
!      CALL WRITETP () FIX THIS

!      call analysis () FIX THIS

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

CALL NEW_LOOP()


!IF (DPDINPUT.EQ.1) THEN    !LOOP FOR DPD: DPD-GW algorithm
!    CALL FORCE()
!    DO I = 1, NSTEP
!        CALL LOOPDPD(I)
!    END DO
!ELSE ! LOOP FOR NVE/NVT/NPT/LA: leap-frog algorithm
!    IF (IBRDESCR .EQ. 0) THEN
!       fcutB = RCUT_BEAD*RCUT_BEAD
!       fcutA = RCUT_ATOM*RCUT_ATOM      
!        IF(MTS_CHECK .EQ. 0) THEN
!            idt=1/DT
!            idt2=1/(DT**2)
!            idt3=1/(DT**3)
!            idt4=1/(DT**4)
!            idt5=1/(DT**5)
!            mtsdt2 = (DT**2)*0.5
!            dt3 = (DT**3)/6.
!            dt4 = (DT**4)/24.
!            dt5 = (DT**5)/120.
!            I0=1
!            I1=2
!            I2=3
!            I3=4
!            if(virtsite .eq. 0)then
!                if(type_mts .eq. 3)then
!                    I4=0
!                    I5=0
!                    TSMTS = 4
!                    Kmts = -TSMTS
!                    MTSPARAM2=5
!                    CALL LOOP_HYBR_MTS3_COM()
!                elseif(type_mts .eq. 4)then
!                    I4=5
!                    I5=0
!                    TSMTS = 5
!                    Kmts = -TSMTS
!                    MTSPARAM2=6
!                    CALL LOOP_HYBR_MTS4_COM()
!                elseif(type_mts .eq. 5)then
!                    I4=5
!                    I5=6
!                    TSMTS = 6
!                    Kmts = -TSMTS
!                    MTSPARAM2=7
!                    CALL LOOP_HYBR_MTS5_COM()
!                end if
!            else
!                if(type_mts .eq. 3)then
!                    I4=0
!                    I5=0
!                    TSMTS = 4
!                    Kmts = -TSMTS
!                    MTSPARAM2=5
!                    CALL LOOP_HYBR_MTS3()
!                elseif(type_mts .eq. 4)then
!                    I4=5
!                    I5=0
!                    TSMTS = 5
!                    Kmts = -TSMTS
!                    MTSPARAM2=6
!                    CALL LOOP_HYBR_MTS4()
!                elseif(type_mts .eq. 5)then
!                    I4=5
!                    I5=6
!                    TSMTS = 6
!                    Kmts = -TSMTS
!                    MTSPARAM2=7
!                    CALL LOOP_HYBR_MTS5()
!                end if
!            end if
!        ELSE IF (MTS_CHECK .EQ. 1) THEN
!            if(VIRTSITE .EQ. 0)then
!                CALL LOOP_HYBR_COM()
!            else
!                CALL LOOP_HYBR()
!            end if                    
!        END IF
!    ELSE IF (IBRDESCR .EQ. 1) THEN
!        DO I = 1, NSTEP
!            CALL LOOP(I)
!        END DO
!    END IF
!END IF

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
      Write(*,*)'  (="',':','" (="',':','")      (="',':','")    '
      Write(*,*)'  (,('')('')     (,('')('')    (,('')('')  '
      Write(*,*)'¤.¸.¤.¸.¤°´¯`´¯`°¤.¸.¤.¸.¤°´¯`´¯`°¤.¸'
      Write(*,*)'**************************************'

    END PROGRAM IBISCO

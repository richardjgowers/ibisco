!	***************************************************************************
SUBROUTINE RDCONTROL ()

  USE MODULEPARSING
  USE VAR
  USE  PAIR
  USE  RNEMD
  IMPLICIT NONE
  INTEGER ALARM
  CHARACTER(80)	TEXT

  !	**************************************************************************
  !Ensemble	        0        (0=NVE, 1=NPT, 2=NVT, 3=EM-V)
  !Temperature	        500      (Temperature/K)
  !Pressure	        100      (Pressure/kPa)
  !Natoms                 2400     (Number of atoms)
  !Nsteps                 1000     (Number of time steps)
  !DT                     1        (Time step in fs)
  !TAUT                   1        (Temperature relaxation time in fs)
  !TAUP                   1        (Pressure relaxation time in fs)
  !Cutoff		        20       (The cutoff distance for non-bonded interaction in nm)
  !Neighbour_list_cutoff  2.3      (Distance in nm up to which pairs of particles are included in the neighbour list.)
  !Update_neighbour_list  15	 (Set intervals of neighbour list updates; 10..20 is usually a good value)
  !Nsampling	        1	 (Interval of sampling of quantities)
  !Ntrajectory            200      (Number of time steps between storing configuration)
  !Halt_Drift	        10       (Interval at which the net drift of the system is reset to zero)
  !Naverage	        100      (Number of time steps between storing average data and restart file)
  !Non-Bonded	        4	 (Use non-bonded potential on 1..4 OR 1..5? 4=1..4, 5=1..5)
  !Interaction	        0	 (Use guassian function for bond and bend interactions or tabels? 0=gaussian,1=table)
  !Lambda                 0.6      (Velocity tunner for DPD thermostat )
  !DPD_cutoff             1.5000    (The cutoff distance for DPD random force and drag force interaction in nm)
  !Weight_type            Linear          (linear: 1-r/rcut or step weight function)
  !Shear_viscosity        y        (N= no calculation of shear visocity, Y= calcuation shear viscosity by RNEMD)
  !Num_RNEMD_slab               20       (number of slabs in RNEMD simulation)
  !Num_RNEMD_exchange            60       ( time step interval of velocity exchange )
  !Num_RNEMD_prof               61       (time step interval of recording RNEMD profile file)
  !Num_RNEMD_trj                61       (time step interval of recording RNEMD trajectory file ) 
  !Poiseuille_flow       y       (Y = Calculate the viscosity by using Periodic Poiseuille Flow (PPF) method, N = not use this method)
  !External_Force        0.02    (External force implimented in PPF method, in pN)
  !Num_slabs_PPF             30      (Num of slabes to collect the velocity profile in PPF method)
  !N_ROLLING_PRO_PPF      10000  (Num of time steps used for rolling average of PPF profiles: Density, Temp, stress, local viscosity)
  !LAFREQ                0.5    (collision frequency of LA thermastat)
  !LA_INPUT              Y      (Y = USING LA, N = NOT USING LA)
  !ISEED                 61909
  !END
  !	**************************************************************************

  OPEN (2, IOSTAT=IOS, FILE='control', STATUS='OLD')

  IF (IOS.NE.0) THEN
     WRITE (1,*)	&
          ' **** FATAL ERROR! File control does not exist ****'
     WRITE (*,*)	&
          ' **** FATAL ERROR! File control does not exist ****'
     ISTOP=1
     RETURN
  END IF

  !	TAKE A INVALID VALUE TO ENSEMBLE AND INTERACT AND RESTART
  ENSEMBLE = 10
  INTERACT = 10
  RESTART  = 10
  VISCINPUT = 10  
  DPDINPUT = 10
  LAINPUT = 10
  ALARM = 10

  DO WHILE (.TRUE.)

     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'ensemble') THEN
        READ (STRNGS(2),*) TEXT
        IF ((TEXT == 'NVE').OR.((TEXT == 'nve'))) ENSEMBLE = 0
        IF ((TEXT == 'NVT').OR.((TEXT == 'nvt'))) ENSEMBLE = 1
        IF ((TEXT == 'NPT').OR.((TEXT == 'npt'))) ENSEMBLE = 2
        IF ((TEXT == 'LA' ).OR.((TEXT == 'la'))) ENSEMBLE = 3
        IF ((TEXT == 'DPD').or.((TEXT == 'dpd'))) ENSEMBLE = 4

        IF ( ENSEMBLE == 10 ) THEN 
           WRITE (1,*)	&
                ' **** FATAL ERROR! You did not input a valid Ensemble in control ****'
           WRITE (1,*)     &
                '*****Possible option: NVT/NVE/NPT/LA/DPD ****************************'
           WRITE (*,*)	&
                ' **** FATAL ERROR! You did not input a valid Ensemble in control ****'
           WRITE (*,*)     &
                '*****Possible option: NVT/NVE/NPT/LA/DPD ****************************'
           ISTOP=1
           RETURN
        END IF
        EXIT
     END IF
     IF (IOS2 /= 0) EXIT
  END DO
  REWIND (2)

  IF (ENSEMBLE==3) LAINPUT = 1
  IF (ENSEMBLE==4) DPDINPUT = 1
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'temperature') THEN
        READ (STRNGS(2),*) TEMP0
        TEMP = TEMP0 / TEMPSCALE
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Temperature// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Temperature// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  REWIND (2)
  ALARM = 10
  DO WHILE (.TRUE.)

     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'pressure') THEN
        READ (STRNGS(2),*) PRESSURE0
        PRESSURE = PRESSURE0/PSCALE
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Pressure// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Pressure// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)

     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'atoms') THEN
        READ (STRNGS(2),*) NATOMS
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Natoms// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Natoms// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)

  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'num_of_time_steps') THEN
        READ (STRNGS(2),*) NSTEP
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Nsteps// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Nsteps// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)

  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'time_step') THEN
        READ (STRNGS(2),*) DT0
        DT = DT0*1.0D-12/TIMESCALE
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //DT// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //DT// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  !*********************************************************************************************

  IBRDESCR = 10
  REWIND (2)

  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'virtual_sites') THEN
        READ (STRNGS(2),*) NVIRTA
        NITEMS = NATOMS + NVIRTA
        IF(NVIRTA .gt. 0) THEN
           IBRDESCR = 0
        ELSE
           IBRDESCR = 1
        END IF
        EXIT
     END IF
     IF (IOS2 /= 0) EXIT
  END DO

  !*********************************************************************************************
  !*********************************************************************************************

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'bead_cutoff') THEN
        READ (STRNGS(2),*) RCUT_BEAD
        RCUTSQ_BEAD = RCUT_BEAD * RCUT_BEAD
        ALARM = 0
        EXIT
     END IF
     IF (IOS2 /= 0) EXIT
  END DO

  IF (ALARM .NE. 0 .AND. IBRDESCR .NE. 1) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //bead_cutoff// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //bead_neighbour_list_cutoff// is missing in //control// file***********'
     ISTOP=1
  END IF
  !****************************************************************************************************************

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'bead_neighbour_list_cutoff') THEN
        READ (STRNGS(2),*) RLIST_BEAD
        ALARM = 0

        IF ( RLIST_BEAD < RCUT_BEAD ) THEN 
           WRITE (1,*)	&
		' **** FATAL ERROR! Neighbour list cutoff must be a bit larger than the cutoff. ****'
           WRITE (*,*)	&
		' **** FATAL ERROR! Neighbour list cutoff must be a bit larger than the cutoff. ****'
           ISTOP=1
           RETURN
        END IF
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM .NE. 0 .AND. IBRDESCR .NE. 1) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //bead_neighbour_list_cutoff// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //bead_neighbour_list_cutoff// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  !****************************************************************************************************************
  IF(IBRDESCR .EQ. 0) THEN
     ALARM = 10
     REWIND (2)
     DO WHILE (.TRUE.)
        READ (2, '(A80)',IOSTAT=IOS2) LINE
        CALL PARSE ()
        IF (STRNGS(1) == 'non_bonded_bead') THEN
           READ (STRNGS(2),*) NONBEXC_BEAD
           ALARM = 0

           IF ((NONBEXC_BEAD .NE. 4).AND.(NONBEXC_BEAD .NE. 5)) THEN
              WRITE (1,*) ' '
              WRITE (1,*) &
                   ' **** FATAL ERROR! FATAL ERROR! ****'
              WRITE (1,*) ' Invalid non-bonded interactions '
              WRITE (1,*) ' **** Check DATA file ****'
              WRITE (*,*) ' '
              WRITE (*,*) &
                   ' **** FATAL ERROR! FATAL ERROR! ****'
              WRITE (*,*) ' Invalid non-bonded interactions '
              WRITE (*,*) ' **** Check DATA file ****'
              ISTOP = 1
              RETURN
           END IF
           EXIT
        END IF

        IF (IOS2 /= 0) EXIT
     END DO
     IF (ALARM.NE.0) THEN
        WRITE (1,*) &
             '*********FATAL ERROR: KEYWORD //Non-Bonded// is missing in //control// file***********'
        WRITE(*,*) &
             '*********FATAL ERROR: KEYWORD //Non-Bonded// is missing in //control// file***********'
        ISTOP=1
     ENDIF
  END IF
  !****************************************************************************************************************

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'update_neighbour_list') THEN
        READ (STRNGS(2),*) NUPDATE
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Update_neighbour_list// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Update_neighbour_list// is missing in //control// file***********'
     ISTOP=1
  ENDIF

  !****************************************************************************************************************

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'com_update') THEN
        READ (STRNGS(2),*) VUPDATE
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF ((ALARM .NE. 0) .AND. (IBRDESCR .NE. 1)) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //com_update// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //com_update// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  IF (nupdate .le. vupdate .and. (IBRDESCR .eq. 0)) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: Neighbour list update must be grater than COM update ***********'
     WRITE(*,*) 
     WRITE(*,*) &
          '*********FATAL ERROR: Neighbour list update must be grater than COM update ***********'
     WRITE(*,*) 
     ISTOP=1
  ENDIF

  !****************************************************************************************************************

  MTS_CHECK = 10
  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF ((STRNGS(1) .EQ. 'MTS') .OR. (STRNGS(1) .EQ. 'mts')) THEN
        READ (STRNGS(2),*) TEXT
        IF ((TEXT(1:1) == 'Y').OR.((TEXT(1:1) == 'y'))) MTS_CHECK = 0
        IF ((TEXT(1:1) == 'N').OR.((TEXT(1:1) == 'n'))) MTS_CHECK = 1
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF ((ALARM .NE. 0) .AND. (IBRDESCR .NE. 1)) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //mts// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //mts// is missing in //control// file***********'
     ISTOP=1
  ENDIF

  !****************************************************************************************************************

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) .EQ. 'tay_ord') THEN
        READ (STRNGS(2),*) type_mts
        ALARM = 0
        if(.not. (type_mts .ge. 3 .and. type_mts .le. 5))then
           WRITE (1,*) &
                '*********FATAL ERROR: tay_ord must be greater than 2 and lower than 6***********'
           WRITE(*,*) &
                '*********FATAL ERROR: tay_ord must be greater than 2 and lower than 6***********'
           ISTOP=1
        end if
        EXIT
     END IF
     IF (IOS2 /= 0) EXIT
  END DO

  IF ((ALARM .NE. 0) .AND. (MTS_CHECK .eq. 0)) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //tay_ord// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //tay_ord// is missing in //control// file***********'
     ISTOP=1
  ENDIF

  !****************************************************************************************************************

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) .EQ. 'approx_ts') THEN
        READ (STRNGS(2),*) mtsparam
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF ((ALARM .NE. 0) .AND. (MTS_CHECK .eq. 0)) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //approx_ts// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //approx_ts// is missing in //control// file***********'
     ISTOP=1
  ENDIF

  !****************************************************************************************************************
  !****************************************************************************************************************
  !****************************************************************************************************************
  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'temperature_coupling_time') THEN
        READ (STRNGS(2),*) TAUT0
        TAUT = TAUT0*1.0D-12/TIMESCALE
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.and.(ENSEMBLE.EQ.1.OR.ENSEMBLE.EQ.2)) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //TAUT// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //TAUT// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'pressure_coupling_time') THEN
        READ (STRNGS(2),*) TAUP0
        TAUP = TAUP0*1.0D-12/TIMESCALE
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.and.ENSEMBLE.EQ.2) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //TAUP// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //TAUP// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)

  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'isothermal_compressibility') THEN
        READ (STRNGS(2),*) BETA0
        BETA = BETA0*PSCALE
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.AND.ENSEMBLE.EQ.2) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //BETA// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //BETA// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'cutoff') THEN
        READ (STRNGS(2),*) RCUT_ATOM
        RCUTSQ_ATOM = RCUT_ATOM * RCUT_ATOM
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //cutoff// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //cutoff// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  !************************************************************************************************************
  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'neighbour_list_cutoff') THEN
        READ (STRNGS(2),*) RLIST_ATOM
        ALARM = 0
        RLISTSQ_ATOM = RLIST_ATOM * RLIST_ATOM
        IF ( RLIST_ATOM < RCUT_ATOM ) THEN 
           WRITE (1,*)	&
		' **** FATAL ERROR! Neighbour list cutoff must be a bit larger than the cutoff. ****'
           WRITE (*,*)	&
		' **** FATAL ERROR! Neighbour list cutoff must be a bit larger than the cutoff. ****'
           ISTOP=1
           RETURN
        END IF
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //neighbour_list_cutoff// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //neighbour_list_cutoff// is missing in //control// file***********'
     ISTOP=1
  ENDIF

  !****************************************************************************************************************

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'sampling') THEN
        READ (STRNGS(2),*) NSAMPLING
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //sampling// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //sampling// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'trajectory') THEN
        READ (STRNGS(2),*) NTRJ
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //trajectory// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //trajectory// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)		
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'halt_drift') THEN
        READ (STRNGS(2),*) HALT_DRIFT
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //halt_Drift// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //halt_Drift// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'rolling_averages') THEN
        READ (STRNGS(2),*) NAVERAGE
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //rolling_averages// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //rolling_averages// is missing in //control// file***********'
     ISTOP=1
  ENDIF

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'non_bonded') THEN
        READ (STRNGS(2),*) NONBEXC_ATOM
        ALARM = 0

        IF ((NONBEXC_ATOM .NE. 4).AND.(NONBEXC_ATOM .NE. 5)) THEN
           WRITE (1,*) ' '
           WRITE (1,*) &
                ' **** FATAL ERROR! FATAL ERROR! ****'
           WRITE (1,*) ' Invalid non-bonded interactions '
           WRITE (1,*) ' **** Check control file ****'
           WRITE (*,*) ' '
           WRITE (*,*) &
                ' **** FATAL ERROR! FATAL ERROR! ****'
           WRITE (*,*) ' Invalid non-bonded interactions '
           WRITE (*,*) ' **** Check control file ****'
           ISTOP = 1
           RETURN
        END IF
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Non-Bonded// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Non-Bonded// is missing in //control// file***********'
     ISTOP=1
  ENDIF

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'interaction') THEN
        READ (STRNGS(2),*) TEXT
        ALARM = 0
        IF ((TEXT(1:1) == 'G').OR.((TEXT(1:1) == 'g'))) INTERACT = 0
        IF ((TEXT(1:1) == 'T').OR.((TEXT(1:1) == 't'))) INTERACT = 1
        IF ( INTERACT == 10 ) THEN 
           WRITE (1,*)	&
                ' **** FATAL ERROR! You did not input a valid Interaction in control ****'
           WRITE (*,*)	&
                ' **** FATAL ERROR! You did not input a valid Interaction in control ****'
           ISTOP=1
           RETURN
        END IF
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Interaction// is missing in //control// file***********'
     write(1,*) &
          '*********The value must be either //GAUSSIAN// or //TABLE//****************************'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Interation// is missing in //control// file***********'
     write(*,*) &
          '*********The value must be either //GAUSSIAN// or //TABLE//****************************'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)


  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'initialize_velocities') THEN
        READ (STRNGS(2),*) TEXT
        ALARM = 0
        IF ((TEXT(1:1) == 'Y').OR.((TEXT(1:1) == 'y'))) RESTART = 1
        IF ((TEXT(1:1) == 'N').OR.((TEXT(1:1) == 'n'))) RESTART = 0

        IF ( RESTART == 10 ) THEN 
           WRITE (1,*)	&
                ' **** FATAL ERROR: KEWORD //Restart_velocities// is missing in //control// file ****'
           WRITE (*,*)	&
                ' **** FATAL ERROR: KEWORD //Restart_velocities// is missing in //control// file ****'
           ISTOP=1
           RETURN
        END IF
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Restart_velocities// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Restart_velocities// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)

  RCUTDPD = RCUT_ATOM
  DO WHILE (.TRUE.)
     READ (2, '(A80)', IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'DPD_LA_cutoff') THEN
        READ (STRNGS(2),*) RCUTDPD
        ALARM = 0
        if (DPDINPUT.EQ.1)  then
           WRITE (*,*) "DPD_cutoff", RCUTDPD
           WRITE (1,*) "DPD_cutoff", RCUTDPD
        endif
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.and.DPDINPUT.EQ.1) THEN
     WRITE (1,*) &
          '*********** NOTE: THE CUTOFF FOR DPD FORCES ARE TOOK AS THE SAME AS NONBONDED INTERACTION BY DEFAULT********* '
     WRITE(*,*) &
          '*********** NOTE: THE CUTOFF FOR DPD FORCES ARE TOOK AS THE SAME AS NONBONDED INTERACTION BY DEFAULT********* '
  ENDIF
  ALARM = 10
  REWIND (2)

  DPD_BONDED = 10
  IF (DPDINPUT.EQ.1.OR.LAINPUT.EQ.1) THEN
     DPD_BONDED = 1       !SWITCH ON THE DPD OR LA THERMOSTAT FOR ALL THE BONDED PAIRS BY DEFAULT
     DO WHILE (.TRUE.)
        READ (2, '(A80)', IOSTAT=IOS2) LINE
        CALL PARSE ()
        IF (STRNGS(1) == 'DPD_LA_on_bonded_pairs') THEN
           READ (STRNGS(2),*) TEXT
           IF ((TEXT(1:1) == 'Y').OR.((TEXT(1:1) == 'y'))) DPD_BONDED = 1
           IF ((TEXT(1:1) == 'N').OR.((TEXT(1:1) == 'n'))) DPD_BONDED = 0
           EXIT
        END IF

        IF (IOS2 /= 0) EXIT
     END DO
     REWIND (2)
  ENDIF

  WDPDT = 2      !DEFAULT TYPE IS LINEAR
  DO WHILE (.TRUE.)  
     READ (2, '(A80)', IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'DPD_weighting_type') THEN
        READ (STRNGS(2),*) TEXT
        ALARM = 0
        IF ((TEXT(1:1) == 'L').OR.((TEXT(1:1) == 'l'))) WDPDT = 2
        IF ((TEXT(1:1) == 'S').OR.((TEXT(1:1) == 's'))) WDPDT = 1
        if (DPDINPUT.EQ.1) WRITE (*,*) "weighting function of DPD: ", TEXT
        EXIT
     END IF
     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.and.DPDINPUT.EQ.1) THEN
     WRITE (1,*) &
          '********* NOTE: KEYWORD //Weight_type// is missing in //control// file***********'
     WRITE (1,*) &
          '********* ITS  VALUE IS SET TO //LINEAR// BY DEFAULT*****************************'
     WRITE(*,*) &
          '********* NOTE: KEYWORD //Weight_type// is missing in //control// file***********'
     WRITE (*,*) &
          '********* ITS VALUE IS SET TO //LINEAR// BY DEFAULT*****************************'
  ENDIF
  ALARM = 10
  REWIND (2)


  DO WHILE (.TRUE.)
     READ (2, '(A80)', IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'RNEMD') THEN
        READ (STRNGS(2),*) TEXT
        !               Write (*,*) "Shear_viscosity", TEXT
        IF ((TEXT(1:1) == 'Y').OR.((TEXT(1:1) == 'y'))) VISCINPUT = 1
        IF ((TEXT(1:1) == 'N').OR.((TEXT(1:1) == 'n'))) VISCINPUT = 0
        IF (VISCINPUT == 1) Write (*,*) "Shear_viscosity", TEXT
        !                     IF ( VISCINPUT  == 10 ) THEN
        !                             WRITE (1,*)     &
        !                        ' **** FATAL ERROR! You did not input a valid VISC file for RNEMD ****'
        !                             WRITE (*,*)     &
        !                        ' **** FATAL ERROR! You did not input a valid VISC file for RNEMD  ****'
        !                         ISTOP=1
        !                         RETURN
        !                     END IF
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  ALARM = 10
  REWIND (2)


  DO WHILE (.TRUE.)
     READ (2, '(A80)', IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'num_RNEMD_slab') THEN
        READ (STRNGS(2),*) NUMSLAB
        ALARM = 0
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.AND.VISCINPUT.EQ.1) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Num_RNEMD_slab// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Num_RNEMD_slab// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)

  DO WHILE (.TRUE.)
     READ (2, '(A80)', IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'num_RNEMD_exchange') THEN
        READ (STRNGS(2),*)  NEXCH
        ALARM = 0
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.AND.VISCINPUT.EQ.1) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Num_RNEMD_exchange// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Num_RNEMD_exchange// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)



  DO WHILE (.TRUE.)
     READ (2, '(A80)', IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'num_RNEMD_prof') THEN
        READ (STRNGS(2),*)  NEMDPROF
        ALARM = 0
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.AND.VISCINPUT.EQ.1) THEN
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Num_RNEMD_prof// is missing in //control// file***********'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Num_RNEMD_prof// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)


  DO WHILE (.TRUE.)
     READ (2, '(A80)', IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'num_RNEMD_trj') THEN
        READ (STRNGS(2),*) NEMDTRAJ
        ALARM = 0
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.AND.VISCINPUT.EQ.1) THEN
     WRITE (1,*) &
          '********* Warning: You are running RNEMD simulation, ************************************'
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Num_RNEMD_trj// is missing in //control// file***********'
     WRITE (*,*) &
          '********* Warning: You are running RNEMD simulation, ************************************'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Num_RNEMD_trj// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)


  !    READ Viscosity input file for RNEMD calculation
  IF ( Viscinput .Eq. 1) then
     write (*,*) '******************************************'
     write (*,*) '*** you are running RNEMD - VISCOSITY ****'
     write (*,*) 'NUMSLAB=',NUMSLAB, 'NEXCH=', NEXCH
     write (*,*) 'NEMDPROF=', NEMDPROF, 'NEMDTRAJ=',NEMDTRAJ
     write (*,*) '******************************************'
  Endif

  !      REAN IN parameters for Periodic Poiseuille Flow method       
  DO WHILE (.TRUE.)
     READ (2, '(A80)', IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'Poiseuille_flow') THEN
        READ (STRNGS(2),*) TEXT
        IF ((TEXT(1:1) == 'Y').OR.((TEXT(1:1) == 'y'))) PPF_INPUT = 1
        IF ((TEXT(1:1) == 'N').OR.((TEXT(1:1) == 'n'))) PPF_INPUT = 0
        if (PPF_INPUT.EQ.1) Write (*,*) "Poiseuille_Flow", TEXT
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  ALARM = 10
  REWIND (2)

  DO WHILE (.TRUE.)
     READ (2, '(A80)', IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'external_force') THEN
        READ (STRNGS(2),*)  EXTER_FX
        ALARM = 0
        if (PPF_INPUT.EQ.1) then
           !   Write (*,*) "External_Force", EXTER_FX, '(pN)'
           !   EXTER_FX = EXTER_FX*1.0e-12/FSCALE
        endif
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.AND.PPF_INPUT.EQ.1) THEN
     WRITE (1,*) &
          '********* Warning: You are running PPF_NEMD simulation, ************************************'
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //EXTER_FX// is missing in //control// file***********'
     WRITE (*,*) &
          '********* Warning: You are running PPF_NEMD simulation, ************************************'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //EXTER_FX// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)

  DO WHILE (.TRUE.)
     READ (2, '(A80)', IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'num_slabs_PPF') THEN
        READ (STRNGS(2),*)  SLIDE
        ALARM = 0
        if (PPF_INPUT.EQ.1) then
           Write (*,*) "num_slabs_PPF", SLIDE
        endif
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.AND.PPF_INPUT.EQ.1) THEN
     WRITE (1,*) &
          '********* Warning: You are running PPF_NEMD simulation, ************************************'
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //Num_slabs_PPF// is missing in //control// file***********'
     WRITE (*,*) &
          '********* Warning: You are running PPF_NEMD simulation, ************************************'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //Num_slabs_PPF// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)

  DO WHILE (.TRUE.)
     READ (2,'(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'rolling_averages_PPF') THEN
        READ (STRNGS(2),*)  AVT
        ALARM = 0
        if (PPF_INPUT.EQ.1) Write (*,*) "N_ROLLING_PRO_PPF", AVT
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.AND.PPF_INPUT.EQ.1) THEN
     WRITE (1,*) &
          '********* Warning: You are running PPF_NEMD simulation, ************************************'
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //N_ROLLING_PRO_PPF// is missing in //control// file***********'
     WRITE (*,*) &
          '********* Warning: You are running PPF_NEMD simulation, ************************************'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //N_ROLLING_PRO_PPF// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)
  IF ( PPF_INPUT .Eq. 1) then
     write (*,*) '*************************************************'
     write (*,*) '* you are using Periodic Poiseuille Flow method *'
     write (*,*) 'Num_slabs_PPF=',SLIDE
     write (*,*) 'External_Force=',EXTER_FX
     write (*,*) 'N_ROLLING_PRO_PPF=',AVT
     write (*,*) '*************************************************'
  Endif

  DO WHILE (.TRUE.)
     READ (2,'(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'LA_collision_frequency') THEN
        READ (STRNGS(2),*)  LAFREQ
        ALARM = 0
        !           Write (*,*) "Frequency_of_LA", LAFREQ
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.AND.LAINPUT.EQ.1) THEN
     WRITE (1,*) &
          '********* Warning: You are running Lowe-Andersen simulation, ************************************'
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //LAFREQ// is missing in //control// file***********'
     WRITE (*,*) &
          '********* Warning: You are running Lowe-Andersen simulation, ************************************'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //LAFREQ// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)



  DO WHILE (.TRUE.)
     READ (2,'(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'iseed') THEN
        READ (STRNGS(2),*)  iseed
        ALARM = 0
        !      Write (*,*) "iseed", iseed
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.AND.(DPDINPUT.EQ.1.OR.LAINPUT.EQ.1)) THEN
     WRITE (1,*) &
          '********* Warning: You are running DPD/Lowe-Andersen simulation, ************************************'
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //ISEED// is missing in //control// file***********'
     WRITE (*,*) &
          '********* Warning: You are running DPD/Lowe-Andersen simulation, ************************************'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //ISEED// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)


  DO WHILE (.TRUE.)
     READ (2,'(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'sigama') THEN
        READ (STRNGS(2),*) sigma
        sigma = sigma*1e-18/(FSCALE*DSQRT(TIMESCALE))
        ALARM = 0
        !     Write (*,*) "sigama (noise strength in DPD) ", sigma ,"pN.ps^(1/2)"
        EXIT
     ENDIF
     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.AND.DPDINPUT.EQ.1) THEN
     WRITE (1,*) &
          '********* Warning: You are running DPD simulation, ************************************'
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //sigama// is missing in //control// file***********'
     WRITE (*,*) &
          '********* Warning: You are running DPD/ simulation, ************************************'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //sigama// is missing in //control// file***********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)

  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'time_steps_before_steady_state') THEN
        READ (STRNGS(2),*) TSET
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.and.(VISCINPUT.EQ.1.OR.PPF_INPUT.EQ.1)) THEN
     WRITE (1,*) &
          '********* FATAL ERROR: YOU ARE RUNNING NEMD (RNEMD/PPF) SIMULATION **************************'
     WRITE (1,*) &
          '*********FATAL ERROR: KEYWORD //T_EQ// is missing in //control// file ***********************'
     WRITE (1,*) &
          '********* NOTE: T_EQ: is a preset number of time-steps needed to reach a stady state ********'
     WRITE (*,*) &
          '********* FATAL ERROR: YOU ARE RUNNING NEMD (RNEMD/PPF) SIMULATION **************************'
     WRITE(*,*) &
          '*********FATAL ERROR: KEYWORD //T_EQ// is missing in //control// file************************'
     WRITE(*,*) &
          '********* NOTE: T_EQ: is a preset number of time-steps needed to reach a stady state ********'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)

  IF (VISCINPUT.EQ.1.AND.PPF_INPUT.EQ.1) THEN
     ISTOP = 1
     WRITE(1,*) &
          '******** FATAL ERROR: RNEMD and PPF are running at the same time, which is absolutly forbidden *******'
     write(1,*) &
          '******** Please switch off one of them ***********************************************************'
     WRITE(*,*) &
          '******** FATAL ERROR: RNEMD and PPF are running at the same time, which is absolutly forbidden ****'
     write(*,*) &
          '******** Please switch off one of them ***********************************************************'
     ISTOP =1
  ENDIF
  ALARM = 10
  REWIND (2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (DPDINPUT.EQ.1.AND.DPD_BONDED.EQ.1)  then
     WRITE (*,*)  '**** DPD THERMOSTAT IS ACTING ON ALL THE POSSIBLE PAIRS ****'
     WRITE (1,*)  '**** DPD THERMOSTAT IS ACTING ON ALL THE POSSIBLE PAIRS ****'
  elseif (DPDINPUT.EQ.1.AND.DPD_BONDED.EQ.0) then
     WRITE (*,*)  '**** DPD THERMOSTAT IS ACTING ONLY ON NONBONDED PAIRS ****'
     WRITE (1,*)  '**** DPD THERMOSTAT IS ACTING ONLY ON NONBONDED PAIRS ****'
  ENDIF
  if (LAINPUT.EQ.1.AND.DPD_BONDED.EQ.1)  then
     WRITE (*,*)  '**** LA THERMOSTAT IS ACTING ON ALL THE POSSIBLE PAIRS ****'
     WRITE (1,*)  '**** LA THERMOSTAT IS ACTING ON ALL THE POSSIBLE PAIRS ****'
  elseif (LAINPUT.EQ.1.AND.DPD_BONDED.EQ.0) then
     WRITE (*,*)  '**** LA THERMOSTAT IS ACTING ONLY ON NONBONDED PAIRS ****'
     WRITE (1,*)  '**** LA THERMOSTAT IS ACTING ONLY ON NONBONDED PAIRS ****'
  ENDIF

  DO  WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'END'.or.STRNGS(1).EQ.'end' &
          .or.STRNGS(1).EQ.'End') THEN
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE (1,*) &
          '********* NOTE: YOU NEED A STRING OF /EDN/ TO END THE /control/ FILE ********'
     WRITE(*,*) &
          '********* NOTE: YOU NEED A STRING OF /EDN/ TO END THE /control/ FILE ********'
     ISTOP=1
  ENDIF
  CLOSE (2)

  RETURN
END SUBROUTINE RDCONTROL

!	*********************************************************************

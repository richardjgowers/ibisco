!> @file
!!
!!
!> @brief Reads the file 'control'
!!
!> @details An example control file is given below \n
!!          The order of the keywords is not important \n
!!
!! @verbatim
!! ensemble                   NVT         (options: NVE, NVT, NPT)
!! temperature                450         (Temperature (K))
!! pressure                   101.3       (Pressure (kPa))
!! atoms                      8184        (Number of particles)
!! num_of_time_steps          1000        (Number of time steps)
!! time_step                  0.001       (Time step in ps)
!! temperature_coupling_time  0.2D0       (Temperature relaxation time in ps)
!! pressure_coupling_time     5.0D0       (Pressure relaxation time in ps)
!! isothermal_compressibility 1.0D-6      (Isothermal compressibility (1/kPa))
!! hybrid_description         Y           (Y or N)
!! virtual_sites              960         (Number of virtual sites in system)
!! MTS                        N           (Options to activate the Multiple Time Step calculation)
!! bead_cutoff                1.200       (Nonbonded cutoff for beads)
!! bead_neighbour_list_cutoff 1.300       (Neighbour list cutoff for beads)
!! non_bonded_bead            4           (Use non-bonded potential for bead on 1..4 OR 1..5)
!! cutoff                     0.900       (regular particle cutoff distance)
!! neighbour_list_cutoff      1.000       (regular particle neighbour list cutoff)
!! update_neighbour_list      20          (Neighbour list update frequency)
!! sampling                   100         (Interval of sampling of quantities)
!! trajectory                 100         (Frequency of saving s-md.trj
!! halt_drift                 100         (How often to reset drift of system)
!! rolling_averages           100         (Frequency of saving s-md.out)
!! non_bonded                 4           (Use non-bonded potential on 1..4 OR 1..5?)
!! initialize_velocities      N           (Generate new random velocities? Y/N)
!! end
!! @endverbatim

SUBROUTINE RDCONTROL ()
  !< \addtogroup read_input

  USE MODULEPARSING
  USE VAR

  IMPLICIT NONE
  INTEGER :: ALARM, IOS, IOS2
  CHARACTER(80)	TEXT

  OPEN (2, IOSTAT=IOS, FILE='control', STATUS='OLD')

  IF (IOS.NE.0) THEN
     WRITE (1,*) ' **** FATAL ERROR! File control does not exist ****'
     WRITE (*,*) ' **** FATAL ERROR! File control does not exist ****'
     ISTOP=1
     RETURN
  END IF

  !	TAKE A INVALID VALUE TO ENSEMBLE AND INTERACT AND RESTART
  ENSEMBLE = 10
  RESTART  = 10 
  ALARM = 10

  DO WHILE (.TRUE.)

     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'ensemble') THEN
        READ (STRNGS(2),*) TEXT
        IF ((TEXT == 'NVE').OR.((TEXT == 'nve'))) ENSEMBLE = 0
        IF ((TEXT == 'NVT').OR.((TEXT == 'nvt'))) ENSEMBLE = 1
        IF ((TEXT == 'NPT').OR.((TEXT == 'npt'))) ENSEMBLE = 2

        IF ( ENSEMBLE == 10 ) THEN 
           WRITE (1,*) ' **** FATAL ERROR! You did not input a valid Ensemble in control ****'
           WRITE (1,*) ' **** Possible options: NVT/NVE/NPT                               ****'
           WRITE (*,*) ' **** FATAL ERROR! You did not input a valid Ensemble in control ****'
           WRITE (*,*) ' **** Possible options: NVT/NVE/NPT                               ****'
           ISTOP=1
           RETURN
        END IF
        EXIT
     END IF
     IF (IOS2 /= 0) EXIT
  END DO
  REWIND (2)

  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'temperature') THEN
        READ (STRNGS(2),*) TEMP0
        TEMP = TEMP0 / TEMPSCALE
        TEMP_IN = TEMP
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

  IF(IBRDESCR .eq. 10) THEN !Error reading this
     WRITE(*,*) '*********FATAL ERROR: KEYWORD //virtual_sites// is missing in //control// file***********'
     WRITE(1,*)'*********FATAL ERROR: KEYWORD //virtual_sites// is missing in //control// file***********'
     ISTOP =1
  END IF

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
     IF (STRNGS(1) == 'initialize_velocities') THEN
        READ (STRNGS(2),*) TEXT
        ALARM = 0
        IF ((TEXT(1:1) == 'Y').OR.((TEXT(1:1) == 'y'))) RESTART = 1
        IF ((TEXT(1:1) == 'N').OR.((TEXT(1:1) == 'n'))) RESTART = 0

        IF ( RESTART == 10 ) THEN 
           WRITE(1,*) '*** FATAL ERROR: KEYWORD //initialize_velocities// is missing in //control// file ***'
           WRITE(*,*) '*** FATAL ERROR: KEYWORD //initialize_velocities// is missing in //control// file ***'
           ISTOP=1
           RETURN
        END IF
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     WRITE(1,*) '*** FATAL ERROR: KEYWORD //Restart_velocities// is missing in //control// file ***'
     WRITE(*,*) '*** FATAL ERROR: KEYWORD //Restart_velocities// is missing in //control// file ***'
     ISTOP=1
  ENDIF
  ALARM = 10
  REWIND (2)

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
     WRITE(1,*) '********* NOTE: YOU NEED A STRING OF /END/ TO END THE /control/ FILE ********'
     WRITE(*,*) '********* NOTE: YOU NEED A STRING OF /END/ TO END THE /control/ FILE ********'
     ISTOP=1
  ENDIF
  CLOSE (2)

  RETURN
END SUBROUTINE RDCONTROL

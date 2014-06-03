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
!! virtual_sites              960         (Number of virtual sites in system)
!! bead_cutoff                1.200       (Nonbonded cutoff for beads)
!! bead_neighbour_list_cutoff 1.300       (Neighbour list cutoff for beads)
!! non_bonded_bead            4           (Use non-bonded potential for bead on 1..4 OR 1..5)
!! cutoff                     0.900       (regular particle cutoff distance)
!! neighbour_list_cutoff      1.000       (regular particle neighbour list cutoff)
!! non_bonded                 4           (Use non-bonded potential on 1..4 OR 1..5?)
!! update_neighbour_list      20          (Neighbour list update frequency)
!! halt_drift                 100         (How often to reset drift of system)
!! trajectory                 100         (Frequency of saving results)
!! initialize_velocities      N           (Generate new random velocities? Y/N)
!! end
!! @endverbatim

SUBROUTINE RDCONTROL ()

  USE MODULEPARSING
  USE VAR
  USE MTS

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
     CALL CONTROL_ERROR('temperature')
     ISTOP=1
  END IF
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
     CALL CONTROL_ERROR('pressure')
     ISTOP=1
  END IF
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
     CALL  CONTROL_ERROR('atoms')
     ISTOP=1
  END IF
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
     CALL CONTROL_ERROR('num_of_time_steps')
     ISTOP=1
  END IF
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
     CALL CONTROL_ERROR('time_step')
     ISTOP=1
  END IF

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

  IF(IBRDESCR .eq. 10) THEN 
     CALL CONTROL_ERROR('virtual_sites')
     ISTOP =1
  END IF

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
     CALL CONTROL_ERROR('bead_cutoff')
     ISTOP=1
  END IF

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'bead_neighbour_list_cutoff') THEN
        READ (STRNGS(2),*) RLIST_BEAD
        ALARM = 0

        IF ( RLIST_BEAD < RCUT_BEAD ) THEN 
           WRITE(1,*) ' **** FATAL ERROR! Neighbour list cutoff must be larger than the cutoff. ****'
           WRITE(*,*) ' **** FATAL ERROR! Neighbour list cutoff must be larger than the cutoff. ****'
           ISTOP=1
           RETURN
        END IF
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM .NE. 0 .AND. IBRDESCR .NE. 1) THEN
     CALL CONTROL_ERROR('bead_neighbour_list_cutoff')
     ISTOP=1
  END IF

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
              WRITE (1,*) ' **** FATAL ERROR! FATAL ERROR! ****'
              WRITE (1,*) ' Invalid non-bonded interactions '
              WRITE (1,*) ' **** Check control file ****'
              WRITE (*,*) ' '
              WRITE (*,*) ' **** FATAL ERROR! FATAL ERROR! ****'
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
        CALL CONTROL_ERROR('non_bonded_bead')
        ISTOP=1
     END IF
  END IF

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
     CALL CONTROL_ERROR('update_neighbour_list')
     ISTOP=1
  END IF

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
     CALL CONTROL_ERROR('temperature_coupling_time')
     ISTOP=1
  END IF
  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'pressure_coupling_time') THEN
        READ (STRNGS(2),*) TAUP0
        TAUP = TAUP0 * 1.0D-12 / TIMESCALE
        ALARM = 0
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0.and.ENSEMBLE.EQ.2) THEN
     CALL CONTROL_ERROR('pressure_coupling_time')
     ISTOP=1
  END IF
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
     CALL CONTROL_ERROR('isothermal_compressibility')
     ISTOP=1
  END IF
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
     CALL CONTROL_ERROR('cutoff')
     ISTOP=1
  END IF

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'neighbour_list_cutoff') THEN
        READ (STRNGS(2),*) RLIST_ATOM
        ALARM = 0
        IF ( RLIST_ATOM < RCUT_ATOM ) THEN 
           WRITE(1,*) ' **** FATAL ERROR! Neighbour list cutoff must be larger than the cutoff. ****'
           WRITE(*,*) ' **** FATAL ERROR! Neighbour list cutoff must be larger than the cutoff. ****'
           ISTOP=1
        END IF
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     CALL CONTROL_ERROR('neighbour_list_cutoff')
     ISTOP=1
  END IF

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
     CALL CONTROL_ERROR('trajectory')
     ISTOP=1
  END IF

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
     CALL CONTROL_ERROR('halt_drift')
     ISTOP=1
  END IF

  ALARM = 10
  REWIND (2)
  DO WHILE (.TRUE.)
     READ (2, '(A80)',IOSTAT=IOS2) LINE
     CALL PARSE ()
     IF (STRNGS(1) == 'non_bonded') THEN
        READ (STRNGS(2),*) NONBEXC_ATOM
        ALARM = 0

        IF ((NONBEXC_ATOM .NE. 4).AND.(NONBEXC_ATOM .NE. 5)) THEN
           WRITE(1,*) ' '
           WRITE(1,*) ' **** FATAL ERROR! FATAL ERROR! ****'
           WRITE(1,*) ' Invalid non-bonded interactions '
           WRITE(1,*) ' **** Check control file ****'
           WRITE(*,*) ' '
           WRITE(*,*) ' **** FATAL ERROR! FATAL ERROR! ****'
           WRITE(*,*) ' Invalid non-bonded interactions '
           WRITE(*,*) ' **** Check control file ****'
           ISTOP = 1
           RETURN
        END IF
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     CALL CONTROL_ERROR('non_bonded')
     ISTOP=1
  END IF

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
           CALL CONTROL_ERROR('initialize_velocities')
           ISTOP=1
           RETURN
        END IF
        EXIT
     END IF

     IF (IOS2 /= 0) EXIT
  END DO
  IF (ALARM.NE.0) THEN
     CALL CONTROL_ERROR('initialize_velocities')
     ISTOP=1
  END IF  

  ALARM = 10
  REWIND(2)
  DO WHILE (.TRUE.)
     READ(2, '(A80)', IOSTAT=IOS2) LINE
     CALL PARSE()
     IF (STRNGS(1) == 'mts') THEN
        READ(STRNGS(2), *) NMTS
        ALARM = 0
     END IF
     IF (IOS2 .ne. 0) EXIT
  END DO
  IF (ALARM .ne. 0) THEN
     CALL CONTROL_ERROR('mts')
     ISTOP = 1
  END IF


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
     CALL CONTROL_ERROR('end')
     ISTOP=1
     RETURN
  END IF

  CLOSE (2)

  RETURN

END SUBROUTINE RDCONTROL

SUBROUTINE CONTROL_ERROR(variable)

  IMPLICIT NONE

  CHARACTER(*) :: variable

  WRITE(1,*) '** Fatal error: Keyword ', variable, ' missing from control file **' 
  WRITE(*,*) '** Fatal error: Keyword ', variable, ' missing from control file **' 
  
  RETURN

END SUBROUTINE CONTROL_ERROR

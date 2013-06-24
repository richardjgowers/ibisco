! Debugging
#define DEBUG_OOP .TRUE.! Select to include the debugging stuff of the OOP
#include "ibi-preprocess.h"
!    *******************************************************************

  PROGRAM IBISCO

    USE MODULEPARSING

    USE VAR

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

open(UNIT=45, FILE='revno')
READ(45,*) revno
close(45)

    WRITE(*,*) 'IBIsCO TIME! Revision ',revno
    OPEN ( 115 , FILE = 's-md.out')
    OPEN ( 116 , FILE = 's-md.tp')
    OPEN ( 113 , FILE = 's-md.trj', form='UNFORMATTED', access='SEQUENTIAL')
    WRITE( 115, *)'IBIsCO Revision', revno
    OPEN (1, FILE='ERROR')
    WRITE(*,*)
    ISTOP = 0

    !DEFINE REDUCED UNITS
    CALL UNIT ()

    !Read control parameters
    NVIRTA = 0 !By default no virtual sites
    CALL RDCONTROL ()
    IF (ISTOP == 1) STOP 'Failed in RDCONTROL'
!rlist_bead=rlist_atom
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
    END IF

    IF(IBRDESCR .EQ. 0) THEN
       CALL MAKE_LISTS()
       IF(ISTOP .eq. 1) STOP 'Failed at MAKE_LISTS'
    ELSE
       !Nonhybrid settings, single list containing all atoms
       NUMATOMS = NATOMS
       ALLOCATE(ATOM(NUMATOMS))
       DO I=1,NATOMS
          ATOM(I) = I
       END DO
    END IF

    IF(IBRDESCR .EQ. 0) THEN
       ! Calculate centers of VS, either COM or adopt an atom's coords
       CALL VIRTUAL_DEF()
    END IF

    !     SHIFT THE ATOMS INSIDE THE BOX
    CALL SHIFT ()

    ALLOCATE(POINT_ATOM(NUMATOMS+1),CELL_ATOM(NUMATOMS),LCLIST_ATOM(NUMATOMS))
    IF(IBRDESCR .eq. 0) THEN
       ALLOCATE(POINT_BEAD(NCOARSE+1),CELL_BEAD(NCOARSE),LCLIST_BEAD(NCOARSE))
    END IF

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
         CALL LINKS (HEAD_BEAD,MAXNUMCELL_BEAD,BEAD,NCOARSE,CELL_BEAD &
              , NCELLX_BEAD, NCELLY_BEAD, NCELLZ_BEAD, LCLIST_BEAD)
      END IF

      CALL UPDATE_NEIGHBOURLIST()

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

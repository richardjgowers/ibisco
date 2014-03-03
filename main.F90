! Debugging
!#define DEBUG_OOP .TRUE.! Select to include the debugging stuff of the OOP
#include "ibi-preprocess.h"
!    *******************************************************************

  PROGRAM IBISCO

    USE MODULEPARSING

    USE VAR

    IMPLICIT NONE

    INTEGER :: I, A

    WRITE(*,*) 'IBIsCO TIME! Revision 77'
    OPEN ( 115 , FILE = 's-md.out')
    OPEN ( 116 , FILE = 's-md.tp')
    OPEN ( 113 , FILE = 's-md.trj', form='UNFORMATTED', access='SEQUENTIAL')
#ifdef TIMING_ON
    OPEN ( 44  , FILE = 'timing.out')
#endif

    WRITE( 115, *)'IBIsCO Revision 77'
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

    MAXNAB_ATOM = 1000 !Max number of neighbours for a single atom
    ALLOCATE(LIST_ATOM(MAXNAB_ATOM,NUMATOMS),CELL_ATOM(NUMATOMS),LCLIST_ATOM(NUMATOMS))
    IF(IBRDESCR .eq. 0) THEN
       MAXNAB_BEAD = 1000
       ALLOCATE(LIST_BEAD(MAXNAB_BEAD,NCOARSE),CELL_BEAD(NCOARSE),LCLIST_BEAD(NCOARSE))
    END IF

    TIN = TEMP
    TFAC = SQRT(TIN) !!double precision
    gamma = 0.5D0*sigma**2/TIN

    !MAKE TABLE FORCES
    CALL FTABLE ()

    !     GIVE THE VELOCITIES IF THE INITIAL TIME IS ZERO
    IF (RESTART == 1) THEN
       CALL COMVEL ()
       CALL SCALEV ()
    END IF

    !     CALCULATE THE INITIAL TEMPERATURE AND KINETIC ENERGY
    EK(1) = 0.0
    EK(2) = 0.0
    WRITE(*,*) NUMATOMS, NUMBEADS
    MKTEMP_ATOM = 2.0D0 / REAL(3.0D0 * (NUMATOMS - 1))
    MKTEMP_BEAD = 2.0D0 / REAL(3.0D0 * (NUMBEADS - 1))
    DO A = 1, NUMATOMS
       I = ATOM(A)
       EK(1) = EK(1) + MASS(ITYPE(I))*(VX(I)**2.0D0 + VY(I)**2.0D0 + VZ(I)**2.0D0)
       TOTMASS = TOTMASS + MASS(ITYPE(I))
    END DO

    DO A=1,NUMBEADS
       I = BEAD(A)
       EK(2) = EK(2) + MASS(ITYPE(I))*(VX(I)**2.0D0 + VY(I)**2.0D0 + VZ(I)**2.0D0)
       TOTMASS = TOTMASS + MASS(ITYPE(I))
    END DO

    WRITE(*,*) MKTEMP, MKTEMP_ATOM, MKTEMP_BEAD
    WRITE(*,*) EK(1), EK(2)

    EK(1) = 0.5D0 * EK(1)
    EK(2) = 0.5D0 * EK(2)
    TEMP = SUM(EK) * MKTEMP	
    TEMP_ATOM = EK(1) * MKTEMP_ATOM
    TEMP_BEAD = EK(2) * MKTEMP_BEAD

    WRITE(*,*) TEMP*TEMPSCALE, TEMP_ATOM*TEMPSCALE, TEMP_BEAD*TEMPSCALE

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
    
    !Read virtangles
    IF(IBRDESCR .eq. 0) THEN
       CALL RDVIRTANGLES()
       IF(ISTOP .eq. 1) STOP 'Failed in RDVIRTANGLES'
    END IF

CONNECTIONS = 0
CONNECTED_TO = 0
MAX_CONTACT = 0 !Is the largest gap between 2 connected things

    CALL BUILD_CONNECTIVITY(NUMATOMS,ATOM,NONBEXC_ATOM)
    IF(ISTOP .eq. 1) STOP 'Failed in BUILD_CONNECTIVITY'

    IF(IBRDESCR .eq. 0) THEN
       CALL BUILD_CONNECTIVITY(NCOARSE,BEAD,NONBEXC_BEAD)
       IF(ISTOP .eq. 1) STOP 'Failed in BUILD_CONNECTIVITY beads'
       CALL RDVIRTBONDS()
       IF(ISTOP .EQ. 1) STOP 'Failed in RDVIRTBONDS'
    END IF

#ifdef REPORT_EXCLUSIONS
    CALL REPORT_EXCLUSIONS() ! Reports all nonbonded exclusions
#endif

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
      CALL WRITEPSF()
!      CALL WRITETP () FIX THIS

!      call analysis () FIX THIS

      !!*************************************************
      !*********** START OF MAIN LOOP*****************
      !*************************************************

WRITE(*,*)
WRITE(*,*) 'Beginning simulation'
WRITE(*,*)

CALL NEW_LOOP()


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

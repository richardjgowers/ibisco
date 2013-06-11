SUBROUTINE TIMING()

  USE VAR
  USE OMP_LIB

  IMPLICIT NONE

  REAL*8 :: t_MOVING !Because of how program is arranged, this has to be calculated separately

  t_MOVING = t_SHIFT(2) - t_SHIFT(1) + &
       t_MOVE(2) - t_MOVE(1) + &
       t_MOMENTUM(2) - t_MOMENTUM(1) + &
       t_SCALEBP(2) - t_SCALEBP(1)

  WRITE(44,*) 'Step            ', STEP
  WRITE(44,*) 'Total_loop_time ', t_LOOP(2) - t_LOOP(1)
  WRITE(44,*) ''
  WRITE(44,*) 'Force_calculation  ', t_FORCECALC(2) - t_FORCECALC(1)
  WRITE(44,*) '  Virtualdef         ', t_VIRTUAL_DEF(2) - t_VIRTUAL_DEF(1)
  WRITE(44,*) '  Neighbourlist      ', t_UPDATE_NEIGHBOURLIST(2) - t_UPDATE_NEIGHBOURLIST(1)
  WRITE(44,*) '    Maps_atom          ',t_MAPS_atom(2) - t_MAPS_atom(1)
  WRITE(44,*) '    Links_atom         ',t_LINKS_atom(2) - t_LINKS_atom(1)
  WRITE(44,*) '    List_atom          ', t_NEW_NEIGHBOUR_atom(2) - t_NEW_NEIGHBOUR_atom(1)
  WRITE(44,*) '    Maps_bead          ',t_MAPS_bead(2) - t_MAPS_bead(1)
  WRITE(44,*) '    Links_bead         ',t_LINKS_bead(2) - t_LINKS_bead(1)
  WRITE(44,*) '    List_bead          ', t_NEW_NEIGHBOUR_bead(2) - t_NEW_NEIGHBOUR_bead(1)
  WRITE(44,*) '  Force              ', t_NEW_FORCE(2) - t_NEW_FORCE(1)
  WRITE(44,*) '    Nonbonded_atom     ',t_NONBONDED_atom(2) - t_NONBONDED_atom(1)
  WRITE(44,*) '    Nonbonded_bead     ',t_NONBONDED_bead(2) - t_NONBONDED_bead(1)
  WRITE(44,*) '    Distribute_vs      ',t_DISTRIBUTE_VSFORCE(2) - t_DISTRIBUTE_VSFORCE(1)
  WRITE(44,*) '    Bonded             ',t_BONDED_FORCE(2) - t_BONDED_FORCE(1)
  WRITE(44,*) ''
  WRITE(44,*) 'Moving             ', t_MOVING
  WRITE(44,*) '  Shift              ',t_SHIFT(2) - t_SHIFT(1)
  WRITE(44,*) '  Move               ',t_MOVE(2) - t_MOVE(1)
  WRITE(44,*) '  Momentum           ',t_MOMENTUM(2) - t_MOMENTUM(1)
  WRITE(44,*) '  ScaleBP            ',t_SCALEBP(2) - t_SCALEBP(1)
  WRITE(44,*) ''
  WRITE(44,*) 'Reporting_results  ',t_REPORTRESULTS(2) - t_REPORTRESULTS(1)
  WRITE(44,*) '  Average            ',t_AVERAGE(2) - t_AVERAGE(1)
  WRITE(44,*) '  Trj_writing        ',t_WRITETRJ(2) - t_WRITETRJ(1)
  WRITE(44,*) '  Output             ',t_OUTPUT(2) - t_OUTPUT(1)
  WRITE(44,*) '________________________'

  FLUSH(44)

  RETURN

END SUBROUTINE TIMING

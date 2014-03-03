try:
    import MDAnalysis as mda
except ImportError:
    raise ImportError("You need MDAnalysis version 0.8 to run tests")

import numpy as np
from numpy.testing import *
from nose.tools import with_setup

class reference_hybrid_data(object):
    """Test case made from hybrid polyamide
    Tests everything except OOP torsions, MTS and COM virtual sites
    """
    #Short is first frames, 10 steps in
    #Long is 4th frames, 40 steps in
    #Coords and velocities of first 3 atoms
    #Box dimensions

    ref_coords_short = np.array([[  67.87813568, -117.15910339,   20.33826447],
                                 [  70.13469696, -115.7314682 ,   19.97340965],
                                 [  70.47128296, -115.85394287,   21.00744247]], dtype=np.float32)

    ref_velos_short = np.array([[ -7.05872011,   2.2907176 ,   3.69597816],
                                [ 10.29509735,   2.6859746 ,   1.91628158],
                                [-18.3086338 ,  -5.09177303,  12.15194225]], dtype=np.float32)

    ref_box_short = np.array([ 55.50342178,  55.50342178,  55.50342178,  90.        , #u.dimensions
                               90.        ,  90.        ], dtype=np.float32)

    ref_coords_long = np.array([[  67.76287842, -117.16596985,   20.47505951],
                                [  70.21298218, -115.93202209,   20.06541824],
                                [  70.69917297, -116.15378571,   21.01545334]], dtype=np.float32)

    ref_velos_long = np.array([[ -4.0683465 ,   5.58288145,   3.43824601],
                               [ -5.11029816,  -8.98984718,   2.34498286],
                               [ 14.03837681, -16.75909233,  -6.12066698]], dtype=np.float32)

    ref_box_long = np.array([ 55.56259155,  55.56259155,  55.56259155,  90.        ,
                              90.        ,  90.        ], dtype=np.float32)



class reference_atomistic_data(object):
    """Test case made from atomistic polyamide
    Tests not using hybrid options
    """
    ref_coords_short = np.array([[  67.27505493, -117.7746582 ,   19.0885849 ], #First three coordinates
                                 [  68.14879608, -118.17668152,   18.57328033], #step 10
                                 [  66.50193787, -117.52475739,   18.3514843 ]],dtype=np.float32)

    ref_velos_short = np.array([[ -3.97902226,   3.21878624,   1.12899768], #First three velocities
                               [-20.3668499 ,  -8.04893303, -19.79080391], #step 10
                               [ 15.88148975,  23.6400032 , -17.29924011]],dtype=np.float32)

    ref_box_short = np.array([ 55.42282867,  55.42282867,  55.42282867,  90.        ,
                               90.        ,  90.        ], dtype=np.float32)


    ref_coords_long = np.array([[  67.18212891, -117.65636444,   19.08372307],
                                [  68.03952026, -117.93842316,   18.47685623],
                                [  66.59827423, -116.88020325,   18.58401489]], dtype=np.float32)

    ref_velos_long = np.array([[ -5.76776314,  12.34260082,   1.68799782],
                               [-10.26216602,  -5.9834342 ,   1.50855136],
                               [  6.8113637 ,   6.97455645, -30.62001038]], dtype=np.float32)

    ref_box_long = np.array([ 55.42282867,  55.42282867,  55.42282867,  90.        ,
                              90.        ,  90.        ], dtype=np.float32)


class AtomisticTest(TestCase, reference_atomistic_data):
    def setUp(self):
        self.universe = mda.Universe('atomisticPA.psf','atomistic.trj',format='TRZ')
        self.trz = self.universe.trajectory
        self.ts = self.universe.trajectory.ts
        self.prec = 5

    def tearDown(self):
        del self.universe
        del self.trz
        del self.ts

    def test_short_times(self):
        self.trz.rewind()
        assert_almost_equal(self.ref_coords_short, self.universe.atoms.coordinates()[0:3], self.prec)
        assert_almost_equal(self.ref_velos_short,  self.universe.atoms.velocities[0:3], self.prec)
        assert_almost_equal(self.ref_box_short,    self.universe.dimensions, self.prec)
        
    def test_long_times(self):
        self.trz.rewind()
        self.trz.next()
        self.trz.next()
        self.trz.next()
        assert_almost_equal(self.ref_coords_long, self.universe.atoms.coordinates()[0:3], self.prec)
        assert_almost_equal(self.ref_velos_long,  self.universe.atoms.velocities[0:3], self.prec)
        assert_almost_equal(self.ref_box_long,    self.universe.dimensions, self.prec)

class HybridTest(TestCase, reference_hybrid_data):
    def setUp(self):
        self.universe = mda.Universe('hybridPA.psf','hybrid.trj',format='TRZ')
        self.trz = self.universe.trajectory
        self.ts = self.universe.trajectory.ts
        self.prec = 5 # desired precision for tests to be ran at

    def tearDown(self):
        del self.universe
        del self.trz
        del self.ts

    def test_short_times(self): #test after 1 step
        self.trz.rewind()
        assert_almost_equal(self.ref_coords_short, self.universe.atoms.coordinates()[0:3], self.prec)
        assert_almost_equal(self.ref_velos_short,  self.universe.atoms.velocities[0:3], self.prec)
        assert_almost_equal(self.ref_box_short,    self.universe.dimensions, self.prec)

    def test_long_times(self): #test after many steps to catch different errors
        self.trz.rewind()
        self.trz.next()
        self.trz.next()
        self.trz.next()
        assert_almost_equal(self.ref_coords_long, self.universe.atoms.coordinates()[0:3], self.prec)
        assert_almost_equal(self.ref_velos_long,  self.universe.atoms.velocities[0:3], self.prec)
        assert_almost_equal(self.ref_box_long,    self.universe.dimensions, self.prec)

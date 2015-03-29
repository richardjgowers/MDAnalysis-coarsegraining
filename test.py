# coarse graining in MDAnalysis
# Copyright (c) 2015 Richard J Gowers
# Released under the GNU Lesser General Public License, version 2 or later.

from numpy.testing import TestCase, assert_equal

from MDAnalysis import Universe
from cguniverse import CGUniverse


class TestCGUniverse(TestCase):
    def setUp(self):
        mapping = ([0, 1, 2],
                   [3, 4, 5],
                   [6, 7, 8])

        self.at_u = Universe('new.gro')
        self.u = CGUniverse('new.gro', mapping=mapping)

    def tearDown(self):
        del self.at_u
        del self.u

    def test_pos1(self):
        res = self.at_u.residues[0]

        assert_equal(res.centerOfMass(), self.u.beads[0].position)

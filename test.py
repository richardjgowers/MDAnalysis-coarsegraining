# coarse graining in MDAnalysis
# Copyright (c) 2015 Richard J Gowers
# Released under the GNU Lesser General Public License, version 2 or later.

from numpy.testing import TestCase, assert_allclose

from MDAnalysis import Universe
from cguniverse import CGUniverse


class TestCGUniverse(TestCase):
    def setUp(self):
        mapping = {'SOL': [[0, 1, 2]]}

        self.at_u = Universe('new.gro')
        self.u = CGUniverse('new.gro', mapping=mapping)

    def tearDown(self):
        del self.at_u
        del self.u

    def test_pos1(self):
        """First "atom" position should be COM of first residue"""
        res = self.at_u.residues[0]

        assert_allclose(res.centerOfMass(), self.u.beads[0].position)

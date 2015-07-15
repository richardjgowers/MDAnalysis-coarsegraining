# coarse graining in MDAnalysis
# Copyright (c) 2015 Richard J Gowers
# Released under the GNU Lesser General Public License, version 2 or later.

from numpy.testing import TestCase, assert_equal, assert_allclose

from MDAnalysis import Universe
from cguniverse import CGUniverse


class TestCGUniverse(object):
    def setUp(self):
        mapping = {'SOL': [[0, 1, 2]]}

        self.at_u = Universe('testdata/new.gro')
        self.u = CGUniverse('testdata/new.gro', mapping=mapping)

    def tearDown(self):
        del self.at_u
        del self.u

    def test_creation(self):
        assert_equal(len(self.u.atoms), 3)

    def test_pos1(self):
        """First "atom" position should be COM of first residue"""
        res = self.at_u.residues[0]

        assert_allclose(res.centerOfMass(), self.u.atoms[0].position)

    def test_bead_len(self):
        assert_equal(len(self.u.atoms[0].atoms), 3)

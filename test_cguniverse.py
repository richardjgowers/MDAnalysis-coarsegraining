# coarse graining in MDAnalysis
# Copyright (c) 2015 Richard J Gowers
# Released under the GNU Lesser General Public License, version 2 or later.

from numpy.testing import TestCase, assert_equal, assert_allclose

import MDAnalysis as mda
from cguniverse import CGUniverse

import pytest


@pytest.fixture
def at_universe():
    return mda.Universe('testdata/new.gro')

@pytest.fixture
def cg_universe():
    mapping = {'SOL': [[0, 1, 2]]}

    return CGUniverse('testdata/new.gro', mapping=mapping)


def test_creation(cg_universe):
    assert_equal(len(cg_universe.atoms), 3)


def test_pos1(at_universe, cg_universe):
    """First "atom" position should be COM of first residue"""
    for res, atom in zip(at_universe.residues, cg_universe.atoms):
        assert_allclose(res.atoms.center_of_mass(), atom.position)

# coarse graining in MDAnalysis
# Copyright (c) 2015, 2019 Richard J Gowers
# Released under the GNU Lesser General Public License, version 2 or later.

from numpy.testing import assert_allclose

import MDAnalysis as mda
from cguniverse import CGUniverse

import pytest


@pytest.fixture
def water_at_universe():
    return mda.Universe('testdata/water.gro')

@pytest.fixture
def water_cg_universe():
    mapping = {'SOL': [[0, 1, 2]]}

    return CGUniverse('testdata/water.gro', mapping=mapping)


def test_creation(water_cg_universe):
    assert len(water_cg_universe.atoms) == 3


def test_pos1(water_at_universe, water_cg_universe):
    """First "atom" position should be COM of first residue"""
    for res, atom in zip(water_at_universe.residues, water_cg_universe.atoms):
        assert_allclose(res.atoms.center_of_mass(), atom.position)


@pytest.fixture
def oct_atu():
    return mda.Universe('testdata/octanol.gro')

@pytest.fixture
def oct_cgu():
    return CGUniverse('testdata/octanol.gro',
                      mapping={'OcOH':[list(range(11)), list(range(11, 27))]})


def test_residue_com(oct_atu, oct_cgu):
    # check that the COM of the first residue is unchanged
    assert_allclose(oct_atu.residues[0].atoms.center_of_mass(),
                    oct_cgu.residues[0].atoms.center_of_mass())

# coarse graining in MDAnalysis
# Copyright (c) 2015 Richard J Gowers
# Released under the GNU Lesser General Public License, version 2 or later.

"""Universe object for coarse graining in MDAnalysis toolkit.

Works like a regular MDAnalysis Universe, only requires a mapping as well

eg:

cgu = CGUniverse('out.gro','out.trr', mapping)

cgu.atoms # returns an AtomGroup of the beads

cgu.atoms.positions # returns the positions of the beads

"""
import MDAnalysis as mda
from MDAnalysis.core.AtomGroup import Atom, AtomGroup

from cgtraj import CGTraj


class Bead(Atom):
    """Building block for CG Universe

    Should act like an Atom, but reports a mass weighted average on the
    atoms within.
    """

    __slots__ = (
        "atoms",
        "number", "id", "name", "type", "resname", "resid", "segid",
        "mass", "charge", "residue", "segment",
        "__universe",
        "radius", "bfactor", "resnum", "serial", "altLoc")

    def __init__(self, atoms, *args, **kwargs):
        """atoms should be an AtomGroup representing the bead"""
        self.atoms = AtomGroup(atoms)
        super(Bead, self).__init__(*args, ** kwargs)

    def __repr__(self):
        return ("<Bead {idx}: {name} of type {t} of resname {rname}, "
                "resid {rid} and segid {sid}{altloc}"
                " representing {natoms} atoms>".format(
                    idx=self.number + 1, name=self.name, t=self.type,
                    rname=self.resname, rid=self.resid, sid=self.segid,
                    altloc="" if not self.altLoc
                    else " and altloc {}".format(self.altLoc),
                    natoms=len(self._atoms)))


class CGUniverse(mda.Universe):
    """A universe of CG particles"""
    def __init__(self, *args, **kwargs):
        """Initialise like a normal MDAnalysis Universe but give the mapping keyword.

        Mapping must be a dictionary with resnames as keys.
        Each resname must then correspond to a list of list of indices,
        signifying how to split up a single residue into many beads.

        eg:
        mapping = {'A':[[0, 1, 2], [3, 4, 5]],
                   'B':[[0, 1], [2, 3], [4, 5]]}

        Would split 'A' residues into 2 beads of 3 atoms,
        and 'B' residues into 3 beads of 2 atoms.

        Note that the indices inside the mapping refer to the atom's position
        inside the residue, not their absolute position in the Universe.

        """
        try:
            mapping = kwargs.pop('mapping')
        except KeyError:
            raise ValueError("CGUniverse requires the mapping keyword")
        # Atomistic Universe
        # TODO: wrap this in try/except and process errors in constructing
        # atomistic universe properly.
        self.atu = mda.Universe(*args, **kwargs)

        # Coarse grained Universe
        # Make a blank Universe for myself.
        super(CGUniverse, self).__init__()

        # Fake up some beads
        # TODO: Names, better topology etc
        beads = self._apply_map(mapping)
#        beads = [AtomGroup(self.atu.atoms[m]) for m in mapping]

        beadgroup = [Bead(b, i, 'BEAD', 'BEAD',
                          b[0].resname, b[0].resid, b[0].segid,
                          1.0, 1.0) for i, b in enumerate(beads)]
        for b in beadgroup:
            b.universe = self
        self.universe = self
        self.atoms = AtomGroup(beadgroup)
        self.beads = beadgroup

        # TODO: Wrap this in try except to catch errors in making fake Reader
        # This replaces load_new in a traditional Universe
        self.trajectory = CGTraj(self.atu.trajectory, self.beads)

    def _apply_map(self, mapping):
        """Apply the mapping scheme to the beads"""
        beads = []
        for r, m in mapping.items():
            for res in self.atu.selectAtoms('resname {}'.format(r)).residues:
                beads.extend([AtomGroup(res.atoms[idx]) for idx in m])

        return beads

    def __repr__(self):
        return "<CG Universe with {} beads>".format(len(self.beads))

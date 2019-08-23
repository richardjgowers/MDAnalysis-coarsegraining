# coarse graining in MDAnalysis
# Copyright (c) 2015, 2019 Richard J Gowers
# Released under the GNU Lesser General Public License, version 2 or later.

"""Universe object for coarse graining in MDAnalysis toolkit.

Works like a regular MDAnalysis Universe, only requires a mapping as well

eg:

cgu = CGUniverse('out.gro','out.trr', mapping)

cgu.atoms # returns an AtomGroup of the beads

cgu.atoms.positions # returns the positions of the beads

"""
import MDAnalysis as mda
from MDAnalysis.coordinates import base
from MDAnalysis.core import topologyattrs as ta
from MDAnalysis.core import topology
import numpy as np


class CGTraj(base.ProtoReader):
    """Fakes a coarse grained trajectory object

    Takes an atomistic trajectory and list of beads and manipulates both
    to recreate a reader for a coarse grained Universe.

    Would also probably work as a standalone thingy for writing out
    coarse grained trajectories.
    """
    def __init__(self, atomistic_traj, mapping, weights):
        """
        Parameters
        ----------
        atomistic_traj : mda.Reader
          atomistic trajectory
        mapping : numpy array
          which bead
        weights : numpy array
          for each atom, what weight
        """
        super(base.ProtoReader, self).__init__()

        self._t = atomistic_traj
        self._mapping = mapping
        self._weights = weights

        self.n_atoms = mapping.max() + 1

        self.ts = base.Timestep(self.n_atoms)
        self._fill_ts(self._t.ts)

    def _read_next_timestep(self):
        """Return the next timestep"""
        # Get the next TS from the atom trajectory
        at_ts = self._t.next()

        self._fill_ts(at_ts)

        return self.ts

    def _read_frame(self, frame):
        """Return a single frame"""
        at_ts = self._t[frame]

        self._fill_ts(at_ts)

        return self.ts

    def _fill_ts(self, other_ts):
        """Rip information from atomistic TS into our ts

        Make positions based on COM of our beads.
        """
        self.ts.frame = other_ts.frame
        self.ts._unitcell = other_ts._unitcell

        other_ts.positions *= self._weights[:, None]
        np.add.at(self.ts._pos, self._mapping, other_ts.positions)

    def _reopen(self):
        self._t._reopen()

    @property
    def n_frames(self):
        return len(self._t)

    def __repr__(self):
        return "<CG Trajectory doing {} beads >".format(len(self.beads))


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
        atu = self.atu = mda.Universe(*args, **kwargs)

        # Coarse grained Universe
        # Make a blank Universe for myself.
        super(CGUniverse, self).__init__()

        # Set up CG Topology
        nbeads = 0
        nbead_res = 0
        # Bead topology
        bead_resnames = []
        bead_atomresindex = []
        bead_masses = []
        bead_atomnames = []
        # for traj mapping
        atom_to_bead = np.zeros(len(atu.atoms), dtype=int)
        weights = np.zeros(len(atu.atoms), dtype=float)
        # Loop over residues in atomistic system and determine beads
        for residue in atu.residues:
            resmap = mapping[residue.resname]
            bead_resnames.append(residue.resname)

            for i, bead in enumerate(resmap):
                atom_to_bead[residue.atoms.ix[bead]] = nbeads
                nbeads += 1  # increment after as we want 0-based
                # mass weighting of positions
                masses = residue.atoms.masses[bead]
                weights[residue.atoms.ix[bead]] = masses / masses.sum()
                # bead attributes
                bead_masses.append(masses.sum())
                bead_atomnames.append('{}_{}'.format(residue.resname, i+1))
                bead_atomresindex.append(nbead_res)
            nbead_res += 1  # post increment again

        self._topology = topology.Topology(nbeads, nbead_res, 1,
                                           atom_resindex=bead_atomresindex,
                                           attrs=[
                                               ta.Masses(bead_masses),
                                               ta.Atomnames(bead_atomnames),
                                               ta.Resnames(bead_resnames),
        ])
        self.filename = atu.filename
        self._generate_from_topology()

        # This replaces load_new in a traditional Universe
        self.trajectory = CGTraj(self.atu.trajectory, atom_to_bead, weights)

    def __repr__(self):
        return "<CG Universe with {} beads>".format(len(self.atoms))

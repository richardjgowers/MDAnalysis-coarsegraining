# coarse graining in MDAnalysis
# Copyright (c) 2015 Richard J Gowers
# Released under the GNU Lesser General Public License, version 2 or later.

import MDAnalysis
from MDAnalysis.coordinates import base


class CGTraj(base.Reader):
    """Fakes a coarse grained trajectory object

    Takes an atomistic trajectory and list of beads and manipulates both
    to recreate a reader for a coarse grained Universe.

    Would also probably work as a standalone thingy for writing out
    coarse grained trajectories.
    """
    def __init__(self, trajectory, beads):
        """
        Arguments:
            universe - the atomistic Universe you start with
            mapping - list of list of indices
        """
        self._t = trajectory
        self.beads = beads

        self.numatoms = len(self.beads)
        self.convert_units = MDAnalysis.core.flags['convert_lengths']

        self.fixed = False
        self.periodic = True
        self.skip = 1

        self.ts = base.Timestep(self.numatoms)

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
        self.ts._pos[:] = [b.atoms.centerOfMass() for b in self.beads]

    def _reopen(self):
        # Rewind my reference trajectory
        self._t[0]

    def __iter__(self):
        self._reopen()
        while True:
            try:
                yield self._read_next_timestep()
            except StopIteration:
                self.rewind()
                raise StopIteration

    def __len__(self):
        return len(self._u.trajectory)

    def __repr__(self):
        return "<CG Trajectory doing {} beads >".format(len(self.beads))

import MDAnalysis as mda
from MDAnalysis.core.AtomGroup import AtomGroup, Atom

class Bead(Atom):
    def __init__(self, atoms):
        self.atoms = atoms

    def position(self):
        return self.atoms.centerOfMass()

u = mda.Universe('out.gro')

mapping = ([0, 1, 2],
           [3, 4, 5],
           [6, 7, 8])

beads = [AtomGroup(u.atoms[m]) for m in mapping]

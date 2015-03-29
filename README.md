# MDAnalysis-coarsegraining
An attempt at making a special subtype of MDAnalysis.Universe which allows a mapping scheme to be applied on the fly to an atomistic trajectory.

Currently, load an atomistic set of results as normal, but supply a mapping scheme (currently a list of lists of indices).
This should then make the CGUniverse supply positions (later velocities and forces?) according to the mapping scheme.

Potential useful for creating reference behaviour of atomistic trajectory when constructing coarse-grained models.

Current state is very WIP.

Requires MDAnalysis.

## Simple example, three water molecules coarse-grained into 3 sites:

mapping = ([0, 1, 2],
           [3, 4, 5],
           [6, 7, 8])

# Make a CGUniverse object.
# Designed to mimic a MDAnalysis Universe, but supplies information on beads not atoms
cgu = CGUniverse('new.gro', mapping=mapping)

# This is now a list of the beads
list(cgu.atoms)

# Center of mass for each water molecule, ie 3 coordinates
cgu.atoms.positions
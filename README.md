# MDAnalysis-coarsegraining
An attempt at making a special subtype of MDAnalysis.Universe which allows a mapping scheme to be applied on the fly to an atomistic trajectory.

Currently, load an atomistic set of results as normal, but supply a mapping scheme.
The mapping scheme is a dictionary which links residue names to a list of list of indices for each bead.
This should then make the CGUniverse supply positions (later velocities and forces?) according to the mapping scheme.

Potential useful for creating reference behaviour of atomistic trajectory when constructing coarse-grained models.

Current state is very WIP.

Requires MDAnalysis.

Simple example, three water molecules coarse-grained into 3 sites:
```
mapping = {'SOL':[[0, 1, 2]]}
```

Make a CGUniverse object.
Designed to mimic a MDAnalysis Universe, but supplies information on beads not atoms
```
cgu = CGUniverse('new.gro', mapping=mapping)
```

The atoms attribute of the Universe now provides information on the beads.
```
list(cgu.atoms)
```

Center of mass for each water molecule.
```
cgu.atoms.positions
```

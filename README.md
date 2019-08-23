# MDAnalysis-coarsegraining

A special subclass of MDAnalysis.Universe which allows a mapping scheme to be applied on the fly to an atomistic trajectory,
giving a virtual coarse-grained representation.

Currently, load an atomistic set of results as normal, but supply a mapping scheme.
The mapping scheme is a dictionary which links residue names to a list of list of indices for each bead.
This should then make the CGUniverse supply positions (later velocities and forces?) according to the mapping scheme.

Potential useful for creating reference behaviour of atomistic trajectory when constructing coarse-grained models.

Requires MDAnalysis v0.19+

Simple example, three water molecules coarse-grained into 3 sites:
```python
from cguniverse import CGUniverse

# for each resname, list of lists, each sublist becomes a bead.
# the indices refer to 
mapping = {'SOL':[[0, 1, 2]]}
cgu = CGUniverse('new.gro', mapping=mapping)

# The CGU has 3 "atoms", each representing a molecule
list(cgu.atoms)

# Each "atom" is a bead which gives the center of mass of each water molecule
cgu.atoms.positions
```

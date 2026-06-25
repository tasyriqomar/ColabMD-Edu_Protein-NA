from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename='RBD_clean.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.0)
PDBFile.writeFile(fixer.topology, fixer.positions, open('fixed.pdb', 'w'))


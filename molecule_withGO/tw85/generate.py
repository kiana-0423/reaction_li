from rdkit import Chem
from rdkit.Chem import AllChem
from ase.io import read, write


mol = Chem.MolFromSmiles('CC(/C=C/CCCCCCCCCCCC)CCCCCCCCCCCC')  # dodecane
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
Chem.MolToMolFile(mol, 'tw85.mol')
ase_mol = read('tw85.mol')
ase_mol.write('4.pdb')
from rdkit import Chem
from rdkit.Chem import AllChem
from ase.io import read, write


mol = Chem.MolFromSmiles('C1(OCCO)C(OCCO)C(C(COCCC(=O)OCC)OCCO)CC1') 
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
Chem.MolToMolFile(mol, 'tw80.mol')
ase_mol = read('tw80.mol')
ase_mol.write('tw80.pdb')
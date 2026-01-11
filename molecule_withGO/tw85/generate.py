from rdkit import Chem
from rdkit.Chem import AllChem
from ase.io import read, write


mol = Chem.MolFromSmiles('C1(OCC(=O)OCC)C(OCCC(OCC)=O)C(C(OCCO)OCCC(OCC)=O)CC1')  # dodecane
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
Chem.MolToMolFile(mol, 'tw85.mol')
ase_mol = read('tw85.mol')
ase_mol.write('tw85.pdb')
from ase.build import graphene, make_supercell
from ase.io import write

# 1) 原始六方 graphene
atoms = graphene(vacuum=18.0)  # 先给足真空

# 2) 变成正交胞（4 个 C/胞）
P = [[0, 1, 0],
     [-2, 1, 0],
     [0, 0, 1]]
atoms = make_supercell(atoms, P)

# 3) 再按需要扩胞（仍保持正交）
atoms *= (8, 8, 1)

# 4) 居中真空
atoms.center(axis=2)

write("POSCAR", atoms, format="vasp")
write("GO_rect.xyz", atoms)

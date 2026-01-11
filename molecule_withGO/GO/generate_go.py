import random
import numpy as np
from ase.build import graphene, make_supercell
from ase.neighborlist import NeighborList
from ase.io import write
from ase import Atom

# ----------------------------
# User controls
# ----------------------------
seed = 0
C_target = 64
C_over_O = 1.63          # target C/O
epoxy_ratio = 0.60       # fraction of O atoms that are epoxy (rest are hydroxyl)
vacuum = 18.0            # Å
z_shift_O = 1.25         # Å, initial O height above/below plane
z_shift_H = 2.15         # Å, initial H height for OH

random.seed(seed)

# ----------------------------
# 1) Build orthorhombic graphene cell (4 C atoms)
#    Transformation: A=a2, B=a2-2a1 (orthogonal)
# ----------------------------
atoms = graphene(vacuum=vacuum)

P = [[ 1,  1, 0],
     [-1,  1, 0],
     [ 0,  0, 1]]
atoms = make_supercell(atoms, P)
# ----------------------------
# 2) Expand to reach ~64 C
#    Each rectangular cell has 4 C, so (4,4,1) => 64 C
# ----------------------------
n = int(round((C_target / 4) ** 0.5))  # for 64 => 4
atoms *= (n, n, 1)
atoms.center(axis=2)

C_indices = [i for i, a in enumerate(atoms) if a.symbol == "C"]
C = len(C_indices)
if C != 64:
    print(f"Warning: got C={C}, expected ~{C_target}. You may adjust n.")

# ----------------------------
# 3) Decide O count from C/O
# ----------------------------
O_target = int(round(C / C_over_O))  # for C=64 => 39
# epoxy and hydroxyl counts
e = int(round(O_target * epoxy_ratio))
h = O_target - e

# sanity: need 2e + h <= C
# if violated, reduce epoxy until feasible
while 2 * e + h > C and e > 0:
    e -= 1
    h = O_target - e

print(f"C={C}, O_target={O_target}, epoxy={e}, hydroxyl={h}, C/O={C/O_target:.3f}, usedC={2*e+h}")

# ----------------------------
# 4) Neighbor list to get C-C bonds (for epoxy)
# ----------------------------
cutoffs = [1.75] * len(atoms)  # >1.42 Å C-C
nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
nl.update(atoms)

C_set = set(C_indices)

# collect unique nearest-neighbor C-C bonds
cc_bonds = []
for i in C_indices:
    neigh, offsets = nl.get_neighbors(i)
    for j, off in zip(neigh, offsets):
        if j in C_set and i < j:
            d = atoms.get_distance(i, j, mic=True)
            if d < 1.65:
                cc_bonds.append((i, j))

random.shuffle(cc_bonds)

# choose epoxy bonds without reusing C
used_C = set()
epoxy_bonds = []
for i, j in cc_bonds:
    if i not in used_C and j not in used_C:
        epoxy_bonds.append((i, j))
        used_C.add(i)
        used_C.add(j)
    if len(epoxy_bonds) == e:
        break

if len(epoxy_bonds) < e:
    raise RuntimeError("Not enough available C-C bonds for the requested epoxy count. Reduce epoxy_ratio or enlarge cell.")

# choose OH sites from remaining C
remaining_C = [i for i in C_indices if i not in used_C]
random.shuffle(remaining_C)
oh_sites = remaining_C[:h]
if len(oh_sites) < h:
    raise RuntimeError("Not enough remaining C sites for OH. Reduce O_target or enlarge cell.")

# ----------------------------
# 5) Add O/H atoms
# ----------------------------
def add_atom(sym, pos):
    atoms.append(Atom(sym, pos))

side = 1  # alternate +z / -z to reduce curvature

# Epoxy: O near bond midpoint
for i, j in epoxy_bonds:
    ri = atoms[i].position
    rj = atoms[j].position
    mid = 0.5 * (ri + rj)
    add_atom("O", [mid[0], mid[1], mid[2] + side * z_shift_O])
    side *= -1

# Hydroxyl: O above C, H further out
for i in oh_sites:
    rc = atoms[i].position
    o_pos = [rc[0], rc[1], rc[2] + side * z_shift_O]
    h_pos = [rc[0], rc[1], rc[2] + side * z_shift_H]
    add_atom("O", o_pos)
    add_atom("H", h_pos)
    side *= -1

# ----------------------------
# 6) Output
# ----------------------------
write("GO_C64_O{}_seed{}.xyz".format(O_target, seed), atoms)
write("POSCAR_GO_C64", atoms, format="vasp")

print("Wrote: POSCAR_GO_C64 and GO_C64_O{}_seed{}.xyz".format(O_target, seed))

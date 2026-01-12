import numpy as np
from ase import Atoms
from ase.io import write

# ====== user settings ======
N_Fe, N_O = 40, 0
rho = 5.7  # g/cm^3 (you can change)
seed = 1

# minimum distance matrix (Å) - can be tuned
dmin = {
    ('Fe','Fe'): 2.0,
    ('Fe','O') : 1.7,
    ('O','O')  : 1.5
}

# ====== constants ======
NA = 6.02214076e23
M_Fe = 55.845  # g/mol
M_O  = 15.999  # g/mol

mass_g = (N_Fe*M_Fe + N_O*M_O)/NA
V_cm3 = mass_g / rho
V_A3 = V_cm3 * 1e24
L = V_A3 ** (1/3)

rng = np.random.default_rng(seed)

symbols = ['Fe']*N_Fe + ['O']*N_O
pos = []

def ok(new_sym, new_pos, symbols_exist, pos_exist):
    for s, p in zip(symbols_exist, pos_exist):
        key = (new_sym, s) if (new_sym, s) in dmin else (s, new_sym)
        if np.linalg.norm(new_pos - p) < dmin[key]:
            return False
    return True

for sym in symbols:
    for _ in range(200000):
        trial = rng.random(3) * L
        if ok(sym, trial, symbols[:len(pos)], pos):
            pos.append(trial)
            break
    else:
        raise RuntimeError("Failed to place atoms. Try larger L or smaller dmin.")

atoms = Atoms(symbols=symbols, positions=np.array(pos), cell=[L, L, L], pbc=True)

# outputs
write("Fe40O20.xyz", atoms)
write("POSCAR", atoms, format="vasp")
write("Fe40O20.data", atoms, format="lammps-data")
print("Box length L (Å) =", L)

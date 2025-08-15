import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, ket2dm, mesolve, Qobj, Options

# ---------- Constants & units ----------
hbar_meV_ps = 0.6582119514     # ħ in meV·ps
k_B_meV_per_K = 0.08617333262  # k_B in meV/K
sim_unit_ps = 0.659            # 1 simulation time unit = 0.659 ps
opts = Options(progress_bar=False)

# ---------- Hamiltonian parameters ----------
# Site energies relative to the mean (meV): Δ = 59.28 meV
eps_D_meV = -29.64
eps_A_meV = +29.64
Delta_meV = eps_A_meV - eps_D_meV

# Electronic coupling (meV) — can be ±; sign flips which eigenstate is lower
J_meV = -34.0

# Convert to simulation units (dimensionless frequencies)
# H_site = [[-Δ/2, J], [J, +Δ/2]]  (units: meV)
scale = sim_unit_ps / hbar_meV_ps
Delta_sim = Delta_meV * scale
J_sim = J_meV * scale

# ---------- Bath / thermal parameters ----------
lambda_fast_meV = 50.0 / 8.065544  # 50 cm^-1 -> meV
tau_c_ps = 1.0
T_K = 300.0

# ---------- Basis states (site basis) ----------
ket1 = basis(2, 0)   # |1>
ket2 = basis(2, 1)   # |2>

# ---------- Hamiltonian in site basis (simulation units) ----------
H_sim = Qobj([[-0.5*Delta_sim,  J_sim],
              [ J_sim,          +0.5*Delta_sim]])

# ---------- Diagonalize to get energy eigenstates ----------
evals, evecs = H_sim.eigenstates()

# We'll keep both eigenstates, but DO NOT sort by energy for labeling.
# Instead, we determine bright/dark by transition dipole strength.
# Assume equal, parallel site transition dipoles (you can edit mu1, mu2 to change geometry):
mu1 = np.array([1.0, 0.0, 0.0])  # arbitrary units
mu2 = np.array([1.0, 0.0, 0.0])  # parallel & equal to mu1

def dipole_strength(evec):
    """
    Compute |μ_total|^2 for an eigenstate in the site basis:
    μ_total = c1 * μ1 + c2 * μ2, where evec = c1|1> + c2|2>.
    Minimal-change fix: read c1,c2 directly from evec.full().
    """
    vec = evec.full().flatten()
    c1 = vec[0]
    c2 = vec[1]
    mu_tot = c1 * mu1 + c2 * mu2
    return float(np.vdot(mu_tot, mu_tot))  # complex-safe magnitude^2

# Compute strengths for both eigenstates as returned
strength_0 = dipole_strength(evecs[0])
strength_1 = dipole_strength(evecs[1])

# Assign bright/dark by dipole strength (physically meaningful, J-like labeling)
if strength_0 >= strength_1:
    ket_bright = evecs[0]
    ket_dark   = evecs[1]
    E_bright_sim = float(evals[0])
    E_dark_sim   = float(evals[1])
else:
    ket_bright = evecs[1]
    ket_dark   = evecs[0]
    E_bright_sim = float(evals[1])
    E_dark_sim   = float(evals[0])

# ---------- Exciton gap in meV (independent of label convention) ----------
exciton_split_meV = np.sqrt(Delta_meV**2 + 4.0*J_meV**2)
omega0_ps_inv = exciton_split_meV / hbar_meV_ps  # transition frequency (ps^-1)

# ---------- Thermal transfer rates (between energy eigenstates) ----------
k_c_ps = 1.0 / tau_c_ps

# Classical Drude–Lorentz spectrum at ω0 (population-only toy model)
S_cl_ps_inv = (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2) \
              * (k_c_ps / (k_c_ps**2 + omega0_ps_inv**2))

beta = 1.0 / (k_B_meV_per_K * T_K)
n_BE = 1.0 / (np.exp(beta * exciton_split_meV) - 1.0)

# Down = higher -> lower; Up = lower -> higher, by detailed balance.
k_down_ps = S_cl_ps_inv * (n_BE + 1.0)
k_up_ps   = S_cl_ps_inv * n_BE

# Convert to simulation units
k_down_sim = k_down_ps * sim_unit_ps
k_up_sim   = k_up_ps   * sim_unit_ps

# ---------- Collapse operators in the energy basis (bright/dark labels) ----------
# We still need to send population from the higher-energy eigenstate to the lower-energy one.
# Determine which labeled state is higher by comparing their eigenvalues:
if E_bright_sim > E_dark_sim:
    # bright is higher
    c_ops = [
        np.sqrt(k_down_sim) * ket_dark   * ket_bright.dag(),  # bright -> dark
        np.sqrt(k_up_sim)   * ket_bright * ket_dark.dag()     # dark   -> bright
    ]
else:
    # dark is higher
    c_ops = [
        np.sqrt(k_down_sim) * ket_bright * ket_dark.dag(),    # dark   -> bright
        np.sqrt(k_up_sim)   * ket_dark   * ket_bright.dag()   # bright -> dark
    ]

# ---------- Initial state: 50/50 mixed in the energy basis ----------
rho0 = 0.5 * ket2dm(ket_bright) + 0.5 * ket2dm(ket_dark)

# ---------- Time grids (coarse-grained over sub-ps dynamics) ----------
# dt_sim = 1 → dt_ps ≈ 0.659 ps  (≥ 0.5 ps as requested)
tlist_sim = np.arange(0.0, 400.0, 1.0)
tlist_ps  = tlist_sim * sim_unit_ps

# ---------- E-ops: populations in the bright/dark energy basis ----------
e_ops = {'p_bright': ket2dm(ket_bright), 'p_dark': ket2dm(ket_dark)}

# ---------- Evolve ----------
res = mesolve(H_sim, rho0, tlist_sim, c_ops=c_ops, e_ops=e_ops, options=opts)
pBB = np.real(res.expect[0])
pDD = np.real(res.expect[1])

# ---------- Plot (energy basis: bright/dark labels preserved) ----------
plt.figure(figsize=(7.5, 4.8))
plt.plot(tlist_ps, pBB, label=r'$\rho_{++}$', color='dodgerblue', linestyle='dashed')
plt.plot(tlist_ps, pDD, label=r'$\rho_{--}$', color='black', linestyle='dashed')
plt.xlabel("Time (ps)", fontsize=24)
plt.ylabel("Population", fontsize=24)
plt.tick_params(labelsize=18)
plt.ylim(-0.02, 1.02)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=18)
plt.tight_layout()
plt.show()

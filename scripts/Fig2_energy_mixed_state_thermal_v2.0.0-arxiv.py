
import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, ket2dm, mesolve, Qobj, Options

# ---------- Constants & units ----------
hbar_meV_ps = 0.6582119514     # ħ in meV·ps
k_B_meV_per_K = 0.08617333262  # k_B in meV/K
def cm_to_meV(x_cm):
    return x_cm / 8.065544
LAMBDA_FAST_CM = 270.0
LAMBDA_FAST_MEV = cm_to_meV(LAMBDA_FAST_CM)
sim_unit_ps = 0.659            # 1 simulation time unit = 0.659 ps
opts = Options(progress_bar=False, atol=1e-10, rtol=1e-8, nsteps=10000)


# --- Dephasing toggle & bath param for site-basis dephasing ---
USE_SITE_DEPHASING = False  # set False for transfer-only toy model
tau_c_ps = 0.1
# ---------- Pure-dephasing rate (classical OU / Drude–Lorentz) ----------
def gamma_phi(lambda_fast_meV, tau_c_ps, T_K):
    k_c = 1.0 / tau_c_ps
    return (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2 * k_c)
# ---------- Bose–Einstein occupancy with numerical guards ----------
def bose_einstein(DeltaE_meV, T_K):
    if T_K <= 0:
        return 0.0
    x = DeltaE_meV / (k_B_meV_per_K * T_K)
    if x > 50.0:
        return 0.0
    if x < 1e-6:
        return 1.0 / x
    return 1.0 / np.expm1(x)

# ---------- Classical Drude–Lorentz transfer prefactor + detailed balance ----------
def energy_transfer_rates(lambda_fast_meV, tau_c_ps, T_K, DeltaE_meV):
    k_c = 1.0 / tau_c_ps
    omega0_ps_inv = DeltaE_meV / hbar_meV_ps
    S_cl = (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2) * (k_c / (k_c**2 + omega0_ps_inv**2))
    n_BE = bose_einstein(DeltaE_meV, T_K)
    k_down_ps = S_cl * (n_BE + 1.0)
    k_up_ps   = S_cl * n_BE
    return k_down_ps, k_up_ps

# ---------- System parameters (heterodimer) ----------
# Site energies relative to the mean (meV): Δ = 59 meV
eps_D_meV = -29.5
eps_A_meV = +29.5
Delta_meV = eps_A_meV - eps_D_meV
J_meV = -34.0 # Electronic coupling (meV) — sign flips which exciton is lower

# Convert to simulation units (dimensionless frequencies)
scale = sim_unit_ps / hbar_meV_ps
Delta_sim = Delta_meV * scale
J_sim = J_meV * scale

# ---------- Hamiltonian in site basis (simulation units) ----------
# H_site = [[-Δ/2, J], [J, +Δ/2]]
H_sim = Qobj([[-0.5*Delta_sim,  J_sim],
              [ J_sim,          +0.5*Delta_sim]])

# ---------- Bath / thermal parameters ----------
lambda_fast_meV = LAMBDA_FAST_MEV  # 270 cm^-1 -> meV
tau_c_ps = 0.1                     # correlation time of fast bath (ps)
T_K = 293.0                        # temperature in Kelvin

# ---------- Basis states (site basis) ----------
ket1 = basis(2, 0)   # |1>
ket2 = basis(2, 1)   # |2>

# ---------- Diagonalize to get energy eigenstates ----------
evals, evecs = H_sim.eigenstates()

# Determine bright/dark by relative energy (bright = lower energy)
if float(evals[0]) <= float(evals[1]):
    ket_plus, ket_minus = evecs[0], evecs[1]
    E_bright_sim, E_dark_sim = float(evals[0]), float(evals[1])
else:
    ket_plus, ket_minus = evecs[1], evecs[0]
    E_bright_sim, E_dark_sim = float(evals[1]), float(evals[0])

# ---------- Exciton gap in meV (independent of label convention) ----------

exciton_split_meV = np.sqrt(Delta_meV**2 + 4.0*J_meV**2)
omega0_ps_inv = exciton_split_meV / hbar_meV_ps  # ps^-1

# ---------- Thermal transfer rates (between energy eigenstates) ----------
k_c_ps = 1.0 / tau_c_ps

# Classical Drude–Lorentz spectrum at ω0 (population-only toy model)
S_cl_ps_inv = (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2) \
              * (k_c_ps / (k_c_ps**2 + omega0_ps_inv**2))

beta = 1.0 / (k_B_meV_per_K * T_K)
n_BE = 1.0 / (np.exp(beta * exciton_split_meV) - 1.0)

# Down = higher -> lower; Up = lower -> higher (detailed balance)
k_down_ps = S_cl_ps_inv * (n_BE + 1.0)
k_up_ps   = S_cl_ps_inv * n_BE

# Convert to simulation units
k_down_sim = k_down_ps * sim_unit_ps
k_up_sim   = k_up_ps   * sim_unit_ps

# ---------- Collapse operators in the energy basis (bright/dark) ----------
# Send population from the higher-energy eigenstate to the lower-energy one.
if E_bright_sim > E_dark_sim:
    # bright is higher
    c_ops = [
        np.sqrt(max(k_down_sim, 0.0)) * ket_minus   * ket_plus.dag(),  # bright -> dark
        np.sqrt(max(k_up_sim,   0.0)) * ket_plus * ket_minus.dag()     # dark   -> bright
    ]
else:
    # dark is higher
    c_ops = [
        np.sqrt(max(k_down_sim, 0.0)) * ket_plus * ket_minus.dag(),    # dark   -> bright
        np.sqrt(max(k_up_sim,   0.0)) * ket_minus   * ket_plus.dag()   # bright -> dark
    ]

# ---------- Initial state: 50/50 mixed in the energy basis ----------
rho0 = 0.5 * ket2dm(ket_plus) + 0.5 * ket2dm(ket_minus)

# --- Grid: heterodimer-based Δt, Nspp = 200 ---
Nspp = 200
DeltaE_ref_meV = exciton_split_meV
Tosc_ps = 2.0 * np.pi * hbar_meV_ps / DeltaE_ref_meV
dt_ps   = Tosc_ps / Nspp
dt_sim  = dt_ps / sim_unit_ps
print(f"[grid:Fig2] Tosc={Tosc_ps*1e3:.2f} fs, Δt={dt_ps*1e3:.3f} fs (heterodimer)")

# --- Time grid (long window, same as original: 5 ps) ---
T_long_ps = 5.0
tlist_sim = np.arange(0.0, T_long_ps/sim_unit_ps + 0.5*dt_sim, dt_sim)
tlist_ps  = tlist_sim * sim_unit_ps

# ---------- E-ops: populations in the bright/dark energy basis ----------
# Optional site-basis pure dephasing
if USE_SITE_DEPHASING:
    gphi_ps  = gamma_phi(LAMBDA_FAST_MEV, tau_c_ps, T_K)
    gphi_sim = gphi_ps * sim_unit_ps
    ket1, ket2 = basis(2,0), basis(2,1)
    c_ops += [np.sqrt(max(gphi_sim, 0.0)) * ket1 * ket1.dag(),
              np.sqrt(max(gphi_sim, 0.0)) * ket2 * ket2.dag()]

e_ops = {'p_bright': ket2dm(ket_plus), 'p_dark': ket2dm(ket_minus)}

# ---------- Evolve ----------
res = mesolve(H_sim, rho0, tlist_sim, c_ops=c_ops, e_ops=e_ops, options=opts)
pBB = np.real(res.expect[0])
pDD = np.real(res.expect[1])

# ---------- Plot ----------
plt.figure(figsize=(7.0, 4.8))
plt.plot(tlist_ps, pBB, label=r'$\rho_{++}$', color='dodgerblue')
plt.plot(tlist_ps, pDD, label=r'$\rho_{--}$', color='black')
plt.xlabel("Time (ps)", fontsize=12)
plt.ylabel("Population", fontsize=12)
plt.ylim(-0.02, 1.02)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

# -------- Git commit hash helper (best-effort) --------
def _get_git_commit_short():
    try:
        import subprocess
        h = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.STDOUT)
        return h.decode().strip()
    except Exception:
        return "unknown"

_commit_hash = _get_git_commit_short()
print("Git commit:", _commit_hash)

plt.savefig("Fig2_energy_mixed_state_thermal_v2.0.0-arxiv.png", dpi=300, bbox_inches="tight")
print("Saved: Fig2_energy_mixed_state_thermal_v2.0.0-arxiv.png")


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
sim_unit_ps = 0.659            # 1 simulation unit = 0.659 ps
opts = Options(progress_bar=False, atol=1e-10, rtol=1e-8, nsteps=10000)

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

# ---------- System parameters (homodimer) ----------
J_meV = -34.0                  # coupling in meV

# ---------- Bath / thermal parameters ----------
lambda_fast_meV = cm_to_meV(LAMBDA_FAST_CM) # 270 cm^-1 -> meV
tau_c_ps = 0.1
T_K = 293.0

# ---------- Dephasing model (pure dephasing on sites) ----------
def gamma_phi(lambda_fast_meV, tau_c_ps, T_K):
    k_c = 1.0 / tau_c_ps  # ps^-1
    return (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2 * k_c)  # ps^-1
gamma_phi_ps_inv = gamma_phi(lambda_fast_meV, tau_c_ps, T_K)  # ps^-1

# ---------- Convert to simulation units ----------
J_sim = J_meV * (sim_unit_ps / hbar_meV_ps)     # dimensionless
gamma_sim = gamma_phi_ps_inv * sim_unit_ps      # 1/sim_unit

# ---------- Hamiltonian in site basis (sim units) ----------
H_sim = Qobj([[0.0, J_sim],
              [J_sim, 0.0]])

# --- Grid: homodimer-based Δt (Δ = 0), Nspp = 200 ---
Nspp = 200
DeltaE_ref_meV = 2.0 * abs(J_meV)                         # fastest oscillation for homodimer: ΔE = 2|J|
Tosc_ps        = 2.0 * np.pi * hbar_meV_ps / DeltaE_ref_meV
dt_ps          = Tosc_ps / Nspp
dt_fs          = dt_ps * 1000.0
dt_sim         = dt_ps / sim_unit_ps

T_win_fs  = 25.0                                          # window in fs
tlist_fs  = np.arange(0.0, T_win_fs + 0.5*dt_fs, dt_fs)    # labels
tlist_sim = (tlist_fs / 1000.0) / sim_unit_ps              # mesolve grid

# ---------- Basis states ----------
ket1 = basis(2, 0)   # |1>
ket2 = basis(2, 1)   # |2>
ket_plus  = (ket1 + ket2).unit()           # |+>
ket_minus = (ket1 - ket2).unit()           # |->

# ---------- Collapse operators: site dephasing + thermal energy transfer ----------
# Site-basis pure dephasing
c_ops = [
    np.sqrt(gamma_sim) * ket1 * ket1.dag(),
    np.sqrt(gamma_sim) * ket2 * ket2.dag()
]

# Thermal |+> <-> |-> population transfer
omega_meV = 2.0 * abs(J_meV)                         # transition energy ΔE (meV)
omega_ps_inv = omega_meV / hbar_meV_ps               # ps^-1
k_c_ps = 1.0 / tau_c_ps                              # ps^-1
S_cl_ps_inv = (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2) \
      * (k_c_ps / (k_c_ps**2 + omega_ps_inv**2))
beta = 1.0 / (k_B_meV_per_K * T_K)
n_th = 1.0 / (np.exp(beta * omega_meV) - 1.0)
gamma_down_sim = (S_cl_ps_inv * (n_th + 1.0)) * sim_unit_ps
gamma_up_sim   = (S_cl_ps_inv * n_th)               * sim_unit_ps
C_down = np.sqrt(max(gamma_down_sim, 0.0)) * (ket_plus  * ket_minus.dag())
C_up   = np.sqrt(max(gamma_up_sim,   0.0)) * (ket_minus * ket_plus.dag())

c_ops.extend([C_down, C_up])

# ---------- Helper ----------
def pops_vs_time(H, psi0, meas0, meas1):
    e_ops = {'p0': ket2dm(meas0), 'p1': ket2dm(meas1)}
    res = mesolve(H, psi0, tlist_sim, c_ops=c_ops, e_ops=e_ops, options=opts)
    return np.real(res.expect[0]), np.real(res.expect[1])

# ---------- Single (ENERGY basis) panel; Initial state: ψ0 = |+> ----------
pPP, pMM = pops_vs_time(H_sim, ket_plus, ket_plus, ket_minus)

# ---------- Plot ----------
fig, ax = plt.subplots(1, 1, figsize=(7.0, 4.8))
ax.plot(tlist_fs, pPP, label=r'$\rho_{++}$', color='dodgerblue')
ax.plot(tlist_fs, pMM, label=r'$\rho_{--}$', color='black')
ax.set_xlabel("Time (fs)", fontsize=12)
ax.set_ylabel("Population", fontsize=12)
ax.set_ylim(-0.02, 1.02)
ax.grid(True, alpha=0.3)
ax.legend()
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

plt.savefig("Fig1_energy_basis_thermal_dephasing_v2.0.0-arxiv.png", dpi=300, bbox_inches="tight")
print("Saved: Fig1_energy_basis_thermal_dephasing_v2.0.0-arxiv.png")

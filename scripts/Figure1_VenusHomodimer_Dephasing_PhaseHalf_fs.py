import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, ket2dm, mesolve, Qobj, Options

# ---------- Constants & unit settings ----------
hbar_meV_ps = 0.6582119514     # ħ in meV·ps
k_B_meV_per_K = 0.08617333262  # k_B in meV/K
J_meV = -34.0                   # coupling in meV
sim_unit_ps = 0.659             # 1 simulation unit = 0.659 ps
opts = Options(progress_bar=False)

# ---------- Dephasing model (user-provided formula) ----------
def gamma_phi(lambda_fast_meV, tau_c_ps, T_K):
    k_c = 1.0 / tau_c_ps  # ps^-1
    return (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2 * k_c)  # ps^-1

# Example parameters (Venus-like)
lambda_fast_meV = 20.0 / 8.065544  # 20 cm^-1 -> meV
tau_c_ps = 1.0
T_K = 293.0
gamma_phi_ps_inv = gamma_phi(lambda_fast_meV, tau_c_ps, T_K)  # ps^-1
T2_ps = 1.0 / gamma_phi_ps_inv if gamma_phi_ps_inv > 0 else np.inf

# ---------- Convert to simulation units ----------
J_sim = J_meV * (sim_unit_ps / hbar_meV_ps)     # dimensionless
gamma_sim = gamma_phi_ps_inv * sim_unit_ps      # 1/sim_unit

# ---------- Time grid (unchanged) ----------
tlist_sim = np.linspace(0.0, 0.2, 9000)   # simulation units
tlist_fs  = tlist_sim * sim_unit_ps * 1000  # fs for plotting

# ---------- Basis states ----------
ket1 = basis(2, 0)   # |1>
ket2 = basis(2, 1)   # |2>
ket_plus  = (ket1 + ket2).unit()           # |+>
ket_minus = (ket1 - ket2).unit()           # |->
ket_phi   = (ket1 + 1j*ket2).unit()        # φ = π/2 site superposition

# ---------- Hamiltonian in site basis (sim units) ----------
H_sim = Qobj([[0.0, J_sim],
              [J_sim, 0.0]])

# ---------- Collapse operators for pure dephasing in SITE basis (unchanged) ----------
c_ops = [
    np.sqrt(gamma_sim) * ket1 * ket1.dag(),
    np.sqrt(gamma_sim) * ket2 * ket2.dag()
]

# ---------- Thermal |+> <-> |-> population transfer (match other scripts) ----------
omega_meV = 2.0 * abs(J_meV)                         # transition energy ΔE (meV)
omega_ps_inv = omega_meV / hbar_meV_ps               # ps^-1
k_c_ps = 1.0 / tau_c_ps                              # ps^-1

# Drude–Lorentz classical spectrum at ω0 (same as other scripts)
S_cl_ps_inv = (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2) \
      * (k_c_ps / (k_c_ps**2 + omega_ps_inv**2))

beta = 1.0 / (k_B_meV_per_K * T_K)
n_th = 1.0 / (np.exp(beta * omega_meV) - 1.0)

gamma_down_ps_inv = S_cl_ps_inv * (n_th + 1.0)       # higher -> lower
gamma_up_ps_inv   = S_cl_ps_inv * n_th               # lower  -> higher

# convert to simulation units
gamma_down_sim = gamma_down_ps_inv * sim_unit_ps
gamma_up_sim   = gamma_up_ps_inv   * sim_unit_ps

# jump operators in the ENERGY basis (|+>, |->)
C_down = np.sqrt(gamma_down_sim) * (ket_plus  * ket_minus.dag())
C_up   = np.sqrt(gamma_up_sim)   * (ket_minus * ket_plus.dag())

# rebuild c_ops with pure dephasing + thermal transfers
c_ops = [
    np.sqrt(gamma_sim) * ket1 * ket1.dag(),
    np.sqrt(gamma_sim) * ket2 * ket2.dag(),
    C_down, C_up
]
# ---------- END thermal block ----------

def populations_vs_time(H, psi0, meas0, meas1):
    e_ops = {'p0': ket2dm(meas0), 'p1': ket2dm(meas1)}
    res = mesolve(H, psi0, tlist_sim, c_ops=c_ops, e_ops=e_ops, options=opts)
    return np.real(res.expect[0]), np.real(res.expect[1])

# Left panel (SITE basis): initialize with φ = π/2 superposition (|1> + i|2>)/sqrt(2)
p11_left, p22_left = populations_vs_time(H_sim, ket_phi, ket1, ket2)

# Right panel (ENERGY basis): ψ0 = |+> (unchanged)
pPP_right, pMM_right = populations_vs_time(H_sim, ket_plus, ket_plus, ket_minus)

# ---------- Plot (unchanged) ----------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), sharex=True, sharey=True)

plt.xticks(np.arange(0, 140, 50))

# Left: Site basis
axes[0].plot(tlist_fs, p11_left, label=r'$\rho_{11}$', color='dodgerblue', linestyle='solid')
axes[0].plot(tlist_fs, p22_left, label=r'$\rho_{22}$', color='black', linestyle='solid')

# |rho_12| curve (kept identical: computed from |1>)
res_left = mesolve(H_sim, ket1, tlist_sim, c_ops=c_ops, e_ops=[ket1*ket2.dag()], options=opts)
rho12_abs_left = np.abs(res_left.expect[0])
axes[0].plot(tlist_fs, rho12_abs_left, label=r'$|\rho_{12}|$', color='red', linestyle='solid')

axes[0].set_xlabel("Time (fs)", fontsize=24)
axes[0].set_ylabel("Population / Coherence", fontsize=24)
axes[0].set_ylim(-0.05, 1.05)
axes[0].grid(True, alpha=0.3)
axes[0].legend(fontsize=18)
axes[0].tick_params(labelsize=18)

# Right: Energy basis
axes[1].plot(tlist_fs, pPP_right, label=r'$\rho_{++}$', color='dodgerblue', linestyle='dashed')
axes[1].plot(tlist_fs, pMM_right, label=r'$\rho_{--}$', color='black', linestyle='dashed')
axes[1].set_xlabel("Time (fs)", fontsize=24)
axes[1].grid(True, alpha=0.3)
axes[1].legend(fontsize=18)
axes[1].tick_params(labelsize=18)

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

# Stamp the figure subtly (bottom-right)
try:
    import matplotlib.pyplot as _plt
    _fig = _plt.gcf()
    _fig.text(0.995, 0.005, f"commit: {_commit_hash}", ha="right", va="bottom", fontsize=8, alpha=0.8)
except Exception:
    pass

    plt.savefig(r'/mnt/data/therm_compare_proper_preps_units_dephasing_VenusTrueHomodimer_fs_rho_12_phase_half_v1.1.0-arxiv.png', dpi=300, bbox_inches='tight')
    print('Saved:', r'/mnt/data/therm_compare_proper_preps_units_dephasing_VenusTrueHomodimer_fs_rho_12_phase_half_v1.1.0-arxiv.png')

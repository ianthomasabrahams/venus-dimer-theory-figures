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
tau_c_ps = 0.5
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

# ---------- ADD: Thermal state transfer in ENERGY basis (|+> <-> |->) ----------
# Detailed-balance rates at temperature T_K for transition energy ΔE = 2|J_meV|
T1_ps = 1.0                               # base relaxation timescale (feel free to tune)
gamma0_ps_inv = 1.0 / T1_ps               # ps^-1
omega_meV = 2.0 * abs(J_meV)              # meV
beta = 1.0 / (k_B_meV_per_K * T_K)        # 1/meV
n_th = 1.0 / (np.exp(omega_meV * beta) - 1.0)

gamma_down_ps_inv = gamma0_ps_inv * (n_th + 1.0)  # |-> → |+>
gamma_up_ps_inv   = gamma0_ps_inv * (n_th)        # |+> → |->

# convert to simulation units
gamma_down_sim = gamma_down_ps_inv * sim_unit_ps
gamma_up_sim   = gamma_up_ps_inv   * sim_unit_ps

# jump operators in the ENERGY basis
C_down = np.sqrt(gamma_down_sim) * (ket_plus  * ket_minus.dag())  # dark to bright
C_up   = np.sqrt(gamma_up_sim)   * (ket_minus * ket_plus.dag())   # bright to dark

# append to existing collapse ops
c_ops += [C_down, C_up]
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
plt.show()

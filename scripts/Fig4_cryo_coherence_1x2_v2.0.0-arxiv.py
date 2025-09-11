
# -----------------------------------------------------------------------------------
# 1x2 figure (single row) for a two-site excitonic dimer.
# Left panel: Im[rho_pm] at 250 mK (energy basis, init (|+> + i|->)/sqrt(2))
# Right panel: Im[rho_pm] at 400 mK (same initial state)
# -----------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, mesolve, Qobj, Options, expect

# ---------- Constants & units ----------
hbar_meV_ps = 0.6582119514
k_B_meV_per_K = 0.08617333262
def cm_to_meV(x_cm):
    return x_cm / 8.065544
LAMBDA_FAST_CM = 270.0
LAMBDA_FAST_MEV = cm_to_meV(LAMBDA_FAST_CM)
sim_unit_ps = 0.659
opts = Options(progress_bar=False, atol=1e-10, rtol=1e-8, nsteps=10000)

# ---------- Bose–Einstein occupancy with numerical guards ----------
def bose_einstein(DeltaE_meV, T_K):
    if T_K <= 0: return 0.0

    x = DeltaE_meV / (k_B_meV_per_K * T_K)
    if x > 50.0: return 0.0
    if x < 1e-6: return 1.0 / x
    return 1.0 / np.expm1(x)

# ---------- Dephasing model (pure dephasing on sites) ----------
def gamma_phi(lambda_fast_meV, tau_c_ps, T_K):
    gamma_c = 1.0 / tau_c_ps
    return (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2 * gamma_c)

# Thermal |+> <-> |-> population transfer
def energy_transfer_rates(lambda_fast_meV, tau_c_ps, T_K, DeltaE_meV):
    gamma_c = 1.0 / tau_c_ps
    omega0_ps_inv = DeltaE_meV / hbar_meV_ps
    S_cl = (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2) * (gamma_c / (gamma_c**2 + omega0_ps_inv**2))
    n_BE = bose_einstein(DeltaE_meV, T_K)
    gamma_down_ps = S_cl * (n_BE + 1.0)
    gamma_up_ps   = S_cl * n_BE
    return gamma_down_ps, gamma_up_ps

# ---------- System parameters (homodimer) ----------
J_meV = -34.0

# ---------- Hamiltonian (homodimer) ----------
J_sim = J_meV * (sim_unit_ps / hbar_meV_ps)
H_sim = Qobj([[0.0, J_sim], [J_sim, 0.0]])

# ---------- Bases ----------
ket1 = basis(2, 0)
ket2 = basis(2, 1)
ket_plus  = (ket1 + ket2).unit()
ket_minus = (ket1 - ket2).unit()
ket_phi = (ket_plus + 1j*ket_minus).unit()

# ---------- Collapse operators builder ----------
def make_c_ops(lambda_fast_meV, tau_c_ps, T_K):
    gphi_ps = gamma_phi(lambda_fast_meV, tau_c_ps, T_K)
    gphi_sim = gphi_ps * sim_unit_ps
    c_ops = [np.sqrt(max(gphi_sim, 0.0)) * ket1 * ket1.dag(),
             np.sqrt(max(gphi_sim, 0.0)) * ket2 * ket2.dag()]
    DeltaE_meV = 2.0 * abs(J_meV)
    gamma_down_ps, gamma_up_ps = energy_transfer_rates(lambda_fast_meV, tau_c_ps, T_K, DeltaE_meV)

    L_down = np.sqrt(max(gamma_down_ps * sim_unit_ps, 0.0)) * ket_plus * ket_minus.dag()
    L_up   = np.sqrt(max(gamma_up_ps   * sim_unit_ps, 0.0)) * ket_minus * ket_plus.dag()
    c_ops += [L_down, L_up]
    return c_ops, gphi_ps

# ---------- Runner (Initial state: ket_phi = (|+>+1j*|->).unit()) ----------
def run_case(lambda_fast_cm, tau_c_ps, T_K, t_ps):
    lam_meV = cm_to_meV(lambda_fast_cm)
    c_ops, gphi_ps = make_c_ops(lam_meV, tau_c_ps, T_K)
    t_sim = t_ps / sim_unit_ps
    res = mesolve(H_sim, ket_phi, t_sim, c_ops=c_ops, options=opts)
    op_pm = ket_plus * ket_minus.dag()
    rho_pm_imag = np.imag(expect(op_pm, res.states))
    return {'t_ps': t_ps, 'rho_pm_imag': np.array(rho_pm_imag),
            'gamma_phi_ps': gphi_ps, 'T_K': T_K}

# --- Grid: homodimer-based Δt (Δ=0), Nspp = 200 ---
DeltaE_ref_meV = 2.0 * abs(J_meV)                        # = 68 meV
Tosc_ps = 2.0 * np.pi * hbar_meV_ps / DeltaE_ref_meV
dt_ps   = Tosc_ps / 200.0                                # <-- FIXED: 200 samples per period
print(f"[grid:Fig4] Tosc={Tosc_ps*1e3:.2f} fs, Δt={dt_ps*1e3:.3f} fs (homodimer)")

# ---------- Cases ----------
tau_left_ps, lambda_left_cm, T_left_K = 0.1, 270.0, 0.25
tau_right_ps, lambda_right_cm, T_right_K = 0.1, 270.0, 0.4

# --- Time grids (independent windows, common Δt) ---
t_left_ps  = np.arange(0.0, 2.0 + 0.5*dt_ps, dt_ps)   # 0–2.0 ps
t_right_ps = np.arange(0.0, 1.3 + 0.5*dt_ps, dt_ps)   # 0–1.3 ps

case_left  = run_case(lambda_left_cm,  tau_left_ps,  T_left_K,  t_left_ps)
case_right = run_case(lambda_right_cm, tau_right_ps, T_right_K, t_right_ps)

# ---------- Plot ----------
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), sharex=False)
for j, (ax, case) in enumerate(zip(axes, [case_left, case_right])):
    t_plot = case['t_ps']
    ax.plot(t_plot, case['rho_pm_imag'], label=r'$\mathrm{Im}[\rho_{+-}]$'
, linewidth=1.4, color='red')
    ax.set_ylim(-0.7, 0.7)
    # Pure-dephasing half-life estimate: ln(2)/(2 * gamma_phi)
    if case['gamma_phi_ps'] > 0:
        t_half_ps = np.log(2.0) / (2.0 * case['gamma_phi_ps'])
        ax.axvline(t_half_ps, linestyle='dashed', alpha=0.7, label=rf"$t_{{1/2}} \approx {t_half_ps:.3g}\,\mathrm{{ps}}$", color='red')
    ax.set_xlabel("Time (ps)", fontsize=12)
    if j == 0:
        ax.set_ylabel(r'Coherence', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best', title=rf"$T$ = {1000*case['T_K']:.4g} mK")

plt.tight_layout()
plt.savefig("Fig4_cryo_coherence_1x2_v2.0.0-arxiv.png", dpi=300, bbox_inches="tight")
print("Saved: Fig4_cryo_coherence_1x2_v2.0.0-arxiv.png")
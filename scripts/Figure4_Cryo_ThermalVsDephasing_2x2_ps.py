# cryo_compare_dephasing_therm_2x2_phase_half_ps_ps_indie_site_energy_legend_min.py
# -----------------------------------------------------------------------------------
# 2x2 figure for a two-site excitonic dimer.
# Left column: site populations and |rho_12|
# Right column: energy (|+>,|->) populations ONLY (no coherence trace)
# Coupling J = -34 meV (fixed).
#
# IMPORTANT (per user request):
# - Legends: show ONLY temperature (T) and gamma_phi in the legend title.
# - Curve labels use LaTeX r'$\rho_{11}$' style for ALL populations/coherences.
# - Left column legend contains population and coherence entries (site basis).
# - Right column legend contains strictly population entries (energy basis).
# - Keep two independent time grids (linspace) for top and bottom rows.
# - TOP ROW is now in **picoseconds** (ps). BOTTOM ROW also in ps (independent grid).
# -----------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, ket2dm, mesolve, Qobj

# ---------- Physical constants (in meV & ps) ----------
hbar_meV_ps = 0.6582119514      # ħ in meV·ps
k_B_meV_per_K = 0.08617333262   # k_B in meV/K

# ---------- Fixed system parameter ----------
J_meV = -34.0                    # site-site electronic coupling (meV); energy gap |+>-|-> = 2|J|
sim_unit_ps = 0.659              # 1 simulation unit = 0.659 ps

# ---------- Solver options ----------
opts = {"progress_bar": False, "atol": 1e-10, "rtol": 1e-8, "nsteps": 10000}

# ---------- Conversions ----------
cm_to_meV = lambda x_cm: x_cm / 8.065544  # 1 cm^-1 ≈ 0.123984 meV

# ---------- Numerically-stable Bose–Einstein helper ----------
def bose_einstein(DeltaE_meV, T_K):
    if T_K <= 0:
        return 0.0
    x = DeltaE_meV / (k_B_meV_per_K * T_K)
    if x > 50.0:
        return 0.0
    if x < 1e-6:
        return 1.0 / x
    return 1.0 / np.expm1(x)

# ---------- Dephasing rate model ----------
def gamma_phi(lambda_fast_meV, tau_c_ps, T_K):
    k_c = 1.0 / tau_c_ps  # ps^-1
    return (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2 * k_c)

# ---------- |+> <-> |-> thermal population transfer ----------
def energy_transfer_rates(lambda_fast_meV, tau_c_ps, T_K, DeltaE_meV):
    k_c = 1.0 / tau_c_ps  # ps^-1
    omega0_ps_inv = DeltaE_meV / hbar_meV_ps
    S_cl = (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps**2) * (k_c / (k_c**2 + omega0_ps_inv**2))
    n_BE = bose_einstein(DeltaE_meV, T_K)
    k_down_ps = S_cl * (n_BE + 1.0)  # |+> -> |->
    k_up_ps   = S_cl * n_BE          # |-> -> |+>
    return k_down_ps, k_up_ps

# ---------- Hamiltonian (site basis) ----------
J_sim = J_meV * (sim_unit_ps / hbar_meV_ps)  # dimensionless in sim-units
H_sim = Qobj([[0.0, J_sim], [J_sim, 0.0]])

# ---------- Bases ----------
ket1 = basis(2, 0)
ket2 = basis(2, 1)
ket_plus  = (ket1 + ket2).unit()
ket_minus = (ket1 - ket2).unit()
ket_phi   = (ket1 + 1j*ket2).unit()  # π/2 phase superposition

# ---------- Collapse operators builder ----------
def make_c_ops(lambda_fast_meV, tau_c_ps, T_K):
    # (A) Pure dephasing on sites
    gphi_ps = gamma_phi(lambda_fast_meV, tau_c_ps, T_K)        # ps^-1
    gphi_sim = gphi_ps * sim_unit_ps                           # 1/sim-unit
    c_ops = [np.sqrt(max(gphi_sim, 0.0)) * ket1 * ket1.dag(),
             np.sqrt(max(gphi_sim, 0.0)) * ket2 * ket2.dag()]

    # (B) Thermal |+> <-> |-> population transfer
    DeltaE_meV = 2.0 * abs(J_meV)
    k_down_ps, k_up_ps = energy_transfer_rates(lambda_fast_meV, tau_c_ps, T_K, DeltaE_meV)
    L_down = np.sqrt(max(k_down_ps * sim_unit_ps, 0.0)) * ket_minus * ket_plus.dag()
    L_up   = np.sqrt(max(k_up_ps   * sim_unit_ps, 0.0)) * ket_plus  * ket_minus.dag()
    c_ops += [L_down, L_up]
    return c_ops, gphi_ps

# ---------- Single-case runner ----------
def run_case(lambda_fast_cm, tau_c_ps, T_K, t_ps):
    """
    Return dict with time (ps), site/energy populations and coherences for plotting.
    SITE-BASIS EVOLUTION IS INITIALIZED IN (|1> + i|2>)/sqrt(2).
    ENERGY-BASIS EVOLUTION USES |+>.
    """
    lam_meV = cm_to_meV(lambda_fast_cm)
    c_ops, gphi_ps = make_c_ops(lam_meV, tau_c_ps, T_K)

    # time arrays in sim units and ps for plotting
    t_sim = t_ps / sim_unit_ps

    # Initial states
    psi_site0 = ket_phi
    psi_energy0 = ket_plus

    # Evolve
    res_site = mesolve(H_sim, psi_site0,   t_sim, c_ops=c_ops, e_ops=None, options=opts)
    res_en   = mesolve(H_sim, psi_energy0, t_sim, c_ops=c_ops, e_ops=None, options=opts)

    # Extract populations and coherences
    p11 = []
    p22 = []
    rho12_abs = []
    pPP = []
    pMM = []

    P1 = ket2dm(ket1)
    P2 = ket2dm(ket2)
    Pp = ket2dm(ket_plus)
    Pm = ket2dm(ket_minus)
    op12 = ket1 * ket2.dag()         # |1><2|

    for st_site, st_en in zip(res_site.states, res_en.states):
        # site
        p11.append(np.real((P1 * st_site).tr()))
        p22.append(np.real((P2 * st_site).tr()))
        rho12_abs.append(abs((op12 * st_site).tr()))
        # energy (populations only)
        pPP.append(np.real((Pp * st_en).tr()))
        pMM.append(np.real((Pm * st_en).tr()))

    return {
        't_ps': t_ps,
        'p11': np.array(p11), 'p22': np.array(p22), 'rho12_abs': np.array(rho12_abs),
        'pPP': np.array(pPP), 'pMM': np.array(pMM),
        'gamma_phi_ps': gphi_ps,
        'lambda_fast_cm': lambda_fast_cm, 'tau_c_ps': tau_c_ps, 'T_K': T_K
    }

# ---------- Two parameter sets (unchanged) ----------
# Row 1:
tau1_ps, lambda1_cm, T1_K = 1.0, 20.0, 0.01
# Row 2:
tau2_ps, lambda2_cm, T2_K = 1.0, 20.0, 0.1

# ---------- Separate time grids (independent) ----------
# Top row: (now in ps for plotting as requested)
t_top_ps    = np.linspace(0.0, 200.0, 9000)    # ps #min_sampling:8350
# Bottom row:
t_bottom_ps = np.linspace(0.0, 20.0, 700)      # ps #min_sampling: 668

# ---------- Run the two cases with their own times ----------
case_top = run_case(lambda1_cm, tau1_ps, T1_K, t_top_ps)      # returns time in ps
case_bot = run_case(lambda2_cm, tau2_ps, T2_K, t_bottom_ps)   # returns time in ps

# ---------- Plot: 2x2 grid (both rows in ps) ----------
fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.8), sharex=False)

rowdata = [
    (case_top, "Row 1", "ps"),
    (case_bot, "Row 2", "ps"),
]

for i, (case, labeltxt, unit) in enumerate(rowdata):
    t_plot = case['t_ps']
    xlab   = "Time (ps)"
    title_unit_tag = "t in ps"

    # Legend title: ONLY T and (NOT) gamma_phi
    legend_title = rf"$T$ = {1000*case['T_K']:.4g} mK"#,  $\gamma_\phi$ = {1000*case['gamma_phi_ps']:.3g}$\times 10^{{-3}}$ ps$^{{-1}}$"

    # Left: site basis + |rho_12|
    axL = axes[i, 0]
    axL.plot(t_plot, case['p11'],      label=r'$\rho_{11}$', color='dodgerblue', linestyle='solid', linewidth=1.5)
    axL.plot(t_plot, case['p22'],      label=r'$\rho_{22}$', color='black', linestyle='solid', linewidth=1.2)
    axL.plot(t_plot, case['rho12_abs'], color='white', linestyle='solid', linewidth = 1.5)
    axL.plot(t_plot, case['rho12_abs'],label=r'$|\rho_{12}|$', color='red', linestyle='solid', linewidth=1.2)
    

    # --- coherence half-life vertical line (site basis only) ---
    t_half_ps = (np.log(2.0) / (2.0 * case['gamma_phi_ps'])) if case['gamma_phi_ps'] > 0 else None
    if t_half_ps is not None:
        axL.axvline(
    t_half_ps, color='red', linestyle='dashed', alpha=0.7,
    label=rf"$t_{{1/2}} = {t_half_ps:.3g}\,\mathrm{{ps}}$"
    )

    axL.set_ylabel('Population / Coherence', fontsize=18)
    axL.grid(True, alpha=0.3)
    #if i == 0:
    #    axL.set_title(f'Site basis (init: (|1>+i|2>)/sqrt(2)); J=-34 meV; {title_unit_tag}')
    axL.legend(fontsize=10, loc='best', title=legend_title)
    axL.set_xlabel(xlab, fontsize=18)

    # Right: energy basis (populations only)
    axR = axes[i, 1]
    axR.plot(t_plot, case['pPP'], label=r'$\rho_{++}$', color='dodgerblue', linestyle='dashed')
    axR.plot(t_plot, case['pMM'], label=r'$\rho_{--}$', color='black', linestyle='dashed')
    axR.grid(True, alpha=0.3)
    #if i == 0:
    #    axR.set_title(f'Energy basis: |+>, |-> (init: |+>); {title_unit_tag}')
    axR.legend(fontsize=10, loc='best', title=legend_title)
    axR.set_xlabel(xlab, fontsize=18)

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

    plt.savefig(r'/mnt/data/cryo_compare_dephasing_therm_2x2_phase_half_ps_ps_indie_site_energy_legend_min_t_1_2_v1.1.0-arxiv.png', dpi=300, bbox_inches='tight')
    print('Saved:', r'/mnt/data/cryo_compare_dephasing_therm_2x2_phase_half_ps_ps_indie_site_energy_legend_min_t_1_2_v1.1.0-arxiv.png')
# Optional logging
print(f"Row 1: T={T1_K:.6g} K (TOP row plotted in ps; {t_top_ps[0]}–{t_top_ps[-1]} ps)")
print(f"Row 2: T={T2_K:.6g} K (BOTTOM row plotted in ps; {t_bottom_ps[0]}–{t_bottom_ps[-1]} ps)")

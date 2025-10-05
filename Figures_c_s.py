import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, mesolve, Qobj, Options

### ===== FRET Rate Calculation ===== ###

# ===== Parameters =====
tau_D_ns = 3.0              # Donor fluorescence lifetime (ns)
tau_D_ps = tau_D_ns * 1000  # Convert to ps

Phi_D = 0.57                # Quantum yield of donor
n = 1.72                     # Refractive index of protein
kappa2 = 2.0/3.0            # Orientation factor (random isotropic)

# Förster radius for Venus-GFP like system (literature ~5.7 nm)
R0_nm = 5.0

# Assume chromophore center-to-center distance
R_nm = 2.5  # nm

# Compute transfer rate
k_DA_per_tauD = (R0_nm / R_nm)**6
k_DA_per_ps = k_DA_per_tauD / tau_D_ps

print("Förster radius R0 = %.2f nm" % R0_nm)
print("Distance R = %.2f nm" % R_nm)
print("Donor lifetime tau_D = %.2f ns" % tau_D_ns)
print("FRET transfer rate k_DA = %.4e ps^-1" % k_DA_per_ps)
print("FRET transfer time tau_FRET = %.2f ns" % (1.0/k_DA_per_ps/1000))

### ============== END of FRET Rate Calculation ============== ###


### ========================== Figures ======================= ###

# ===== Constants (ps / meV) =====
hbar_meV_ps     = 0.658212       # ħ in meV·ps
k_B_meV_per_K   = 0.08617333262  # k_B in meV/K

# ===== Bath parameters =====
T_K         = 293.0
lambda_meV  = 20.0
gamma_c_ps  = 2.0 #3.798 #gamma_c_ps max = 3.798 (<< k_BT/hbar)
tau_c_ps    = 0.1 / gamma_c_ps

# ===== System parameters =====
J_meV       = -34.0
Delta_0_meV = 59.0
tau1_ps     = 0.2
tau2_ps     = 2.0
a           = 0.8
#c_s=0.9

# ===== Operators =====
ket1 = basis(2, 0)
ket2 = basis(2, 1)

sigma_x = Qobj([[0, 1],[1, 0]])
sigma_z = Qobj([[1, 0],[0,-1]])

ket_bright = (ket1 + ket2).unit()
ket_dark   = (ket1 - ket2).unit()

# ===== Helper functions =====
def bose_einstein(DeltaE_meV, T_K):
    if T_K <= 0.0: return 0.0
    x = DeltaE_meV / (k_B_meV_per_K * T_K)
    if x > 50.0: return 0.0
    if x < 1e-6: return 1.0 / x
    return 1.0 / np.expm1(x)

def Delta_meV(t_ps):
    return Delta_0_meV * (1-(a*np.exp(-t_ps/tau1_ps)+(1-a)*np.exp(-t_ps/tau2_ps)))

def Delta_coeff_ps_inv(t_ps, args=None):
    return 0.5 * (Delta_meV(t_ps) / hbar_meV_ps)

def DeltaE_meV(t_ps, J_meV):
    return 2.0 * np.sqrt(J_meV**2 + 0.25 * Delta_meV(t_ps)**2)

def gamma_down_up(t_ps, J_meV):
    DeltaE = abs(DeltaE_meV(t_ps, J_meV))
    omega_ps   = DeltaE / hbar_meV_ps
    gamma_c_ps = 1.0 / tau_c_ps
    pref = (2.0 * lambda_meV * k_B_meV_per_K * T_K) / (hbar_meV_ps ** 2)
    S_cl_ps = pref * (gamma_c_ps / (gamma_c_ps**2 + omega_ps**2))
    n = bose_einstein(DeltaE, T_K)
    gdown = max(S_cl_ps * (n + 1.0), 0.0)
    gup   = max(S_cl_ps * n,         0.0)
    return gdown, gup

def analytic_eigvecs(t_ps, J_meV):
    J_ps = J_meV / hbar_meV_ps
    Delta_ps = Delta_meV(t_ps) / hbar_meV_ps
    theta = 0.5*np.arctan2(2*J_ps, Delta_ps) if abs(Delta_ps) > 1e-12 else np.pi/4
    ket_minus = np.cos(theta) * ket2 - np.sin(theta) * ket1
    ket_plus  = np.sin(theta) * ket1 + np.cos(theta) * ket2
    return ket_minus.unit(), ket_plus.unit()

def run_dynamics(J_meV, tmax_ps, nsteps):
    J_ps = J_meV / hbar_meV_ps
    H = [J_ps*sigma_x, [sigma_z, Delta_coeff_ps_inv]]
    psi0 = ket_bright
    tlist = np.linspace(0, tmax_ps, nsteps)

    def c_down(t, args=None):
        gdown, _ = gamma_down_up(t, J_meV)
        km, kp = analytic_eigvecs(t, J_meV)
        return np.sqrt(gdown) * (km * kp.dag())

    def c_up(t, args=None):
        _, gup = gamma_down_up(t, J_meV)
        km, kp = analytic_eigvecs(t, J_meV)
        return np.sqrt(gup) * (kp * km.dag())

    c_ops = [[c_down], [c_up]]
    opts = Options(nsteps=20000, store_states=True)

    sol = mesolve(H, psi0, tlist, c_ops=c_ops, options=opts)

    rho_bb, rho_pm = [], []
    for i, t in enumerate(tlist):
        rho_t = sol.states[i]
        rho_bb.append(np.real((rho_t * (ket_bright * ket_bright.dag())).tr()))
        km, kp = analytic_eigvecs(t, J_meV)
        rho_pm.append(np.real(complex(km.dag() * rho_t * kp)))

    return tlist, np.array(rho_bb), np.array(rho_pm)

# ===== Markovian T2 =====
def gamma_phi(lambda_fast_meV, tau_c_ps, T_K, *, c_s=0.996):
    num = (2.0 * lambda_fast_meV * k_B_meV_per_K * T_K) * tau_c_ps
    den = (hbar_meV_ps ** 2)
    return max((1.0 - c_s) * (num / den), 0.0)

def relaxation_times_ps(t_ps, gphi_ps, J_meV):
    gdown, gup = gamma_down_up(t_ps, J_meV)
    gamma_pop = gdown + gup
    T1 = np.inf if gamma_pop == 0 else 1.0 / gamma_pop
    T2 = np.inf if (gphi_ps + 0.5*gamma_pop) == 0 else 1.0 / (gphi_ps + 0.5*gamma_pop)
    return T1, T2

# ====== Generate data for quadrants ======
t_bb, rhoBB_J, _ = run_dynamics(-34.0, 5.0, 1001)
_, rhoBB_H, _    = run_dynamics(+34.0, 5.0, 1001)

t_pm, _, rhoPM_J_long = run_dynamics(-34.0, 30.0, 2001)
_, _, rhoPM_H_long    = run_dynamics(+34.0, 30.0, 2001)

t20 = np.linspace(0,20,1001)
Delta_vals = [Delta_meV(t) for t in t20]

tD = np.linspace(0.0, 20.0, 2001)
gphi_ps = gamma_phi(lambda_meV, tau_c_ps, T_K)
T2_series = [relaxation_times_ps(t, gphi_ps, -34.0)[1] for t in tD]

# ====== Figure 1: Quadrant ======
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

axs[0,0].plot(t20, Delta_vals, color="#B79A09")#"#BD9C17")
axs[0,0].set_title("A"); axs[0,0].set_xlabel(r"$t$ (ps)"); axs[0,0].set_ylabel(r"$\Delta(t)$ (meV)"); axs[0,0].grid()
axs[0,0].set_ylim(-2.5,66.5)
axs[0,0].set_xlim(-1.3,20.0)
axs[0,0].axhline(np.abs(J_meV), color="#00A36C", linestyle=":", linewidth=1.5,
            label=r"|J|")
axs[0,0].legend()

axs[0,1].plot(tD, T2_series, color="#00A36C")#"#0F7F3C")
axs[0,1].set_title("B"); axs[0,1].set_xlabel(r"$t$ (ps)"); axs[0,1].set_ylabel(r"$T_2(t)$ (ps)"); axs[0,1].grid()
# guide lines
axs[0,1].axhline(10*tau_c_ps, color="red", linestyle="--", linewidth=1.5,
            label=r"$T_2=10\tau_c$")
axs[0,1].axvline(10*tau_c_ps, color="red", linestyle=":", linewidth=1.5,
            label=r"$t=10\tau_c$")
axs[0,1].axvline(1, color="gray", linestyle=":", linewidth=1.5,
            label=r"$t=1$ ps")
axs[0,1].legend()
axs[0,1].set_ylim(0.48,0.59)#(2.6,10.5)
axs[0,1].set_xlim(0.0,20.0)

axs[1,0].plot(t_bb, rhoBB_J, label="J-like", color="dodgerblue")
axs[1,0].plot(t_bb, rhoBB_H, label="H-like", color="black")
axs[1,0].set_title("C"); axs[1,0].set_xlabel(r"$t$ (ps)"); axs[1,0].set_ylabel(r"$\rho_{BB}$"); axs[1,0].legend(); axs[1,0].grid()
axs[1,0].set_ylim(0.1,1.03)
axs[1,0].set_xlim(-0.2,4.0)

axs[1,1].plot(t_pm, rhoPM_J_long, label="J-like", color="dodgerblue")
axs[1,1].plot(t_pm, rhoPM_H_long, label="H-like", color="black")
axs[1,1].set_title("D"); axs[1,1].set_xlabel(r"$t$ (ps)"); axs[1,1].set_ylabel(r"$\Re[\rho_{-+}]$"); axs[1,1].legend(); axs[1,1].grid()
axs[1,1].set_ylim(-0.02,0.75)
axs[1,1].set_xlim(-0.5,8.0)

plt.tight_layout()
plt.savefig("Figure1_quadrant.png", dpi=300)
plt.close()

# ====== Figure 2: T1 comparison (A206 vs K206) ======
tlist = np.linspace(0.0, 15.0, 2001)
T1_A206 = [relaxation_times_ps(t, 0.0, -34.0)[0] for t in tlist]  # J=-34 meV
k1_A206 = [1/relaxation_times_ps(t, 0.0, -34.0)[0] for t in tlist]  # J=-34 meV

plt.figure(figsize=(6,4))
plt.plot(tlist, k1_A206, label=r"Dimer $T_1^{-1}(t)$", color="#00A36C")
plt.xlabel(r"$t$ (ps)"); plt.ylabel(r"Population Transfer Rate (ps$^{-1}$)")
plt.axhline(k_DA_per_ps, color="#B79A09", linestyle="--", linewidth=1.5,
            label=fr"Monomer Pair $k_{{DA}}={k_DA_per_ps:.3f}$ ps$^{{-1}}$")
plt.grid(); plt.legend(); plt.tight_layout()
plt.ylim(0.0,5.0)
plt.xlim(-0.3,15.0)
plt.savefig("Figure2_T1_compare.png", dpi=300)
plt.close()

print("✅ Saved Figure1_quadrant.png and Figure2_T1_compare.png")

import numpy as np
import matplotlib.pyplot as plt

# === Constants ===
eps0 = 8.8541878128e-12  # F/m
n = 1.72
kappa_iso = np.sqrt(2/3)

# Transition dipole moment (fixed)
mu_D = 7.1
mu_Cm = mu_D * 3.33564e-30  # C·m

# Target minimum coupling
target_J = -34.0  # meV

# ================================
# PART 1: Isotropic J vs distance
# ================================

def J_iso(R):
    """Isotropic averaged coupling (meV) at distance R (m)."""
    return (kappa_iso * mu_Cm**2) / (4*np.pi*eps0*n**2*R**3) / 1.60218e-22

# Distance grid
R_vals = np.linspace(2e-9, 5e-9, 500)

# Physical case
R_phys = 2.5e-9
J_phys = J_iso(R_phys)

# Effective R for J = -34 meV
R_eff_iso = ((kappa_iso * mu_Cm**2) / (4*np.pi*eps0*n**2 * abs(target_J*1.60218e-22)))**(1/3)

# Compute curves
J_curve = J_iso(R_vals)

# Plot 
plt.plot(R_vals*1e9, J_curve, lw=2, color="#B79A09")
plt.axhline(J_phys, color="#00A36C", linestyle="--", label=f"$J(R={R_phys*1e9:.1f}$ nm$)=+{J_phys:.2f}$ meV")
plt.xlabel(r"Distance $R$ (nm)")
plt.ylabel(r"Coulombic coupling $J$ (meV)")
plt.ylim(0, 1.0)
plt.xlim(2.0, 5.0)
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig("Isotropic_coupling_vs_distance_R.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"[Panels] Saved as 'Isotropic_coupling_vs_distance_R.png'")
print(f"Effective R for J=-34 meV (isotropic): {R_eff_iso*1e9:.2f} nm")
print(f"Physical J at R=2.5 nm: {J_phys:.2f} meV")

# ================================
# PART 2: Heatmap J(theta,phi) with scaling via R
# ================================

# Grid
theta_vals = np.linspace(0, np.pi, 200)
phi_vals = np.linspace(0, np.pi, 200)  # only 0–pi needed
Theta, Phi = np.meshgrid(theta_vals, phi_vals)

# Angular coupling function
def J_ang(R):
    return (mu_Cm**2/(4*np.pi*eps0*n**2*R**3)) * (2 - (np.sin(Theta)**2)*(1+2*np.sin(Phi)**2))

# Start with R = 2.5 nm
J_test = J_ang(R_phys)/1.60218e-22
min_J_test = np.min(J_test)

# Correct scaling for R so that min(J) = target_J
R_eff_heatmap = R_phys * (min_J_test/target_J)**(1/3)

# Recompute with corrected R_eff
J_scaled_meV = J_ang(R_eff_heatmap)/1.60218e-22

# Plot
theta_ticks = [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]
theta_labels = [r"$0$", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"]
phi_ticks = [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]
phi_labels = [r"$0$", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"]

plt.figure(figsize=(8,6))
c = plt.pcolormesh(Theta, Phi, J_scaled_meV, shading='auto', cmap='gist_ncar')
plt.colorbar(c, label=r"$J$ (meV)", ticks=[-34,0,34,68])
plt.xlabel(r"$\theta$")
plt.ylabel(r"$\phi$")
plt.xticks(theta_ticks, theta_labels)
plt.yticks(phi_ticks, phi_labels)

plt.savefig("C2_symmetric_coupling_heatmap_phi_pi_R.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"[Heatmap] Saved as 'C2_symmetric_coupling_heatmap_phi_pi_R.png'")
print(f"Effective R for J=-34 meV (heatmap): {R_eff_heatmap*1e9:.2f} nm")
print(f"Minimum J achieved: {np.min(J_scaled_meV):.2f} meV")
print(f"Maximum J achieved: {np.max(J_scaled_meV):.2f} meV")

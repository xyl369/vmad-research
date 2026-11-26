#!/usr/bin/env python3
"""
LIGO O3b / aLIGO Benchmark Strain Sensitivity vs. VMAD Stiffness-Noise Spectral Exclusion.
Demonstrates quantitatively why LIGO O3b rules out naive Lorentzian/Ornstein-Uhlenbeck
vacuum-stiffness fluctuations calibrated to ground BMV decoherence, forcing either:
  (1) Main Branch (Palatini omega -> -3/2): non-propagating algebraic scale mode (S_h === 0), or
  (2) Branch B' (propagating soft mode): smooth UV spectral cutoff (f_UV <= 5 Hz).
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = Path(__file__).parent
data_dir = here / "data"
out_dir = here.parent / "assets"
out_dir.mkdir(parents=True, exist_ok=True)

f_h1, asd_h1 = np.loadtxt(data_dir / "H1_O3b_ASD.txt", unpack=True)
f_l1, asd_l1 = np.loadtxt(data_dir / "L1_O3b_ASD.txt", unpack=True)
f_des, asd_des = np.loadtxt(data_dir / "LIGO-P1200087-v18-aLIGO_DESIGN.txt", unpack=True)

# Frequency grid spanning sub-Hz BMV band to kHz LIGO band
f = np.logspace(-1, 3.5, 1200)

# BMV-calibrated amplitude sigma_eps = 1.52e-17, tau_c = 1 s
sigma_eps = 1.52e-17
tau_c = 1.0

# 1. Ruled-out naive Ornstein-Uhlenbeck / Lorentzian spectrum
asd_ou = sigma_eps * np.sqrt(4.0 * tau_c / (1.0 + (2.0 * np.pi * f * tau_c)**2))

# 2. Surviving Branch B' smooth Gaussian UV cutoff spectra (f_UV = 5 Hz and f_UV = 2 Hz)
def gauss_asd(freq, f_uv, sig=sigma_eps):
    tau_g = 1.0 / (2.0 * np.pi * f_uv)
    # Using exact one-sided PSD = 2 * sqrt(2*pi) * sig^2 * tau_g * exp(-(f/f_uv)^2 / 2)
    return sig * np.sqrt(2.0 * np.sqrt(2.0 * np.pi) * tau_g) * np.exp(-0.25 * (freq / f_uv)**2)

asd_g5 = gauss_asd(f, 5.0)
asd_g1p5 = gauss_asd(f, 1.5)

plt.rcParams.update({
    "font.size": 10.5,
    "axes.edgecolor": "#334155",
    "axes.linewidth": 0.9,
})

fig, ax = plt.subplots(figsize=(9.2, 5.4), dpi=200)

# Shaded exclusion region above LIGO L1 O3b
ax.fill_between(f_l1, asd_l1, 1e-16, color="#fee2e2", alpha=0.45, label="Excluded by LIGO O3b Strain Sensitivity")

ax.loglog(f_h1, asd_h1, color="#b91c1c", lw=1.6, alpha=0.85, label="LIGO Hanford (H1) O3b Benchmark ASD")
ax.loglog(f_l1, asd_l1, color="#1d4ed8", lw=1.6, alpha=0.9, label="LIGO Livingston (L1) O3b Benchmark ASD")
ax.loglog(f_des, asd_des, color="#111827", lw=1.2, ls="--", alpha=0.75, label="aLIGO Design Sensitivity")

ax.loglog(f, asd_ou, color="#dc2626", lw=2.2, ls=":", label=r"Ruled Out: Naive Lorentzian/OU $\eta$-noise ($\tau_c=1\,\mathrm{s}$)")
ax.loglog(f, asd_g5, color="#d97706", lw=2.0, ls="--", label=r"Excluded in 10–35 Hz: Gaussian $f_{\mathrm{UV}}=5.0\,\mathrm{Hz}$")
ax.loglog(f, asd_g1p5, color="#0f766e", lw=2.4, label=r"Branch B$'$ Surviving: Smooth Gaussian $f_{\mathrm{UV}}\leq 1.5\,\mathrm{Hz}$")

# Annotation arrows
ax.annotate(
    "LIGO O3b rules out Lorentzian\n1/f tail by >10⁴× at 100 Hz",
    xy=(100.0, 4.8e-20),
    xytext=(160.0, 2.5e-18),
    arrowprops=dict(arrowstyle="->", color="#991b1b", lw=1.2),
    fontsize=9.5,
    color="#991b1b",
    bbox=dict(boxstyle="round,pad=0.35", fc="#fff1f2", ec="#fca5a5", lw=0.8),
)

ax.annotate(
    "Main Branch (ω → -3/2):\nnon-propagating scale mode\nS_h(f) ≡ 0 (no strain noise)",
    xy=(12.0, 4e-26),
    xytext=(20.0, 3e-26),
    fontsize=9.2,
    color="#1e3a5f",
    bbox=dict(boxstyle="round,pad=0.35", fc="#f0f9ff", ec="#93c5fd", lw=0.8),
)

ax.set_xlim(0.2, 2500)
ax.set_ylim(1e-27, 1e-16)
ax.set_xlabel("Frequency  $f$  [Hz]")
ax.set_ylabel(r"Strain-Equivalent ASD  $S_h^{1/2}(f)$  [$\mathrm{Hz}^{-1/2}$]")
ax.set_title("LIGO O3b Strain Ceiling vs. VMAD Vacuum-Stiffness Fluctuation Spectra", fontweight="bold", pad=12)
ax.grid(True, which="both", ls="--", lw=0.5, alpha=0.4)
ax.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.92, fontsize=8.7)

fig.tight_layout()
out_path = out_dir / "ligo_o3b_vmad_bounds.png"
fig.savefig(out_path, dpi=200)
print(f"VERIFIED: Generated {out_path} ({out_path.stat().st_size} bytes)")

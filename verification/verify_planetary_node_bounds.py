#!/usr/bin/env python3
"""
VMAD Academicized Planetary & Solar-System Physical Information Node Verification Suite
Aligned strictly to the 2026-11-21 Unified Arbitration (zero retracted claims).

Verifies all quantitative mathematical-physical bounds across seven Solar-System
environments under combined Information Thermodynamics (Landauer, Stefan-Boltzmann,
Arrhenius), Geophysical/Astrophysical Observability (Diviner, ICESat-2, IGETS, LIGO),
Relativistic Causality, Electromagnetic/Cosmic-Ray Shielding, and Emergent-Spacetime
Dual-Channel Anomalous Decoherence (Branch B' g^2 & R^2 vs Main Branch omega -> -3/2).

Generates Figure 2: assets/solar_system_node_pareto.png
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Fundamental & geophysical constants (SI)
c = 2.99792458e8            # m / s
G0 = 6.67430e-11            # m^3 / (kg s^2)
hbar = 1.054571817e-34      # J s
kB = 1.380649e-23           # J / K
sigma_SB = 5.670374419e-8   # W / (m^2 K^4)
eV = 1.602176634e-19        # J
eps0 = 8.8541878128e-12     # F / m
mu0 = 4.0 * np.pi * 1e-7    # H / m
yr_s = 365.25 * 86400.0
AU = 1.495978707e11         # m

print("=== [1] Thermodynamic & Radiation Ceilings (Landauer / Stefan-Boltzmann / Convective) ===")
j_30K = sigma_SB * (30.0**4)
area_100MW_30K_km2 = (1.0e8 / j_30K) / 1e6
print(f"30 K blackbody flux = {j_30K:.4f} W/m^2  => 100 MW requires {area_100MW_30K_km2:.1f} km^2")
assert 0.045 < j_30K < 0.047
assert area_100MW_30K_km2 > 2000.0

eps_rad = 0.9
j_18p1K = eps_rad * sigma_SB * (18.10**4 - 3.0**4)
area_10mW_18K = 10e-3 / j_18p1K
COP_carnot = 30.0 / (300.0 - 30.0)
COP_real = 0.20 * COP_carnot
P_elec_100MW_lunar = 100.0e6 / COP_real
P_total_rej_lunar = 100.0e6 + P_elec_100MW_lunar
j_300K = eps_rad * sigma_SB * (300.0**4 - 3.0**4)
area_rim_km2 = (P_total_rej_lunar / j_300K) / 1e6
print(f"18.10 K PSR flux (eps=0.9) = {j_18p1K*1e3:.2f} mW/m^2 => 10 mW passive plate = {area_10mW_18K:.2f} m^2")
print(f"Veto of 100 MW lunar PSR active node: needs {P_elec_100MW_lunar/1e9:.2f} GW power & {area_rim_km2:.2f} km^2 rim radiators")
assert 5.3e-3 < j_18p1K < 5.6e-3
assert 1.75 < area_10mW_18K < 1.95
assert abs(P_elec_100MW_lunar / 1e9 - 4.50) < 0.05
assert 11.0 < area_rim_km2 < 11.5

cp_sw = 4200.0
mdot_1MW = 1.0e6 / (cp_sw * 1.0)
print(f"Deep-ocean convective cooling (dT <= 1 K): 1 MW requires {mdot_1MW:.1f} kg/s cold water throughput")
assert 237.0 < mdot_1MW < 239.0

L_ice = 3.34e5
melt_100MW_t_per_yr = (1.0e8 / L_ice) * yr_s / 1e3
print(f"Antarctic 100 MW basal melt = {melt_100MW_t_per_yr/1e6:.2f} million tonnes/yr => ICESat-2 subsidence veto -> stealth cap <= 1-10 kW")
assert 9.3e6 < melt_100MW_t_per_yr < 9.6e6

print("\n=== [2] EM Opacity, Cosmic-Ray Attenuation & Deep-Time Thermal/Arrhenius Physics ===")
sigma_sw = 4.0
tau_e_sw = (80.0 * eps0) / sigma_sw
f_radar = 3.3e9
delta_skin_3p3GHz = np.sqrt(2.0 / (2.0 * np.pi * f_radar * mu0 * sigma_sw))
f_elf = 76.0
delta_skin_elf = np.sqrt(2.0 / (2.0 * np.pi * f_elf * mu0 * sigma_sw))
atten_elf_4km = np.exp(-2.0 * 4000.0 / delta_skin_elf)
print(f"Seawater tau_e = {tau_e_sw:.2e} s; 3.3 GHz skin depth = {delta_skin_3p3GHz*1e3:.2f} mm; 76 Hz ELF 4 km power atten = {atten_elf_4km:.2e}")
assert 1.5e-10 < tau_e_sw < 2.0e-10
assert 4.2e-3 < delta_skin_3p3GHz < 4.6e-3
assert atten_elf_4km < 1e-100

col_50m_g_cm2 = 5000.0 * 3.0
gcr_50m = np.exp(-col_50m_g_cm2 / 160.0)
P_synodic = 29.530589 * 86400.0
alpha_basalt = 7.1e-7
l_skin_lunar = np.sqrt(alpha_basalt * P_synodic / np.pi)
dT_atten_50m = np.exp(-50.0 / l_skin_lunar)
q_cond_50m = 1.8 * 10.0 / 50.0
print(f"Schroter 50 m basalt roof: GCR atten = {gcr_50m:.2e}, thermal skin = {l_skin_lunar:.3f} m (50 m damping = {dT_atten_50m:.2e}), q_max = {q_cond_50m:.2f} W/m^2")
assert 1.5e-41 < gcr_50m < 2.5e-41
assert 0.72 < l_skin_lunar < 0.80
assert abs(q_cond_50m - 0.36) < 1e-9

Ea_defect = 0.35 * eV
arrhenius_ratio = np.exp(-(Ea_defect / kB) * (1.0 / 18.0 - 1.0 / 250.0))
print(f"Arrhenius atomic diffusion ratio (18 K / 250 K, Ea=0.35 eV) = {arrhenius_ratio:.2e} (~10^-91 freezing!)")
assert 1e-95 < arrhenius_ratio < 1e-88

A_diviner = 240.0 * 240.0
T_hotspot = (40.0**4 + 1.0e4 / (sigma_SB * A_diviner))**0.25
dT_hotspot = T_hotspot - 40.0
print(f"Lunar PSR Diviner 240 m pixel (T0=40 K) with 10 kW thermal leak: T = {T_hotspot:.2f} K (dT = +{dT_hotspot:.2f} K)")
assert 48.0 < T_hotspot < 49.5

print("\n=== [3] Relativistic Light-Time Causality Stratification ===")
D_EM = 3.844e8
R_moon = 1.7374e6
t_EM_oneway = D_EM / c
t_EM_rtt = 2.0 * t_EM_oneway
arc_ar_sp = R_moon * np.deg2rad(113.7)
t_ar_sp = arc_ar_sp / c
theta_earth_deg = np.rad2deg(2.0 * np.arcsin(6.371e6 / D_EM))
t_jup_min_lo = (4.2 * AU / c) / 60.0
t_jup_min_hi = (6.2 * AU / c) / 60.0
print(f"Earth-Moon one-way = {t_EM_oneway:.3f} s (RTT = {t_EM_rtt:.3f} s); Earth angular dia from Moon = {theta_earth_deg:.4f} deg")
print(f"Aristarchus <-> Lunar South Pole arc = {arc_ar_sp/1e3:.1f} km ({t_ar_sp*1e3:.2f} ms)")
print(f"Earth-Jupiter (Europa/Callisto) one-way = {t_jup_min_lo:.1f}-{t_jup_min_hi:.1f} min")
assert abs(t_EM_oneway - 1.282) < 0.002
assert abs(theta_earth_deg - 1.8993) < 0.002
assert 3440.0 < arc_ar_sp / 1e3 < 3460.0

print("\n=== [4] Cross-Body Anomalous Decoherence Spatial Gradient (Branch B' g^2 & R^2 Regimes) ===")
g_earth, R_earth = 9.80665, 6.371e6
bodies = [
    ("Earth Surface (Abyssal / Antarctica)", 9.80665, 6.371e6),
    ("Moon (PSR / Aristarchus)",             1.625,   1.7374e6),
    ("Europa (Subsurface Ocean)",            1.315,   1.5608e6),
    ("Enceladus (Low-g Ice Shell)",          0.113,   0.2521e6),
]
for name, g_b, R_b in bodies:
    rg2 = (g_b / g_earth)**2
    rR2 = (R_b / R_earth)**2
    print(f"  {name:38s}: g^2 ratio = {rg2:.4e},  R^2 ratio = {rR2:.4e}")

assert abs((1.625 / g_earth)**2 - 0.0275) < 0.001
assert abs((1.7374e6 / R_earth)**2 - 0.0744) < 0.001
assert abs((1.315 / g_earth)**2 - 0.0180) < 0.001
assert abs((1.5608e6 / R_earth)**2 - 0.0600) < 0.001
assert abs((0.113 / g_earth)**2 - 1.33e-4) < 1e-5
assert abs((0.2521e6 / R_earth)**2 - 1.57e-3) < 5e-5

print("\n=== [5] Generating Publication-Grade Solar-System Node Pareto Phase-Space Figure ===")
out_dir = Path(__file__).parent.parent / "assets"
out_dir.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(10.6, 6.4), dpi=200)

sites = [
    # (label, P_stealth_W, coh_gain_R2, coh_gain_g2, color, marker, tx, ty)
    ("1. Lunar South Pole PSR (18.1 K)\nPassive Deep-Time Seed Vault\n[Arrhenius 10⁻⁹¹; Diviner <10 kW]",
     1e-2, 1.0 / 0.0744, 1.0 / 0.0275, "#1d4ed8", "D", 2.2e-3, 1.8e5),
    ("2. East Antarctic Craton (Subglacial)\nLow-Seismic Passive Archive\n[ICESat-2 melt cap ≤1–10 kW]",
     5e3, 1.0, 1.0, "#0284c7", "s", 1.2e-1, 1.6e2),
    ("3. Aristarchus Lava Tube (250 K)\nThermal-Buffered Periodic Forge\n[50 m basalt GCR 10⁻⁴¹; q≤0.36 W/m²]",
     1e3, 1.0 / 0.0744, 1.0 / 0.0275, "#7c3aed", "^", 4.5e0, 2.2e6),
    ("4. Earth Abyssal Basalt / Hydrothermal\nLow-Latency High-Power Heat Sink\n[238 kg/s/MW; EM atten 10⁻¹²¹]",
     1e8, 1.0, 1.0, "#b91c1c", "o", 2.5e5, 6.0e1),
    ("5. Europa Subsurface Ocean\nAutonomous Heavy-Archive Ark\n[10¹¹–10¹² W tidal; 35–52 min lag]",
     1e10, 1.0 / 0.0600, 1.0 / 0.0180, "#0f766e", "P", 2.5e7, 2.5e5),
    ("6. Enceladus / Small Icy Bodies\nHigh-Coherence Low-g Node\n[g² gain ~7,530×; small heat sink]",
     1e5, 1.0 / 1.57e-3, 1.0 / 1.33e-4, "#059669", "v", 1.2e4, 1.5e8),
    ("7. L4/L5 & DRO Free-Fall Orbits\nCompact Relic & Zero-AD Limit\n[10¹²–10¹⁷ g window; no surface support]",
     1e2, 1e14, 1e22, "#d97706", "*", 2.2e-3, 1.5e13),
]

for label, p_w, g_r2, g_g2, col, mk, tx, ty in sites:
    ms = 14 if mk == "*" else 9
    ax.loglog(p_w, g_r2, marker=mk, markersize=ms, color=col, markeredgecolor="#0f172a", markeredgewidth=0.8, zorder=5)
    if g_g2 > g_r2 * 1.1:
        y_top = min(g_g2, 3e15)
        ax.plot([p_w, p_w], [g_r2, y_top], color=col, lw=1.8, alpha=0.55, ls="-", zorder=4)
    ax.annotate(
        label,
        xy=(p_w, g_r2),
        xytext=(tx, ty),
        fontsize=8.1,
        color="#0f172a",
        arrowprops=dict(arrowstyle="->", color=col, lw=1.0, alpha=0.85, connectionstyle="arc3,rad=0.06"),
        bbox=dict(boxstyle="round,pad=0.32", fc="#ffffff", ec=col, lw=0.85, alpha=0.95),
        zorder=6,
    )

# Shade thermodynamically forbidden deep-space corner
ax.axvspan(1e6, 1e12, ymin=0.68, ymax=1.0, color="#fee2e2", alpha=0.35, zorder=1)
ax.text(
    3e8, 3e13,
    "Thermodynamically Forbidden in Deep Space\n(Stefan–Boltzmann 30 K flux = 0.046 W/m²\n→ 100 MW requires >2,100 km² radiators)",
    fontsize=8.4, color="#991b1b", ha="center", va="center",
    bbox=dict(boxstyle="round,pad=0.35", fc="#fff1f2", ec="#fca5a5", lw=0.9),
)

# Pareto frontier connecting optimal complementary niches (L4/L5 -> Enceladus -> Europa)
pf_x = [1e2, 1e5, 1e10]
pf_y = [1e14, 1.0 / 1.57e-3, 1.0 / 0.0600]
ax.loglog(pf_x, pf_y, color="#334155", ls="--", lw=1.4, alpha=0.75, label="Solar-System Pareto Trade-Off Frontier (Coherence vs. Stealth Power)")

ax.set_xlim(1e-3, 1e12)
ax.set_ylim(0.2, 2e17)
ax.set_xlabel("Stealth Steady-State Thermal Dissipation Ceiling  $P_{\\max}$  [W]", fontsize=10.5)
ax.set_ylabel(r"Branch B$'$ Decoherence Suppression  $(\Gamma_{\mathrm{AD}} / \Gamma_{\oplus})^{-1}$", fontsize=10.5)
ax.set_title("VMAD Solar-System Information Nodes: Thermodynamic vs. Decoherence Pareto Frontier", fontsize=11.5, fontweight="bold", pad=12)
ax.grid(True, which="both", ls="--", lw=0.5, alpha=0.38)
ax.legend(loc="upper center", bbox_to_anchor=(0.35, 0.99), frameon=True, facecolor="white", framealpha=0.95, fontsize=8.4)

fig.tight_layout()
out_path = out_dir / "solar_system_node_pareto.png"
fig.savefig(out_path, dpi=200)
print(f"VERIFIED: Generated {out_path} ({out_path.stat().st_size} bytes)")
print("\nALL PLANETARY NODE THERMODYNAMIC, SHIELDING, CAUSAL & DECOHERENCE ASSERTIONS PASSED.")

#!/usr/bin/env python3
"""
VMAD Planetary & Solar-System Information-Node Numerical Consistency Suite (Final Audited 2026-11-26)
====================================================================================================
Executable numerical consistency checks for the 7 Solar-System physical information & metrology
environments audited in the VMAD submission dossier, generating the 4D-encoded 2D trade-off
projection figure (assets/solar_system_node_pareto.png).

Covers:
  1. Stefan-Boltzmann cryogenic (30 K & 18 K) vs room-temp (300 K) radiator walls, convective seawater
     heat rejection, Antarctic zero-melt conductive cap & ICESat-2 melt veto, and LRO Diviner pixel
     observability ceilings (T0 = 20 K & 40 K).
  2. Shielding & channel asymmetries: seawater EM skin depth / ELF / PLOB vs Earth-Moon vacuum PLOB;
     distinguishing primary GCR hadronic cascade extinction from secondary muon & KREEP dose floors;
     Arrhenius solid-state diffusion freezing (18 K vs 250 K).
  3. Relativistic light-time causality stratification & Earth-Moon geometry.
  4. Per-body Branch B' tidal crossover scales (1.141 - 2.112 AU), causal xi <= 0.002 AU R^2 regime,
     and cross-body parameter-free AD ratios (g^2 & R^2).
  5. Nuclear-density compact relic radiation-pressure acceleration (Jacobs-Starkman-Lynn 2015 window).
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Fundamental constants (SI)
c = 2.99792458e8
G = 6.67430e-11
k_B = 1.380649e-23
eV = 1.602176634e-19
sigma_SB = 5.670374419e-8
eps0 = 8.8541878128e-12
mu0 = 4.0 * np.pi * 1e-7
AU = 1.495978707e11
yr = 3.15576e7

print("=== [1] Thermodynamic, Conduction & Radiometric Observability Ceilings ===")
# 1a. Stefan-Boltzmann cryogenic (30 K) vs room-temp (300 K) radiator wall
j_30K = sigma_SB * 30.0**4
A_100MW_30K_km2 = (1e8 / j_30K) / 1e6
A_100MW_300K_km2 = (1e8 / (sigma_SB * 300.0**4)) / 1e6
j_18K = 0.9 * sigma_SB * (18.10**4 - 3.0**4)  # eps=0.9, T_bg=3 K
A_10mW_18K = 0.01 / j_18K
print(f"30 K blackbody flux = {j_30K:.5f} W/m^2 -> 100 MW requires {A_100MW_30K_km2:.1f} km^2 (vs {A_100MW_300K_km2:.4f} km^2 at 300 K)")
print(f"18.1 K passive plate (eps=0.9) flux = {j_18K*1e3:.2f} mW/m^2 -> 10 mW requires {A_10mW_18K:.2f} m^2")
assert abs(j_30K - 0.04593) < 1e-4
assert abs(A_100MW_30K_km2 - 2177.2) < 2.0
assert abs(A_100MW_300K_km2 - 0.2177) < 1e-3
assert abs(j_18K * 1e3 - 5.47) < 0.05

# 1b. Lunar PSR active heat-pump veto (100 MW from 30 K -> 300 K at 20% Carnot COP)
COP_carnot = 30.0 / (300.0 - 30.0)
COP_real = 0.20 * COP_carnot
P_elec_GW = (100e6 / COP_real) / 1e9
Q_reject_W = 100e6 + P_elec_GW * 1e9
A_rim_km2 = (Q_reject_W / (0.9 * sigma_SB * 300.0**4)) / 1e6
print(f"Lunar PSR 100 MW heat pump (30->300 K, 20% Carnot): P_elec = {P_elec_GW:.2f} GW, 300 K rim radiator = {A_rim_km2:.2f} km^2")
assert abs(P_elec_GW - 4.50) < 0.01
assert abs(A_rim_km2 - 11.13) < 0.1

# 1c. LRO Diviner 240 m x 240 m pixel observability thresholds (T0 = 20 K & 40 K)
A_div = 240.0 * 240.0
P_1K_20K = A_div * sigma_SB * (21.0**4 - 20.0**4)
P_1K_40K = A_div * sigma_SB * (41.0**4 - 40.0**4)
T_div_10kW_40K = (40.0**4 + 1e4 / (A_div * sigma_SB))**0.25
print(f"LRO Diviner 240 m pixel dT=+1 K ceiling: {P_1K_20K:.1f} W (at 20 K), {P_1K_40K:.1f} W (at 40 K)")
print(f"LRO Diviner 240 m pixel (40 K) + 10 kW steady leak -> {T_div_10kW_40K:.2f} K (+{T_div_10kW_40K-40.0:.2f} K)")
assert abs(P_1K_20K - 112.6) < 1.0
assert abs(P_1K_40K - 868.0) < 2.0
assert abs(T_div_10kW_40K - 48.69) < 0.05

# 1d. Seawater convective rejection (Delta T <= 1 K) & Mariana pressure
c_p_water = 4200.0
mdot_1MW = 1e6 / (c_p_water * 1.0)
P_mariana_MPa = 1025.0 * 9.80665 * 10900.0 / 1e6
print(f"Seawater convective flow (1 MW, dT=1 K) = {mdot_1MW:.1f} kg/s; 10.9 km pressure = {P_mariana_MPa:.1f} MPa")
assert abs(mdot_1MW - 238.1) < 0.2
assert abs(P_mariana_MPa - 109.6) < 0.5

# 1e. East Antarctic cold-based zero-melt conductive cap & ICESat-2 basal melt veto
kappa_eff = 0.5 * (2.3 + 2.8)  # ice (2.3) + cratonic granite (2.8) W/(m K)
dT_bed = 10.0                  # Gamburtsev cold-based margin below pressure melting point
P_nomelt_10m_kW = 4.0 * np.pi * kappa_eff * 10.0 * dT_bed / 1e3
P_nomelt_100m_kW = 4.0 * np.pi * kappa_eff * 100.0 * dT_bed / 1e3
L_f_ice = 3.34e5
rho_ice = 917.0
melt_100MW_Mt_yr = (1e8 * yr / L_f_ice) / 1e9
V_100MW_m3_yr = (1e8 * yr / L_f_ice) / rho_ice
dh_3km_cm_yr = (V_100MW_m3_yr / (np.pi * 3000.0**2)) * 100.0
dh_12km_cm_yr = (V_100MW_m3_yr / (np.pi * 12000.0**2)) * 100.0
print(f"Antarctic cold-based zero-melt cap (dT=10 K): {P_nomelt_10m_kW:.2f} kW (R=10 m) to {P_nomelt_100m_kW:.2f} kW (R=100 m)")
print(f"Antarctic 100 MW basal melt = {melt_100MW_Mt_yr:.2f} Mt/yr -> subsidence {dh_12km_cm_yr:.2f}-{dh_3km_cm_yr:.1f} cm/yr (R=12-3 km)")
assert abs(P_nomelt_10m_kW - 3.20) < 0.05
assert abs(P_nomelt_100m_kW - 32.04) < 0.1
assert abs(melt_100MW_Mt_yr - 9.45) < 0.05
assert 2.0 < dh_12km_cm_yr < 3.0 and 35.0 < dh_3km_cm_yr < 38.0

print("\n=== [2] Shielding, Arrhenius Retention & Quantum Channel Asymmetries ===")
# 2a. Seawater EM isolation vs Earth-Moon vacuum PLOB
sigma_sw, eps_r_sw = 4.0, 80.0
tau_e = eps0 * eps_r_sw / sigma_sw
delta_33GHz = 1.0 / np.sqrt(np.pi * 3.3e9 * mu0 * sigma_sw)
delta_76Hz = 1.0 / np.sqrt(np.pi * 76.0 * mu0 * sigma_sw)
atten_ELF_4km = np.exp(-2.0 * 4000.0 / delta_76Hz)
eta_ocean_4km = np.exp(-0.04 * 4000.0)                    # k_abs = 0.04 /m (480 nm)
plob_ocean_4km = eta_ocean_4km / np.log(2.0)              # exact -log2(1 - eta) for eta << 1
eta_vac_EM = (np.pi * 0.5**2 / (1550e-9 * 3.844e8))**2    # 1 m circular aperture diffraction limit
plob_vac_EM = -np.log1p(-eta_vac_EM) / np.log(2.0)
print(f"Seawater tau_e = {tau_e:.2e} s; 3.3 GHz approx skin = {delta_33GHz*1e3:.2f} mm; 76 Hz 4 km atten = {atten_ELF_4km:.2e}")
print(f"PLOB 480 nm 4 km ocean = {plob_ocean_4km:.2e} vs 1550 nm Earth-Moon vacuum = {plob_vac_EM:.2e} ebits/mode")
assert abs(tau_e - 1.77e-10) < 1e-11
assert 1e-122 < atten_ELF_4km < 1e-120
assert 1e-70 < plob_ocean_4km < 1e-69
assert 2.0e-6 < plob_vac_EM < 3.0e-6

# 2b. Lunar lava-tube 50 m basalt: hadronic cascade extinction vs thermal skin depth (note: total dose floor is muons + KREEP)
X_50m = 50.0 * 3000.0 * 0.1  # 15,000 g/cm^2
hadronic_50m = np.exp(-X_50m / 160.0)
alpha_th = 7.1e-7  # lunar basalt thermal diffusivity k/(rho*c_p) ~ 1.8/(3000*840) m^2/s
P_synodic = 29.53 * 86400.0
l_skin = np.sqrt(alpha_th * P_synodic / np.pi)
thermal_damp_50m = np.exp(-50.0 / l_skin)
arrhenius_18_vs_250 = np.exp(-(0.35 * eV / k_B) * (1.0 / 18.0 - 1.0 / 250.0))
print(f"50 m basalt (15,000 g/cm^2): primary GCR hadronic exp(-X/160) = {hadronic_50m:.2e} [total dose floor set by secondary muons & KREEP]")
print(f"50 m basalt thermal skin = {l_skin:.3f} m -> diurnal damping = {thermal_damp_50m:.2e}")
print(f"Arrhenius diffusion ratio (18 K vs 250 K, Ea=0.35 eV) = {arrhenius_18_vs_250:.2e}")
assert 1e-41 < hadronic_50m < 3e-41
assert 0.74 < l_skin < 0.78
assert 1e-91 < arrhenius_18_vs_250 < 2e-91

print("\n=== [3] Relativistic Light-Time Causality Stratification ===")
D_EM = 3.844e8
R_moon = 1.7374e6
t_EM_oneway = D_EM / c
t_EM_rtt = 2.0 * t_EM_oneway
arc_ar_sp = R_moon * np.deg2rad(113.7)
t_ar_sp = arc_ar_sp / c
theta_earth_deg = np.rad2deg(2.0 * np.arcsin(6.371e6 / D_EM))
t_jup_lo = (4.2 * AU / c) / 60.0
t_jup_hi = (6.2 * AU / c) / 60.0
print(f"Earth-Moon one-way = {t_EM_oneway:.3f} s (RTT = {t_EM_rtt:.3f} s); Earth angular dia = {theta_earth_deg:.4f} deg")
print(f"Aristarchus <-> Lunar South Pole arc = {arc_ar_sp/1e3:.1f} km ({t_ar_sp*1e3:.2f} ms); Earth-Jupiter = {t_jup_lo:.1f}-{t_jup_hi:.1f} min")
assert abs(t_EM_oneway - 1.282) < 0.002
assert abs(theta_earth_deg - 1.8993) < 0.002

print("\n=== [4] Per-Body Branch B' Tidal Crossovers & Parameter-Free AD Ratios ===")
g_earth, R_earth = 9.80665, 6.371e6
bodies = [
    ("Earth Surface", 5.9722e24, 9.80665, 6.371e6),
    ("Moon",          7.3420e22, 1.625,   1.7374e6),
    ("Europa",        4.7998e22, 1.315,   1.5608e6),
    ("Enceladus",     1.0802e20, 0.113,   0.2521e6),
]
xi_causal_1s = c * 1.0 / AU
print(f"Causal correlation limit at tau_c=1 s: c*tau_c = {xi_causal_1s:.5f} AU (<< all body crossovers -> universal R^2 regime)")
for name, M_b, g_b, R_b in bodies:
    rs_b = 2.0 * G * M_b / c**2
    xi_b = R_b * np.sqrt(R_b / rs_b) / AU
    print(f"  {name:14s}: xi_crit = {xi_b:.4f} AU | g^2 ratio (xi>2.11 AU) = {(g_b/g_earth)**2:.4e} | R^2 ratio (xi<1.14 AU) = {(R_b/R_earth)**2:.4e}")
assert xi_causal_1s < 0.003
assert abs((1.625 / g_earth)**2 - 0.0275) < 0.001
assert abs((1.7374e6 / R_earth)**2 - 0.0744) < 0.001
assert abs((1.315 / g_earth)**2 - 0.0180) < 0.001
assert abs((1.5608e6 / R_earth)**2 - 0.0600) < 0.001
assert abs((0.113 / g_earth)**2 - 1.33e-4) < 1e-5
assert abs((0.2521e6 / R_earth)**2 - 1.57e-3) < 5e-5

print("\n=== [5] Compact Macro Relic Non-Gravitational Acceleration (Jacobs-Starkman-Lynn Window) ===")
M_macro = 1e14  # 10^17 g (upper edge of 55 g - 10^17 g nuclear-density macro window)
rho_nuc = 1e17  # kg/m^3
R_macro = (3.0 * M_macro / (4.0 * np.pi * rho_nuc))**(1.0 / 3.0)
a_rad_1AU = (1361.0 * np.pi * R_macro**2 / c) / M_macro
print(f"10^17 g nuclear-density compact relic: R = {R_macro*100:.2f} cm, solar radiation a_rad(1 AU) = {a_rad_1AU:.3e} m/s^2")
assert a_rad_1AU < 1e-21

print("\n=== [6] Generating Audited 4D-Encoded 2D Trade-Off Projection Figure ===")
out_dir = Path(__file__).parent.parent / "assets"
out_dir.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(10.8, 6.6), dpi=200)

# Nodes plotted on unified ambient-background / observability thermal ceiling P_max [W]
# Horizontal range bars indicate design/no-melt operating point <-> observational detection ceiling
sites = [
    # (label, P_ref_W, P_lo_W, P_hi_W, coh_R2, coh_g2, color, marker, tx, ty)
    ("1. Lunar South Pole PSR (18–20 K)\n3rd/4th-Axis Optimum: Cryo Seed & Metrology\n[Arrhenius 10⁹¹×; 10 mW plate ↔ 0.1–0.9 kW Diviner]",
     5e2, 1e-2, 8.7e2, 1.0 / 0.0744, 1.0 / 0.0275, "#1d4ed8", "D", 2.2e-3, 3.8e5),
    ("2. East Antarctic Craton (Cold-Based)\nLow-Seismic Passive Cold Archive\n[3–32 kW zero-melt ↔ 1–10 MW ICESat-2 veto]",
     1.5e4, 3.2e3, 3.0e6, 1.0, 1.0, "#0284c7", "s", 2.2e-3, 1.6e2),
    ("3. Aristarchus Lava Tube (250 K)\nEarth-Facing Baseline & Periodic Forge\n[50 m hadronic shield; 2.56 s RTT; q≤0.36 W/m²]",
     1e3, None, None, 1.0 / 0.0744, 1.0 / 0.0275, "#7c3aed", "^", 2.5e0, 5.5e6),
    ("4. Earth Abyssal Basalt / Hydrothermal\n<1 ms High-Power Convective Heat Sink\n[238 kg/s/MW; 100–1,000 MW vents; EM 10⁻¹²¹]",
     1e8, 1e8, 1e9, 1.0, 1.0, "#b91c1c", "o", 2.0e5, 8.0e1),
    ("5. Europa Subsurface Ocean\n2D (Power × AD) Corner · Thermally Unfalsifiable\n[10¹¹–10¹² W tidal; 35–52 min lag; ~10⁷ yr ice lag]",
     1e10, None, None, 1.0 / 0.0600, 1.0 / 0.0180, "#0f766e", "P", 1.8e7, 2.8e5),
    ("6. Enceladus / Small Low-g Icy Bodies\nHigh-Coherence Low-g Pareto Bridge\n[g² gain ~7,530× (ξ>2.11 AU); R² ~639×]",
     1e5, None, None, 1.0 / 1.57e-3, 1.0 / 1.33e-4, "#059669", "v", 8.0e3, 1.8e8),
    ("7. Sun–Earth Trojan / Lunar DRO Free-Fall\nZero-AD Limit & Macro-Relic Necessary Habitat\n[55 g–10¹⁷ g macro window; Gyr stability open]",
     1e2, None, None, 1e14, 1e22, "#d97706", "*", 1.8e-3, 1.6e13),
]

for label, p_w, p_lo, p_hi, g_r2, g_g2, col, mk, tx, ty in sites:
    ms = 14 if mk == "*" else 9
    ax.loglog(p_w, g_r2, marker=mk, markersize=ms, color=col, markeredgecolor="#0f172a", markeredgewidth=0.85, zorder=6)
    # Horizontal span bar (design/zero-melt point <-> observability ceiling)
    if p_lo is not None and p_hi is not None:
        ax.plot([p_lo, p_hi], [g_r2, g_r2], color=col, lw=2.2, alpha=0.55, ls="-", zorder=4)
        ax.plot([p_lo, p_hi], [g_r2, g_r2], marker="|", markersize=7, color=col, ls="none", alpha=0.85, zorder=5)
    # Vertical span bar (R^2 -> g^2 AD suppression range)
    if g_g2 > g_r2 * 1.1:
        y_top = min(g_g2, 3e15)
        ax.plot([p_w, p_w], [g_r2, y_top], color=col, lw=1.8, alpha=0.55, ls="-", zorder=4)
    rad_val = -0.14 if mk == "D" else 0.06
    ax.annotate(
        label,
        xy=(p_w, g_r2),
        xytext=(tx, ty),
        fontsize=7.9,
        color="#0f172a",
        arrowprops=dict(arrowstyle="->", color=col, lw=1.0, alpha=0.85, connectionstyle=f"arc3,rad={rad_val}"),
        bbox=dict(boxstyle="round,pad=0.32", fc="#ffffff", ec=col, lw=0.85, alpha=0.95),
        zorder=7,
    )

# Shade cryogenic radiative forbidden corner (explicitly qualified for <=30 K operation!)
ax.axvspan(1e6, 1e12, ymin=0.67, ymax=1.0, color="#fee2e2", alpha=0.35, zorder=1)
ax.text(
    3.0e8, 3.0e13,
    "Cryogenic Radiative Wall (For ≤30 K Operation)\nStefan–Boltzmann 30 K = 0.046 W/m² → 100 MW\nneeds 2,177 km² (300 K still needs 0.22 km² exposed)",
    fontsize=8.1, color="#991b1b", ha="center", va="center",
    bbox=dict(boxstyle="round,pad=0.35", fc="#fff1f2", ec="#fca5a5", lw=0.9),
)

# 2D (Power vs AD) Pareto frontier
pf_x = [1e2, 1e5, 1e10]
pf_y = [1e14, 1.0 / 1.57e-3, 1.0 / 0.0600]
ax.loglog(pf_x, pf_y, color="#334155", ls="--", lw=1.4, alpha=0.75,
          label="2D (Power × AD) Pareto Frontier  [Horizontal bars: design/no-melt ↔ detection ceiling]")

ax.set_xlim(1e-3, 1e12)
ax.set_ylim(0.2, 2e17)
ax.set_xlabel("Background-Consistent Steady-State Thermal Ceiling  $P_{\\max}$  [W]", fontsize=10.3)
ax.set_ylabel(r"Branch B$'$ Decoherence Suppression  $(\Gamma_{\mathrm{AD}} / \Gamma_{\oplus})^{-1}$", fontsize=10.3)
ax.set_title("VMAD Solar-System Information Nodes: 2D Projection of the 4D Trade-Off Space\n(3rd Axis: 18 K Arrhenius Retention 10⁹¹× · 4th Axis: Earth–Node Light-Time & Observability)", fontsize=10.8, fontweight="bold", pad=10)
ax.grid(True, which="both", ls="--", lw=0.5, alpha=0.38)
ax.legend(loc="upper center", bbox_to_anchor=(0.39, 0.99), frameon=True, facecolor="white", framealpha=0.95, fontsize=8.1)

fig.tight_layout()
out_path = out_dir / "solar_system_node_pareto.png"
fig.savefig(out_path, dpi=200)
print(f"VERIFIED: Generated assets/{out_path.name} ({out_path.stat().st_size} bytes)")
print("\nALL AUDITED PLANETARY NODE CONSISTENCY ASSERTIONS PASSED.")

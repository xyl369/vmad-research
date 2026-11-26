# VMAD Research Portal & Reproducible Consistency Suite

**Author:** Xiaoyi Li (Unaffiliated Independent Theoretical Physics Researcher · Non-Commercial)  
**ORCID:** [0009-0000-0482-2376](https://orcid.org/0009-0000-0482-2376)  
**Live Portal:** [xyl369.github.io/vmad-research](https://xyl369.github.io/vmad-research/)  
**Audited Release:** `v2026.11.26` (26 November 2026 · Primary-Literature Audited Edition)

---

## Overview

This repository hosts the public research portal, submission-ready 4-page research statement, curated 1,247-record foundational physics bibliography, chronological self-falsification audit trail, and five executable Python/SymPy/SciPy numerical and symbolic consistency suites for the **VMAD (Vacuum Metastability & Anomalous Decoherence)** foundational physics programme.

The programme investigates two coupled foundational problems:

1. **Emergent-Spacetime Controllability, Vacuum Metastability & Unitary Anomalous Decoherence:**
   - Rigidity No-Go Triad (UV firewall heuristic $\sim 7.49 \times 10^{80}\ \mathrm{J}$; Dirichlet gradient wall $\delta^2 \le \frac{1}{2\omega_{\min}+3}\frac{r_s}{R}$ falsifying canonical soft-shell metric reduction at $\ge 7.22 \times 10^{47}\ \mathrm{J}$; surface weight bounds $\le 8.48\ \mathrm{keV}\text{–}0.431\ \mathrm{MeV}$).
   - **Shape-Stiff, Scale-Soft bifurcation:** Main Branch (Palatini $\omega \to -3/2$ non-propagating scale mode requiring steep normal-phase pinning $\partial^2 U_E/\partial\delta^2 \gg |T_{\mathrm{SM}}|$ to evade Flanagan / Barausse–Sotiriou–Miller singularities, with 1-loop Standard Model radiative survival as the Months 1–2 decisive target) vs. Branch B′ ($2\omega+3>0$ propagating dual-field channel).
   - Endothermic dark-charged cavities ($E''(R_*) > 0$), four-phase Morphology Diagram, and exact $I_h$ icosahedral vanishing ($l=1\dots 5$, $S_6 = 8.09607$) eliminating leading $l=2$ dynamical quadrupolar gravitational radiation.
   - Dual-channel anomalous decoherence: Main Branch zero ambient AD (optical lattice clock probe $\Delta\delta = 3\times 10^{-16} \to \sim 1\ \mathrm{rad/s}$); Branch B′ per-body tidal crossovers (`1.141–2.112 AU`), `60–120×` Hahn spin-echo suppression, `2.53×10⁻¹²` oscillator heating suppression, IceCube *Nat. Phys.* 2024 & KM3NeT clearance, and conditional full-band LIGO O3b exclusion tightening smooth Gaussian cutoffs to **$f_{\mathrm{UV}} \le 1.5\ \mathrm{Hz}$**.

2. **Solar-System Physical Information-Node 4D Trade-Off & Observability Analysis:**
   - Anchored to Landauer (1961), Bekenstein (1981), Lloyd (*Nature* 2000), Sandberg–Armstrong–Ćirković (2017), Dvali & Osmanov (2023), and Wright et al. (*ApJ* 2014).
   - Quantifies steady-state thermal capacity, Branch B′ decoherence suppression, $18\ \mathrm{K}$ Arrhenius retention ($1.15 \times 10^{-91}$), and geophysical/radiometric observability ceilings (LRO Diviner `113–868 W` 1 K pixel ceiling, East Antarctic `3.2–32 kW` zero-melt / `1–10 MW` ICESat-2 veto, `55 g – 10¹⁷ g` non-evaporating compact macro DM window) across seven Solar-System environments, **without claiming detected anomalous artifacts**.

---

## Repository Contents

| Path | Description |
| :--- | :--- |
| [`VMAD_Research_Statement_2026-11-26.pdf`](VMAD_Research_Statement_2026-11-26.pdf) | Submission-ready 4-page research statement & 6-month computational plan |
| [`index.html`](index.html) / [`research-statement.html`](research-statement.html) | Bilingual (EN/ZH) overview portal & research statement |
| [`unified-synthesis.html`](unified-synthesis.html) | Full 20-background audit, `A1–A6` postulates & Propositions 1–7 arbitration |
| [`evolution-and-audit.html`](evolution-and-audit.html) | Chronological logic evolution (`15 Oct – 26 Nov 2026`), red-team retraction ledger, full `G01–G24` status table & dated GitHub snapshots |
| [`bibliography.html`](bibliography.html) | Searchable working bibliography of 1,247 published foundational physics sources (Groups `01–11` & `13`; collected sources, not authored publications) |
| [`verification/verify_nogo_bounds.py`](verification/verify_nogo_bounds.py) | Rigidity no-go walls, SymPy Palatini $\omega=-3/2$ check, tidal crossover, optical clock & full-band LIGO O3b exclusion |
| [`verification/verify_planetary_node_bounds.py`](verification/verify_planetary_node_bounds.py) | 7-site thermodynamic, shielding, Arrhenius, causal light-time, cross-body AD ratios & Figure 2 generator |
| [`verification/verify_closure_theorems.py`](verification/verify_closure_theorems.py) | SymPy Jordan↔Einstein equivalence, ADM projection, Rayleigh $l=2$, $I_h$ $l=1\dots 5$ vanishing & Verlinde $F_{11}=89$ |
| [`verification/check_Ih_Y6.py`](verification/check_Ih_Y6.py) & [`checks_constants_rates.py`](verification/checks_constants_rates.py) | $I_h$ vertex $Y_6$ invariance, fundamental constants, DSW background & Sakharov $(\ell_P/a)^2$ barrier |
| [`verification/plot_ligo_o3b_bounds.py`](verification/plot_ligo_o3b_bounds.py) | Generates Figure 1 from 2,500-point LIGO H1/L1 O3b & aLIGO benchmark ASD files in `verification/data/` |
| [`verification/verification_output_2026-11-26.txt`](verification/verification_output_2026-11-26.txt) | Complete live execution log (5/5 suites + 2 figures, exit code 0) |

---

## Quick Reproduction

Requires Python 3.10+ with `numpy`, `sympy`, and `matplotlib`:

```bash
python3 verification/verify_nogo_bounds.py
python3 verification/verify_planetary_node_bounds.py
python3 verification/verify_closure_theorems.py
python3 verification/check_Ih_Y6.py
python3 verification/checks_constants_rates.py
python3 verification/plot_ligo_o3b_bounds.py
```

---

## Epistemic Status & Licensing

- **Epistemic Status:** Independent unpublished theoretical research and executable consistency suites (`G01–G24`: 5 closed, 11 in active progress, 8 open/conditional). All earlier hypotheses overturned by symbolic or empirical checks are transparently logged on [`evolution-and-audit.html`](evolution-and-audit.html).
- **License:** All Python scripts in `verification/` are released under the **MIT License** ([`LICENSE`](LICENSE)). Original text, PDF dossier, and generated figures by Xiaoyi Li are released under **CC BY 4.0**. Benchmark sensitivity curves in `verification/data/` follow standard public LIGO/GWOSC attribution.

#!/usr/bin/env python3
"""
VMAD 深海五大物理学前沿赛道与智慧节点量子物理颗粒度对齐 · 数值复算总脚本
对应专题文档（路线 B · 双轨并行）：
- 01-VMAD描述与分析/02-物理审计与验证方案/深海五大物理学前沿赛道与VMAD智慧节点量子对齐冷酷审计_2026-11-17.md
- 01-VMAD描述与分析/06-节点架构与场景推演/VMAD_深海五大物理学前沿与智慧节点量子物理极限推演报告.md

五大定量核验模块（含三级筛选边界与 4 项主动量子推演）：
1. 粒子与天体物理：KM3-230213A (220 PeV, DOC-0379) 双零锥 LIV 限制与定理七零退相干铁律、海冰分辨率/弦长对比、宇宙线 mu 子屏蔽本底；
2. 非线性层化流体力学：芒克 2 TW 配方红线与“底边界层正级联稀释 + 准地转逆级联 2.85 万年热隐身”闭合计算；
3. 深海物理声学与波导混沌：SOFAR 薛定谔同构、定理七位移退相干选址（SOFAR 水层 vs 深海玄武岩 + 阿基米德悬浮）与 10 Hz 肖尔特界面波屏蔽；
4. 极端温压凝聚态：超临界热液威多姆线卡诺自由功 (~89.4 MW) 与介电坍缩反射、万米海沟 110 MPa 静水围压外置箍缩红利；
5. 深海电动力学与光学：趋肤深度 10^-121 法拉第屏蔽、480 nm 蓝绿光窗 PLOB 量子容量归零定理与海底 MHD 动生感应反演。
"""

from __future__ import annotations

import json
import math
from pathlib import Path

# ---------------------------------------------------------------------------
# 基本物理常数（与 verify_closure_theorems.py / verify_antarctic_node_model.py 一致）
# ---------------------------------------------------------------------------
HBAR = 1.054571817e-34       # J*s
C = 2.99792458e8             # m/s
G0 = 6.67430e-11             # m^3/(kg*s^2)
KB = 1.380649e-23            # J/K
EV_TO_J = 1.602176634e-19    # J/eV
G_EARTH = 9.81               # m/s^2
MU0 = 4.0e-7 * math.pi       # H/m
YEAR_S = 365.25 * 86400.0    # s/yr
MPC_M = 3.085677581491367e22 # m/Mpc
M_PROTON = 1.67262192369e-27 # kg

LP = math.sqrt(HBAR * G0 / C**3)
MP = math.sqrt(HBAR * C / G0)
EP_J = MP * C**2


def gamma_ad_eta(m: float, dphi_ext: float, sigma_eps: float, tau_c: float, duration: float) -> float:
    """定理七式 (4-4)：近临界 eta 软模涨落导致的等效原理型异常退相干累积指数。"""
    x = duration / tau_c
    if x < 1e-6:
        bracket = 0.5 * x**2 - (x**3) / 6.0
    else:
        bracket = x - 1.0 + math.exp(-x)
    return ((m * dphi_ext / HBAR) ** 2) * (sigma_eps**2) * (tau_c**2) * bracket


def module_1_astroparticle_km3net_and_muon_shield() -> dict:
    print("=== [模块 1] 粒子与天体物理：KM3-230213A (220 PeV) 双零锥与深海宇宙线屏蔽 ===")
    E_nu_eV = 2.2e17  # 220 PeV (DOC-0379)
    E_nu_J = E_nu_eV * EV_TO_J
    m_nu_eV = 0.1
    m_nu_kg = m_nu_eV * EV_TO_J / C**2

    # 1.1 真空切伦科夫阈值限制超光速 LIV 参数与 LHC 质心能量对比
    delta_liv_kin_limit = (m_nu_eV / E_nu_eV) ** 2
    sqrt_s_TeV = math.sqrt(2.0 * M_PROTON * C**2 * E_nu_J) / (1.0e12 * EV_TO_J)

    # 1.2 定理七 220 PeV 自由落体中微子 (L=100 Mpc) 零退相干 vs 传统泡沫模型
    L_src = 100.0 * MPC_M
    T_flight = L_src / C
    M_gal = 1.5e12 * 1.98847e30
    R_gal = 1.0e4 * 3.08567758e16
    omega_tidal_sq = G0 * M_gal / (R_gal**3)
    dphi_tidal = omega_tidal_sq * (1.0**2)  # 夸张取 1 m 波包尺度
    gamma_vmad_km3 = gamma_ad_eta(m_nu_kg, dphi_tidal, sigma_eps=1.5e-17, tau_c=1.0, duration=T_flight)
    gamma_foam_km3 = (E_nu_J**2 / (HBAR * EP_J)) * T_flight

    # 1.3 海水 vs 南极冰光学散射长度与宇宙线 mu 子深度屏蔽本底 (Tier C 本底账本)
    lambda_s_sea_m = 60.0
    lambda_s_ice_m = 25.0
    scattering_ratio = lambda_s_sea_m / lambda_s_ice_m

    # 经验水深衰减模型 I_mu(h) / I_mu(0) = (1 + h/h0)^(-2) * exp(-h / Lambda_mu)
    h0 = 20.0
    lambda_mu = 1200.0
    i_surface = 100.0  # m^-2 s^-1 (垂直参考通量量级)
    depths = [0.0, 1000.0, 4000.0, 10000.0]
    muon_profile = []
    for h in depths:
        atten = ((1.0 + h / h0) ** -2.0) * math.exp(-h / lambda_mu)
        flux = i_surface * atten
        muon_profile.append({"depth_m": h, "attenuation_factor": atten, "approx_flux_m2_s": flux})
        print(f"  水深 {h:6.0f} m -> 宇宙线 mu 子衰减比={atten:.3e}, 通量~{flux:.3e} m^-2 s^-1")

    print(f"KM3-230213A (220 PeV): 质心能量 sqrt(s)={sqrt_s_TeV:.2f} TeV (超 LHC 14 TeV), 动力学 LIV 上限~{delta_liv_kin_limit:.2e}")
    print(f"  100 Mpc 自由落体定理七退相干={gamma_vmad_km3:.3e} vs 普朗克泡沫={gamma_foam_km3:.3e}")

    assert sqrt_s_TeV > 14.0
    assert gamma_vmad_km3 < 1e-75
    assert gamma_foam_km3 > 1e30
    assert muon_profile[2]["attenuation_factor"] < 1e-5
    assert muon_profile[3]["attenuation_factor"] < 1e-8

    return {
        "km3net_energy_PeV": 220.0,
        "cm_energy_TeV": sqrt_s_TeV,
        "kinematic_liv_limit": delta_liv_kin_limit,
        "vmad_theorem7_gamma_100Mpc": gamma_vmad_km3,
        "qg_foam_gamma_100Mpc": gamma_foam_km3,
        "scattering_length_sea_over_ice": scattering_ratio,
        "muon_shielding_profile": muon_profile,
    }


def module_2_stratified_turbulence_and_thermal_cloak() -> dict:
    print("\n=== [模块 2] 非线性层化流体：芒克红线与“正-逆双级联”2.85 万年热隐身 ===")
    kappa_munk = 1.0e-4   # m^2/s, 维持全球 AMOC 所需全水柱平均扩散率
    kappa_bg = 1.0e-5     # m^2/s, 深海内部实测背景扩散率
    P_munk_tw = 2.1       # TW, 潮汐/风场提供的机械功红线（严禁归因于节点）

    # 底边界层 (BBL) 内潮波 PSI 破碎区正级联快速稀释
    eps_bbl = 1.0e-6      # W/kg
    N_buoy = 1.0e-3       # s^-1 (深海典型浮力频率)
    gamma_mix = 0.2       # Osborn 混合效率系数
    L_ozmidov = math.sqrt(eps_bbl / (N_buoy**3))
    kappa_bbl = gamma_mix * eps_bbl / (N_buoy**2)

    # 100 MW 兰道尔排热在 500m x 500m 底边界层流场 (v=0.05 m/s) 下的近场温升
    P_node_w = 100.0e6
    rho_sw = 1025.0
    cp_sw = 3993.0
    area_bbl = 500.0 * 500.0
    v_bbl = 0.05
    dT_bbl = P_node_w / (rho_sw * cp_sw * area_bbl * v_bbl)

    # 离开 BBL 进入准地转层化区 (Ro<<1, Fr<<1): 动能逆级联水平聚集，垂直扩散被 kappa_bg 锁死
    H_thermo = 3000.0     # m
    tau_vert_s = (H_thermo**2) / kappa_bg
    tau_vert_yr = tau_vert_s / YEAR_S

    print(f"芒克红线: kappa_Munk={kappa_munk:.0e} m^2/s (~{P_munk_tw:.1f} TW 潮汐机械功) vs 内部背景 kappa_bg={kappa_bg:.0e} m^2/s")
    print(f"  底边界层 PSI 破碎区: L_O={L_ozmidov:.2f} m, kappa_BBL={kappa_bbl:.2f} m^2/s -> 100 MW 温升 dT={dT_bbl:.3e} K")
    print(f"  离开 BBL 后受 QG 逆级联与层化锁死，穿越 {H_thermo:.0f} m 水柱垂直扩散时标 = {tau_vert_yr:.3e} yr (~2.85 万年热隐身)")

    assert L_ozmidov > 20.0
    assert dT_bbl < 5e-3
    assert tau_vert_yr > 2.5e4

    return {
        "munk_required_diffusivity_m2_s": kappa_munk,
        "interior_bg_diffusivity_m2_s": kappa_bg,
        "munk_tidal_power_TW": P_munk_tw,
        "bbl_ozmidov_scale_m": L_ozmidov,
        "bbl_local_diffusivity_m2_s": kappa_bbl,
        "bbl_100MW_temp_rise_K": dT_bbl,
        "vertical_diffusion_3000m_yr": tau_vert_yr,
    }


def module_3_sofar_chaos_theorem7_and_scholte() -> dict:
    print("\n=== [模块 3] 物理声学与【主动推演 2】：SOFAR 退相干禁区 vs 基岩+肖尔特波 ===")
    # 3.1 SOFAR 薛定谔同构约化波长 (f=250 Hz, c0=1480 m/s)
    f_sofar = 250.0
    c0_sofar = 1480.0
    hbar_eff_m = c0_sofar / (2.0 * math.pi * f_sofar)

    # 3.2 [主动推演 2] 定理七位移退相干选址生死裁决
    dz_gamburtsev = 1.0e-11  # 南极甘布尔泽夫冷焊克拉通基准
    dz_basalt = 1.0e-9       # 深海玄武岩洋壳固体微震
    dz_sofar = 2.0           # 1,000 m SOFAR 水层内波垂直位移
    ratio_basalt_vs_gam = (dz_basalt / dz_gamburtsev) ** 2
    ratio_sofar_vs_gam = (dz_sofar / dz_gamburtsev) ** 2
    ratio_sofar_vs_basalt = (dz_sofar / dz_basalt) ** 2

    # 嵌套阿基米德超流/磁悬浮等密度平衡 (Delta_rho / rho = 1e-6) 对静态势差的压制
    buoy_mismatch = 1.0e-6
    archimedes_suppression = buoy_mismatch**2

    # 3.3 超低频 10 Hz 肖尔特界面波 (Scholte Wave) 水体侧指数衰减屏蔽
    f_sch = 10.0
    c_sch = 440.0    # m/s (~0.88 * c_s 沉积层剪切波速)
    c_w = 1500.0     # m/s
    k_sch = 2.0 * math.pi * f_sch / c_sch
    k_w = 2.0 * math.pi * f_sch / c_w
    gamma_w = math.sqrt(k_sch**2 - k_w**2)
    penetration_depth_m = 1.0 / gamma_w
    amp_at_100m = math.exp(-gamma_w * 100.0)
    amp_at_4000m = math.exp(-gamma_w * 4000.0)

    print(f"SOFAR 薛定谔同构 (250 Hz): 有效约化波长 lambdabar_0 = {hbar_eff_m:.3f} m")
    print(f"[主动推演 2] 定理七位移退相干对比: SOFAR水层 ({dz_sofar:.1f} m) 较玄武岩基岩 ({dz_basalt:.0e} m) 恶化 {ratio_sofar_vs_basalt:.2e} 倍 (禁区!)")
    print(f"  内部阿基米德等密度平衡 (drho/rho={buoy_mismatch:.0e}) 提供额外 {archimedes_suppression:.0e} 倍静态势差压制")
    print(f"10 Hz 肖尔特界面波 (c={c_sch:.0f} m/s): 水体侧 1/e 穿透尺度 = {penetration_depth_m:.2f} m -> 距底 100 m 衰减 {amp_at_100m:.2e}, 海面 4 km 为 {amp_at_4000m:.1e}")

    assert ratio_sofar_vs_gam > 1e22
    assert ratio_sofar_vs_basalt > 1e18
    assert penetration_depth_m < 10.0
    assert amp_at_4000m < 1e-200

    return {
        "sofar_effective_lambdabar_m": hbar_eff_m,
        "basalt_vs_gamburtsev_ratio": ratio_basalt_vs_gam,
        "sofar_vs_gamburtsev_ratio": ratio_sofar_vs_gam,
        "sofar_vs_basalt_ratio": ratio_sofar_vs_basalt,
        "archimedes_static_suppression": archimedes_suppression,
        "scholte_freq_Hz": f_sch,
        "scholte_speed_m_s": c_sch,
        "scholte_water_decay_scale_m": penetration_depth_m,
        "scholte_amp_at_100m": amp_at_100m,
        "scholte_amp_at_4000m": amp_at_4000m,
    }


def module_4_supercritical_widom_and_hydrostatic_confinement() -> dict:
    print("\n=== [模块 4] 极端凝聚态与【主动推演 1】：威多姆线卡诺机 & 110 MPa 围压箍缩 ===")
    # 4.1 超临界热液卡诺机 (Piccard 类环境: P=50 MPa > 22.06 MPa, T_H=680 K > 647.1 K)
    T_H = 680.0
    T_C = 274.65
    eta_carnot = 1.0 - (T_C / T_H)
    Q_vent_mw = 300.0
    realization = 0.50
    P_work_mw = Q_vent_mw * eta_carnot * realization

    # 跨越威多姆线介电坍缩 (eps_r: 80 -> 2.5)
    eps_liquid = 80.0
    eps_supercrit = 2.5
    n_liquid = math.sqrt(eps_liquid)
    n_supercrit = math.sqrt(eps_supercrit)
    R_fresnel = ((n_liquid - n_supercrit) / (n_liquid + n_supercrit)) ** 2

    # 4.2 [主动推演 1] 万米海沟 110 MPa 静水围压对定理四/五相腔的外置箍缩红利
    h_trench = 10900.0
    rho_sw_deep = 1025.0  # 仅取静水压近线性积分参考基准 (高压真实密度略高于 1025)
    P_hydro = rho_sw_deep * G_EARTH * h_trench  # ~1.096e8 Pa
    R_cav = 10.0
    F_confine_N = 4.0 * math.pi * (R_cav**2) * P_hydro

    # 若相腔内维持平衡体积储能密度的一个亚稳能标基准 Delta_rho_V = 1.0e8 J/m^3,
    # 在定理五 l=2 Rayleigh 极限 (Delta_rho_V > 0 时流体腔必失稳) 下，外部静水围压 P_hydro
    # 分担了总外约束压 (Delta_rho_eff = Delta_rho_V + P_hydro)，而在保持同样总约束压下，
    # 所需内禀相壁体张力/弹性模量负荷降低；反之若对比 1 MPa 浅海与 110 MPa 海沟，静水压差提供 100x 机械围压
    ratio_trench_to_100m = P_hydro / (rho_sw_deep * G_EARTH * 100.0)

    print(f"超临界热液卡诺机 (T_H={T_H:.0f} K, T_C={T_C:.2f} K): 卡诺效率={eta_carnot*100:.2f}%, 300 MW 热液可提自由功={P_work_mw:.2f} MW")
    print(f"  跨越威多姆线介电坍缩 (eps_r: {eps_liquid:.0f} -> {eps_supercrit:.1f}): 折射率 {n_liquid:.2f} -> {n_supercrit:.2f}, 单界面菲涅耳反射率={R_fresnel*100:.1f}%")
    print(f"[主动推演 1] 万米海沟静水围压 P_hydro={P_hydro/1e6:.2f} MPa (为 100 m 浅海的 {ratio_trench_to_100m:.0f} 倍):")
    print(f"  对 R=10 m 相腔直接施加各向同性机械箍缩力 F={F_confine_N:.3e} N, 等效抗膨胀能量密度红利 {P_hydro:.3e} J/m^3")

    assert 0.55 < eta_carnot < 0.65
    assert P_work_mw > 80.0
    assert P_hydro > 1.0e8

    return {
        "hydrothermal_TH_K": T_H,
        "deep_sea_TC_K": T_C,
        "carnot_efficiency": eta_carnot,
        "extractable_free_power_MW": P_work_mw,
        "widom_eps_liquid": eps_liquid,
        "widom_eps_supercritical": eps_supercrit,
        "widom_fresnel_power_reflection": R_fresnel,
        "trench_depth_m": h_trench,
        "hydrostatic_pressure_MPa": P_hydro / 1e6,
        "confinement_force_on_10m_cavity_N": F_confine_N,
    }


def module_5_mhd_skin_depth_and_plob_blue_green_limit() -> dict:
    print("\n=== [模块 5] 深海电动力学与【主动推演 3 & 4】：趋肤屏蔽、PLOB 归零与 MHD ===")
    sigma_sw = 4.0  # S/m
    z_deep = 4000.0

    def skin_depth(freq_hz: float) -> float:
        return 1.0 / math.sqrt(math.pi * freq_hz * MU0 * sigma_sw)

    # 5.1 趋肤深度法拉第笼铁律
    delta_spy1 = skin_depth(3.3e9)
    delta_elf = skin_depth(76.0)
    power_atten_1m_spy1 = math.exp(-2.0 * 1.0 / delta_spy1)
    power_atten_4km_elf = math.exp(-2.0 * z_deep / delta_elf)

    # 5.2 [主动推演 3] 480 nm 蓝绿光窗 PLOB 无中继量子纠缠容量极限
    c_ext = 0.04  # m^-1 (最纯净大洋深水衰减系数)
    plob_cases = []
    for L_m in [100.0, 1000.0, 4000.0]:
        eta_opt = math.exp(-c_ext * L_m)
        # 对极小 eta_opt 避免 1 - eta 下溢成 1.0
        q2_ebits_per_mode = eta_opt / math.log(2.0) if eta_opt < 1e-8 else -math.log2(1.0 - eta_opt)
        rate_at_1e15 = q2_ebits_per_mode * 1.0e15
        years_per_ebit = (1.0 / rate_at_1e15) / YEAR_S if rate_at_1e15 > 0 else float("inf")
        print(f"  蓝绿光窗 L={L_m:4.0f} m: 透射率={eta_opt:.3e} | PLOB Q2={q2_ebits_per_mode:.3e} ebits/mode | 10^15 Hz 下 1 ebit 需时={years_per_ebit:.2e} yr")
        plob_cases.append({
            "distance_m": L_m,
            "optical_transmission": eta_opt,
            "plob_q2_ebits_per_mode": q2_ebits_per_mode,
            "years_per_ebit_at_1e15Hz": years_per_ebit,
        })

    # 5.3 [主动推演 4] 深海动生 MHD 感应场量级 (关联 DOC-0356 - DOC-0363)
    v_flow = 0.1      # m/s
    L_flow = 1000.0   # m
    B0_earth = 4.5e-5 # T (45,000 nT)
    R_m = MU0 * sigma_sw * v_flow * L_flow
    B_ind_nT = R_m * B0_earth * 1e9

    print(f"导电海水趋肤深度: SPY-1B (3.3 GHz) delta={delta_spy1*1e3:.2f} mm (1 m 功率衰减={power_atten_1m_spy1:.2e}) | ELF (76 Hz) delta={delta_elf:.2f} m (4 km 功率衰减={power_atten_4km_elf:.2e})")
    print(f"[主动推演 4] 深海 MHD 动生感应 (v=0.1 m/s, L=1 km): 磁雷诺数 R_m={R_m:.3e}, 感生磁场幅值 B_ind~{B_ind_nT:.2f} nT (电中性 Fibonacci 任意子对此磁通免疫)")

    assert delta_spy1 < 0.01
    assert power_atten_4km_elf < 1e-100
    assert plob_cases[0]["plob_q2_ebits_per_mode"] > 1e-2
    assert plob_cases[2]["plob_q2_ebits_per_mode"] < 1e-65
    assert 10.0 < B_ind_nT < 50.0

    return {
        "seawater_conductivity_S_m": sigma_sw,
        "skin_depth_3p3GHz_m": delta_spy1,
        "power_atten_1m_3p3GHz": power_atten_1m_spy1,
        "skin_depth_76Hz_m": delta_elf,
        "power_atten_4000m_76Hz": power_atten_4km_elf,
        "blue_green_ext_coeff_m_minus_1": c_ext,
        "plob_optical_cases": plob_cases,
        "mhd_magnetic_reynolds": R_m,
        "mhd_induced_field_nT": B_ind_nT,
    }


def main() -> None:
    results = {
        "constants": {
            "lP_m": LP,
            "mP_kg": MP,
            "EP_J": EP_J,
        },
        "module_1_astroparticle": module_1_astroparticle_km3net_and_muon_shield(),
        "module_2_stratified_turbulence": module_2_stratified_turbulence_and_thermal_cloak(),
        "module_3_acoustics_and_theorem7": module_3_sofar_chaos_theorem7_and_scholte(),
        "module_4_supercritical_and_confinement": module_4_supercritical_widom_and_hydrostatic_confinement(),
        "module_5_mhd_and_plob": module_5_mhd_skin_depth_and_plob_blue_green_limit(),
    }
    out_path = Path(__file__).resolve().parent / "deepsea_physics_verification_results.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nALL ASSERTIONS PASSED. Results saved to: {out_path}")


if __name__ == "__main__":
    main()

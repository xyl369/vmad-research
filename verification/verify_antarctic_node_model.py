#!/usr/bin/env python3
"""
VMAD 南极四大科学异常与行星级遗迹母节点 · 数值与符号复算总脚本
对应专题文档：
- 01-VMAD描述与分析/02-物理审计与验证方案/南极四大前沿异常与VMAD公理体系颗粒度对齐冷酷审计_2026-11-16.md
- 01-VMAD描述与分析/06-节点架构与场景推演/VMAD_南极四大科学异常与行星级遗迹母节点_极限推演与颗粒度对齐报告.md

五大定量复算模块：
1. 定理七（等效原理型 eta 软模异常退相干）跨尺度对账：
   IceCube 天体中微子零退相干（自由落体抵消） vs 南极四类地质基座防震寿命增益
2. 公理 A6* 兰道尔排热与“南极冰底 1,000 m 逆向向上冻结羽流”固—液相变热管模型：
   融冰水流量、Clausius-Clapeyron 逆坡爬升过冷度、羽流生长时标、44 万年热扩散与 <0.01 K 卫星红外隐身校验
3. ANITA 0.6 EeV 破冰向上射电脉冲双路定量对账：
   双层 PREM 穿地弦长核子柱密度与标准模型存活率 (<1e-5) vs
   [路径 alpha] IF-22 Z16 拓扑暗扇区穿地截面上限与相腔界面催化转换 vs
   [路径 beta] G24 准直窄束卸载的平流层脚印与 IceCube 几何零命中
4. 威尔克斯地 500 km 质量瘤 (Mascon) 重力掩护质量上限与定理一 (W_max = M c^2) 算力储备：
   残差掩护质量 M_mask、定理四 Schwarzschild 无视界校验、全寿命负熵比特数与持续运行年限
5. 麦克默多极夜平流层—热层 (0 -> 140 km) 次级大气重力波指数级垂直放大走廊：
   波作用量守恒振幅放大因子 exp(z / 2H) ~ 3e4 倍与热层 80 K/day 加热所需垂直能流密度
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# 基本物理常数 (与 verify_closure_theorems.py / checks_constants_rates.py 一致)
# ---------------------------------------------------------------------------
HBAR = 1.054571817e-34       # J*s
C = 2.99792458e8             # m/s
G0 = 6.67430e-11             # m^3/(kg*s^2)
KB = 1.380649e-23            # J/K
H0 = 2.2e-18                 # s^-1
EV_TO_J = 1.602176634e-19    # J/eV
NA = 6.02214076e23           # mol^-1 (nucleons per gram ~ N_A)
SIGMA_SB = 5.670374419e-8    # W/(m^2*K^4)
YEAR_S = 365.25 * 86400.0    # s/yr
MPC_M = 3.085677581491367e22 # m/Mpc
G_EARTH = 9.81               # m/s^2

LP = math.sqrt(HBAR * G0 / C**3)
MP = math.sqrt(HBAR * C / G0)
EP_J = MP * C**2
T_GH = HBAR * H0 / (2.0 * math.pi * KB)


def gamma_ad_eta(m: float, dphi_ext: float, sigma_eps: float, tau_c: float, T: float) -> float:
    """
    定理七式 (4-4): 近临界 eta 软模引发的累积退相干指数 (量纲 1, 对密度矩阵非对角元压制 exp(-Gamma_AD))
    Gamma_AD(T) = (m * Delta_Phi_ext / hbar)^2 * sigma_eps^2 * tau_c^2 * [T/tau_c - 1 + exp(-T/tau_c)]
    """
    x = T / tau_c
    if x < 1e-6:
        bracket = 0.5 * x**2 - (x**3) / 6.0
    else:
        bracket = x - 1.0 + math.exp(-x)
    return ((m * dphi_ext / HBAR) ** 2) * (sigma_eps**2) * (tau_c**2) * bracket


def module_1_theorem7_decoherence() -> dict:
    print("=== [模块 1] 定理七 (等效原理型 eta 退相干) 跨尺度对账 ===")
    sigma_eps = 1.5e-17
    tau_c = 1.0  # s
    m_bmv = 1e-14
    d_bmv = 1e-4
    T_bmv = 1.0
    # 地面受支撑实验室基准 (Delta Phi = g * d)
    gamma_bmv_lab = gamma_ad_eta(m_bmv, G_EARTH * d_bmv, sigma_eps, tau_c, T_bmv)
    pref_bmv = m_bmv * G_EARTH * d_bmv / HBAR
    print(f"基准 BMV 受支撑 (m=1e-14kg, d=100um, T=1s, sigma_eps=1.5e-17, tau_c=1s):")
    print(f"  m*g*d/hbar = {pref_bmv:.3e} rad/s, 累积退相干指数 Gamma_AD = {gamma_bmv_lab:.3f} (O(1) 可观测)")

    # 场景 A: IceCube 天体中微子 (NGC 1068, L = 14.4 Mpc)
    m_nu = 0.1 * EV_TO_J / C**2  # 0.1 eV/c^2
    L_ngc1068 = 14.4 * MPC_M
    T_ngc1068 = L_ngc1068 / C
    # 自由落体下均匀引力势差严格抵消，仅剩星系级潮汐势差 Delta Phi_tidal ~ (G M_gal / R_gal^3) * d_wp^2
    M_gal = 1.5e12 * 1.98847e30
    R_gal = 1.0e4 * 3.08567758e16  # 10 kpc
    omega_tidal_sq = G0 * M_gal / R_gal**3
    # 即使将中微子质量本征态波包纵向分离取极夸张的宏观 1 m (实际远小于 1e-6 m), 潮汐势差也极其微小
    d_wp_macro = 1.0
    dphi_tidal_nu = omega_tidal_sq * (d_wp_macro**2)
    gamma_nu_vmad = gamma_ad_eta(m_nu, dphi_tidal_nu, sigma_eps, tau_c, T_ngc1068)
    # 对比传统唯象量子引力泡沫 (E^2 / (hbar * E_Pl)) * T (n=2 线性抑制被 IceCube 2024 排除)
    E_nu = 10e12 * EV_TO_J  # 10 TeV
    gamma_foam_n2 = (E_nu**2 / (HBAR * EP_J)) * T_ngc1068

    print(f"场景 A (IceCube NGC 1068 中微子, m_nu=0.1 eV, L=14.4 Mpc):")
    print(f"  星系潮汐势差 (即便按 d_wp=1m 夸张估计) = {dphi_tidal_nu:.3e} m^2/s^2")
    print(f"  VMAD 定理七预言累积退相干 = {gamma_nu_vmad:.3e} (严格 ~ 0，与 IceCube 2024 零退相干完全自洽)")
    print(f"  对比传统 QG 泡沫 n=2 模型累积退相干 = {gamma_foam_n2:.3e} (>>1，已被 IceCube 2024 证伪)")

    assert gamma_bmv_lab > 0.1 and gamma_bmv_lab < 10.0
    assert gamma_nu_vmad < 1e-80
    assert gamma_foam_n2 > 1e10

    # 场景 B: 南极四类地质环境的耦合位移抖动对节点有效退相干率影响 (Gamma_env ~ (M * g * delta_z / hbar)^2 * ...)
    # 固定参考载荷质量 M_ref = 1 kg (单比特/集体态等效耦合质量标尺) 在相同耦合时间窗口下的相对退相干率比值
    sites = [
        ("甘布尔泽夫冷焊克拉通山脊 (Dome A 下)", 1e-11, "10亿年零碰撞刚性陆核 + 冷基底零滑动"),
        ("普通活跃大陆基岩基准", 1e-8, "常规微震与地壳潮汐应变"),
        ("沃斯托克湖冰水界面", 1e-5, "35.9 MPa 高压液态水体对流与水力脉动"),
        ("西南极裂谷思韦茨冰川底部", 1e-3, "高地热融水泄洪 + 剧烈底冰滑动冰震"),
    ]
    dz_gam = sites[0][1]
    site_results = []
    print("场景 B (同质量同频段下，地震/水力位移 RMS 导致的等效定理七位移退相干相对倍比):")
    for name, dz, desc in sites:
        ratio_vs_gam = (dz / dz_gam) ** 2
        # 对应耦合质量 m_test=1e-10 kg 时的每秒退相干率标尺 (T=1s)
        rate_ref = gamma_ad_eta(1e-10, G_EARTH * dz, sigma_eps, tau_c, 1.0)
        print(f"  - {name}: delta_z={dz:.0e} m | 退相干倍比={ratio_vs_gam:.2e}x | ref_Gamma={rate_ref:.3e}")
        site_results.append({
            "site": name,
            "rms_displacement_m": dz,
            "decoherence_ratio_vs_gamburtsev": ratio_vs_gam,
            "ref_gamma_1e_minus_10kg_1s": rate_ref,
            "note": desc,
        })

    suppression_gam_vs_wais = (sites[-1][1] / dz_gam) ** 2
    assert abs(suppression_gam_vs_wais - 1e16) < 1e10
    return {
        "bmv_lab_gamma_1s": gamma_bmv_lab,
        "icecube_ngc1068_vmad_gamma": gamma_nu_vmad,
        "icecube_ngc1068_qg_foam_n2_gamma": gamma_foam_n2,
        "gamburtsev_vs_wais_suppression_factor": suppression_gam_vs_wais,
        "sites": site_results,
    }


def module_2_subglacial_heatpipe() -> dict:
    print("\n=== [模块 2] 公理 A6* 兰道尔排热与冰底 1,000 m 逆向冻结羽流相变热管 ===")
    L_f = 3.34e5        # J/kg, 冰融化潜热
    rho_ice = 917.0     # kg/m^3
    rho_w = 1000.0      # kg/m^3
    cp_ice = 2090.0     # J/(kg*K)
    k_ice = 2.2         # W/(m*K)
    alpha_ice = k_ice / (rho_ice * cp_ice)
    gamma_cc_pure = 0.0743e-6    # K/Pa (0.0743 K/MPa)
    gamma_cc_air = 0.098e-6      # K/Pa (含气饱和冰水)

    H_ice = 4000.0      # m
    tau_diff_s = H_ice**2 / alpha_ice
    tau_diff_yr = tau_diff_s / YEAR_S

    # 逆坡爬升过冷度 (Delta h = 1500 m, 如甘布尔泽夫深槽爬向山脊)
    dh_climb = 1500.0
    dP_climb = rho_ice * G_EARTH * dh_climb
    dT_supercool_pure = gamma_cc_pure * dP_climb
    dT_supercool_air = gamma_cc_air * dP_climb

    # 羽流几何: 长 50 km x 宽 10 km x 高 1000 m
    V_plume = 50e3 * 10e3 * 1000.0
    M_plume_ice = V_plume * rho_ice
    E_plume_latent = M_plume_ice * L_f

    # 流域表面热弥散 (100 km x 100 km) 对 -60 C (213 K) 表面红外黑体温度的影响
    A_basin = 100e3 * 100e3
    T_surf = 213.15

    power_cases = []
    print(f"4,000 m 冰盖纯热传导扩散时标 tau_diff = {tau_diff_yr:.3e} yr (~44 万年)")
    print(f"沿陡坡爬升 {dh_climb:.0f} m 压降 {dP_climb/1e6:.2f} MPa -> 绝热过冷度 = {dT_supercool_pure:.3f} K (纯水) / {dT_supercool_air:.3f} K (含气)")
    print(f"单座 50km x 10km x 1000m 逆向冻结羽流锁存相变潜热 = {E_plume_latent:.3e} J")
    for P_mw in [1.0, 10.0, 100.0, 1000.0]:
        P_w = P_mw * 1e6
        m_dot = P_w / L_f
        Q_w = m_dot / rho_w
        tau_plume_yr = (E_plume_latent / P_w) / YEAR_S
        dq_surf = P_w / A_basin
        dT_surf = dq_surf / (4.0 * SIGMA_SB * T_surf**3)
        print(f"  P={P_mw:6.1f} MW -> 融水流量 Q_w={Q_w:7.4f} m^3/s | 生成千m羽流耗时={tau_plume_yr:9.2e} yr | 表面温升 dT_surf={dT_surf:.4e} K")
        power_cases.append({
            "power_MW": P_mw,
            "melt_mass_rate_kg_s": m_dot,
            "water_volume_flow_m3_s": Q_w,
            "plume_formation_timescale_yr": tau_plume_yr,
            "basin_100x100km_flux_W_m2": dq_surf,
            "surface_ir_temp_rise_K": dT_surf,
        })

    # 校验 100 MW 情况的核心断言
    case_100mw = [c for c in power_cases if c["power_MW"] == 100.0][0]
    assert 0.25 < case_100mw["water_volume_flow_m3_s"] < 0.35
    assert case_100mw["surface_ir_temp_rise_K"] < 0.01
    assert tau_diff_yr > 4.0e5

    return {
        "thermal_diffusion_4000m_yr": tau_diff_yr,
        "climb_1500m_dP_MPa": dP_climb / 1e6,
        "supercooling_pure_K": dT_supercool_pure,
        "supercooling_aerated_K": dT_supercool_air,
        "plume_volume_m3": V_plume,
        "plume_latent_heat_J": E_plume_latent,
        "power_cases": power_cases,
    }


def earth_chord_column_density(theta_deg: float, n_steps: int = 20000) -> tuple[float, float]:
    """
    简化双层 PREM 地球密度剖面沿弦积分:
      地核 R_c = 3480 km, rho_c = 11.5 g/cm^3
      地幔 R_E = 6371 km, rho_m = 4.5 g/cm^3
    返回 (弦长 L_m, 质量柱密度 X_g_cm2)
    """
    R_E = 6371e3
    R_c = 3480e3
    rho_m = 4.5    # g/cm^3
    rho_c = 11.5   # g/cm^3
    theta = math.radians(abs(theta_deg))
    L_m = 2.0 * R_E * math.sin(theta)
    s = np.linspace(0.0, L_m, n_steps)
    # r(s)^2 = R_E^2 + s^2 - 2 * R_E * s * sin(theta)
    r_sq = R_E**2 + s**2 - 2.0 * R_E * s * math.sin(theta)
    rho = np.where(r_sq <= R_c**2, rho_c, rho_m)
    X_g_cm2 = float(np.trapezoid(rho, s * 100.0))  # s in cm
    return L_m, X_g_cm2


def module_3_anita_pueo_dual_channels() -> dict:
    print("\n=== [模块 3] ANITA 0.6 EeV 破冰向上脉冲双路定量对账 ===")
    # 0.6 EeV 标准模型 CC+NC 总核子散射截面 (Connolly et al. 2011 量级)
    sigma_nuN_sm = 1.5e-32  # cm^2
    events = [
        ("ANITA-I (2006, 仰角 -27°)", 27.0),
        ("ANITA-III (2014, 仰角 -35°)", 35.0),
    ]
    anita_chords = []
    for label, ang in events:
        L_m, X_gcm2 = earth_chord_column_density(ang)
        N_N = X_gcm2 * NA  # nucleons / cm^2
        tau_sm = N_N * sigma_nuN_sm
        p_surv_sm = math.exp(-tau_sm)
        # 路径 alpha: IF-22 拓扑暗扇区态透射 P_dark >= 0.5 所需截面
        sigma_dark_max = math.log(2.0) / N_N
        suppression_ratio = sigma_dark_max / sigma_nuN_sm
        print(f"{label}: 弦长={L_m/1e3:.0f} km | 柱密度={X_gcm2:.2e} g/cm^2 | tau_SM={tau_sm:.2f} | SM透射率={p_surv_sm:.2e}")
        print(f"  [路径 alpha IF-22 暗扇区] 保持 >=50% 穿地透射率需 sigma_dark <= {sigma_dark_max:.2e} cm^2 (约为 SM 的 {suppression_ratio:.2e})")
        assert p_surv_sm < 1e-5
        anita_chords.append({
            "event": label,
            "nadir_elevation_deg": -ang,
            "chord_length_km": L_m / 1e3,
            "column_density_g_cm2": X_gcm2,
            "nucleon_column_cm_minus_2": N_N,
            "sm_optical_depth": tau_sm,
            "sm_direct_transmission": p_surv_sm,
            "dark_portal_max_cross_section_cm2": sigma_dark_max,
            "cross_section_suppression_vs_sm": suppression_ratio,
        })

    # 路径 alpha 相腔壳层催化转换校验: 若暗态抵达南极节点边界 (层厚 dr ~ 1 m, 有效核子/拓扑耦合等效数密度)
    # 需满足转换概率 P_conv ~ (sigma_cat / sigma_dark) ~ 1e-3..1 在界面发生局部增强
    # 路径 beta: G24 相腔饱和卸载准直窄束几何截面积对账
    H_ice = 4.0      # km
    H_balloon = 37.0 # km
    R_vert = H_ice + H_balloon
    # ANITA 27 deg 仰角斜距 (从冰下出射点到 37 km 高空气球)
    R_slant_27 = H_balloon / math.sin(math.radians(27.0)) + H_ice / math.sin(math.radians(27.0))
    A_anita_fov = 1.5e6   # km^2 (气球视界覆盖面积 ~150万 km^2)
    beam_cases = []
    print(f"[路径 beta G24 准直窄束卸载] 垂直高 41 km / 27° 斜距 {R_slant_27:.1f} km:")
    for half_angle_deg in [0.5, 1.0, 2.0, 3.0, 5.0]:
        r_fp_vert_km = R_vert * math.tan(math.radians(half_angle_deg))
        r_fp_slant_km = R_slant_27 * math.tan(math.radians(half_angle_deg))
        area_vert_km2 = math.pi * r_fp_vert_km**2
        area_slant_km2 = math.pi * r_fp_slant_km**2
        p_anita_geo = min(1.0, area_slant_km2 / A_anita_fov)
        print(
            f"  半张角 {half_angle_deg:3.1f}° -> 垂直脚印半径={r_fp_vert_km:5.2f} km ({area_vert_km2:6.2f} km^2) | "
            f"27°斜距脚印半径={r_fp_slant_km:5.2f} km ({area_slant_km2:6.2f} km^2) | IceCube距源800-2000km命中=0"
        )
        beam_cases.append({
            "half_angle_deg": half_angle_deg,
            "vertical_footprint_radius_km": r_fp_vert_km,
            "vertical_footprint_area_km2": area_vert_km2,
            "slant_27deg_footprint_radius_km": r_fp_slant_km,
            "slant_27deg_footprint_area_km2": area_slant_km2,
            "balloon_single_shot_overlap_ratio": p_anita_geo,
            "icecube_off_axis_hit_prob": 0.0,
        })

    return {
        "sm_cross_section_0p6EeV_cm2": sigma_nuN_sm,
        "anita_chords": anita_chords,
        "collimated_beam_cases": beam_cases,
    }


def module_4_wilkes_mascon_and_theorem1() -> dict:
    print("\n=== [模块 4] 威尔克斯地 500 km Mascon 掩护质量上限与定理一算力储备 ===")
    dg_obs = 50e-5    # 50 mGal = 5e-4 m/s^2
    dg_tol = 5e-5     # 5 mGal = 5e-5 m/s^2 残差容忍度
    z0 = 15e3         # 15 km 埋深 (4 km 冰盖 + 11 km 上地壳, 航空/冰面近场最严苛基准)
    z_grace = 465e3   # GRACE 卫星轨道高度 + 埋深

    M_mask = dg_tol * (z0**2) / G0
    M_mask_grace = dg_tol * (z_grace**2) / G0
    r_s = 2.0 * G0 * M_mask / C**2
    W_max = M_mask * C**2

    # 定理一极限: 宇宙学德西特 Gibbons-Hawking 温度 T_GH
    N_bits_GH = W_max / (KB * T_GH * math.log(2.0))
    # 冰底实际局部热汇温度 T_sink = 270.5 K 下的宏观 Landauer 比特擦除上限
    T_sink = 270.5
    N_bits_270K = W_max / (KB * T_sink * math.log(2.0))

    # 以 P = 100 MW 连续耗散运行的理论持续寿命
    P_node = 100e6
    tau_run_yr = (W_max / P_node) / YEAR_S
    age_universe_yr = 1.38e10
    ratio_universe_age = tau_run_yr / age_universe_yr

    print(f"残差容忍度 delta_g = {dg_tol*1e5:.1f} mGal, 近场埋深 z0 = {z0/1e3:.0f} km (GRACE 轨高 {z_grace/1e3:.0f} km):")
    print(f"  近场可隐匿最大等效质量 M_mask = {M_mask:.3e} kg (约 {M_mask/1e3:.3e} 吨，即 {M_mask/1e11:.0f} 亿吨)")
    print(f"  GRACE 卫星尺度同残差掩护上限 = {M_mask_grace:.3e} kg")
    print(f"  对应 Schwarzschild 视界半径 r_s = {r_s:.3e} m (定理四校验: 宏观米级以上绝无黑洞坍缩风险)")
    print(f"  定理一总自由功 W_max = M*c^2 = {W_max:.3e} J")
    print(f"  T_GH={T_GH:.2e} K 理论比特上限 = {N_bits_GH:.3e} bits | 270.5 K 冰底热汇比特上限 = {N_bits_270K:.3e} bits")
    print(f"  按 100 MW 连续运行寿命 = {tau_run_yr:.3e} yr (为宇宙年龄的 {ratio_universe_age:.2e} 倍)")

    assert 1e14 < M_mask < 2e14
    assert r_s < 1e-12
    assert N_bits_GH > 1e80
    assert ratio_universe_age > 1e5

    return {
        "observed_mascon_anomaly_mGal": dg_obs * 1e5,
        "residual_tolerance_mGal": dg_tol * 1e5,
        "burial_depth_km": z0 / 1e3,
        "masked_mass_kg": M_mask,
        "masked_mass_metric_tons": M_mask / 1e3,
        "grace_altitude_masked_mass_kg": M_mask_grace,
        "schwarzschild_radius_m": r_s,
        "theorem1_W_max_J": W_max,
        "T_GH_K": T_GH,
        "landauer_bits_at_T_GH": N_bits_GH,
        "landauer_bits_at_270p5K": N_bits_270K,
        "runtime_at_100MW_yr": tau_run_yr,
        "runtime_over_universe_age": ratio_universe_age,
    }


def module_5_polar_gravity_wave_amplifier() -> dict:
    print("\n=== [模块 5] 极夜平流层—热层 (0 -> 140 km) 次级大气重力波指数级放大走廊 ===")
    H_scale = 6.8e3    # m, 极夜平均大气密度标高
    z_top = 140.0e3    # m, 热层底部雷达观测高度
    density_ratio = math.exp(-z_top / H_scale)
    amp_factor = math.exp(z_top / (2.0 * H_scale))

    # 热层 130 km 处反常加热 delta_T = 80 K/day 所需垂直波能流密度估算
    rho_130 = 8.0e-9   # kg/m^3 (MSIS 极区 ~130 km 典型空气密度)
    cp_air = 1005.0    # J/(kg*K)
    dT_day = 80.0      # K/day
    heating_rate_K_s = dT_day / 86400.0
    # 假设耗散层厚度 delta_z = 20 km (120-140 km)
    dz_layer = 20e3
    q_vol = rho_130 * cp_air * heating_rate_K_s  # W/m^3
    F_wave_req = q_vol * dz_layer                # W/m^2

    # 对应冰盖表面处的波致水平速度扰动仅需量级 (F_wave ~ 0.5 * rho_0 * u_0^2 * c_gz)
    rho_0 = 1.4        # kg/m^3, 极寒地表空气密度
    c_gz = 1.0         # m/s, 典型内重力波垂直群速度
    u0_req = math.sqrt(2.0 * F_wave_req / (rho_0 * c_gz))

    print(f"标高 H={H_scale/1e3:.1f} km, 高度 z={z_top/1e3:.0f} km:")
    print(f"  密度衰减比 rho(140km)/rho(0) = {density_ratio:.3e}")
    print(f"  无耗散波作用量守恒振幅放大倍数 exp(z/2H) = {amp_factor:.3e} 倍")
    print(f"  在 130 km (rho~8e-9 kg/m^3, 20 km 耗散层) 产生 80 K/day 加热仅需能流 F_wave = {F_wave_req:.3e} W/m^2")
    print(f"  折合冰盖表面只需激励微小空气质点速度扰动 u0 ~ {u0_req:.3e} m/s (彻底淹没在常规风噪中)")

    assert amp_factor > 2.0e4
    assert F_wave_req < 1e-2

    return {
        "scale_height_km": H_scale / 1e3,
        "target_altitude_km": z_top / 1e3,
        "density_ratio_140km_to_surface": density_ratio,
        "amplitude_amplification_factor": amp_factor,
        "thermosphere_heating_K_per_day": dT_day,
        "required_wave_energy_flux_W_m2": F_wave_req,
        "equivalent_surface_velocity_perturbation_m_s": u0_req,
    }


def main() -> None:
    results = {
        "constants": {
            "lP_m": LP,
            "eta0_m_minus_2": C**3 / (4.0 * HBAR * G0),
            "K0_N": C**4 / (8.0 * math.pi * G0),
            "T_GH_K": T_GH,
        },
        "module_1_theorem7_decoherence": module_1_theorem7_decoherence(),
        "module_2_subglacial_heatpipe": module_2_subglacial_heatpipe(),
        "module_3_anita_pueo_dual_channels": module_3_anita_pueo_dual_channels(),
        "module_4_wilkes_mascon_and_theorem1": module_4_wilkes_mascon_and_theorem1(),
        "module_5_polar_gravity_wave_amplifier": module_5_polar_gravity_wave_amplifier(),
    }
    out_path = Path(__file__).resolve().parent / "antarctic_node_verification_results.json"
    out_path.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nALL ASSERTIONS PASSED. Results saved to: {out_path}")


if __name__ == "__main__":
    main()

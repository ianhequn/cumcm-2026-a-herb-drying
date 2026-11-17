"""问题四：附录 4 物性 + 附件 2 实测收缩 R(t) 下的烘干时长与水分场。

相对问题三只新增两件事：物性闭合换成附录 4、半径由附件 2 插值（不自拟收缩律）。
s=(r/R(t))^2 坐标把收缩全部收进几何系数 4/R(t)^2，无需移动网格重映射。
产出：t_dry、表 6（末列"药材表面"）、result4.xlsx、四组合消融、体积反演一致性核查。
"""
from __future__ import annotations

import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import numpy as np

from params import *  # noqa: F401,F403  口径常数（显式声明来源，勿依赖 utils 转出）
from utils import *


def volume_inversion(radius, times_s):
    """体积反演一致性核查：由附件 2 的 R(t) 与附录 4 的密度式反推平均含水率。

    干质量线密度 K = M_s/(pi L) = rho_4(C0) R0^2/(1+C0) 恒定，故
      (b0 + b1 C) R^2 = (1+C) K  =>  C = (K - b0 R^2)/(b1 R^2 - K)，
    系数 b0/b1 直接取自 PROBLEM_FACTS 的附录 4 密度式，不在此重打字面量。
    """
    a4 = FACTS['appendix4_params']
    b0, b1 = float(a4['rho_intercept']), float(a4['rho_slope_on_C'])
    K = DRY_MASS_LINE_A4
    R = np.array([float(radius.R_of(t)) for t in times_s])
    C_vol = (K - b0 * R ** 2) / (b1 * R ** 2 - K)
    V_ratio = (R / R0_M) ** 2
    return R, C_vol, V_ratio


def volume_additive_C(radius, t_query):
    """体积可加性对照闭合：V/V0 = (1+kappa C)/(1+kappa C0)，kappa=rho_s/rho_w。"""
    V_ratio = (float(radius.R_of(t_query)) / R0_M) ** 2
    kappa = VOL_ADD_KAPPA
    return (V_ratio * (1.0 + kappa * C0_KGKG) - 1.0) / kappa, kappa, V_ratio


def validate_capability(payload, res, radius, abl):
    """P4-C1 / P4-C2 / P4-C3 / P4-C4 / P4-C5 的可伪证检查。"""
    props4 = props_appendix4()
    assert radius.n_points == N_A2_POINTS, \
        f'R7/P4-C1：附件 2 应读入 {N_A2_POINTS} 点，实为 {radius.n_points}'
    assert abs(float(radius.R_of(0.0)) - R0_M) < EXACT_TOL, 'P4-C1：R(0) 须等于 R0'
    assert float(radius.R_of(Q34_TABLE_DT_H * SEC_PER_HOUR)) < R_AT_6H_MAX_M, \
        f'P4-C1：R(6 h) 应已收缩到 {R_AT_6H_MAX_M * CM_PER_M:.2f} cm 以下'
    assert float(np.ptp(res['R'])) > 0.0, 'P4-C1：R(t) 未随时间变化，收缩没接进求解'
    assert float(props4.D(C0_KGKG, T0_C)) < float(props_appendix3().D(C0_KGKG, T0_C)), \
        'P4-C1：附录 4 的 D 应小于附录 3'
    assert float(np.max(res['C'][-1])) < C_TH, 'P4-C2：末步 max_r C 未低于阈值'
    assert float(res['Cmax'][IDX_PREV_STEP]) >= C_TH, 'P4-C2：前一输出步应仍未达标'
    assert res['max_pos_grad'] < GRAD_POS_TOL, f'R12：出现正径向梯度 {res["max_pos_grad"]:.3e}'
    assert DURATION_PRIOR_H[0] <= payload['t_dry_h'] <= DURATION_PRIOR_H[1], \
        f'R10：t_dry={payload["t_dry_h"]:.4f} h 落在 2-3 天先验之外'
    assert payload['t_dry_h'] < payload['ablation']['baseline_h'], \
        'P4-C4：问题四时长应短于问题三基准'
    tab = payload['table6']
    assert tab['columns'][-1] == Q4_SURFACE_LABEL, 'P4-C2：表 6 末列须为"药材表面"'
    assert all(r <= min(tab['R_cm']) + GRAD_POS_TOL for r in Q4_TABLE_RADII_CM), \
        'P4-C2：表 6 出现超过当时 R(t) 的固定半径列'
    assert tab['rows'][-1]['label'].startswith('烘干结束'), 'P4-C2：表 6 末行须为烘干结束时间'
    assert payload['xlsx_blank_cells'] > 0, 'P4-C3：超出 R(t) 的位置应留空，不得域外插值'
    assert payload['xlsx_last_header'] == Q4_SURFACE_LABEL, 'P4-C3：result4.xlsx 末列须为药材表面'
    assert abl['props_only_pct'] > ABL_PROPS_ONLY_MIN_PCT, \
        f'P4-C4：仅换物性 {abl["props_only_pct"]:+.1f}% 未超 +{ABL_PROPS_ONLY_MIN_PCT}%'
    assert abl['shrink_only_pct'] < ABL_SHRINK_ONLY_MAX_PCT, \
        f'P4-C4：仅加收缩 {abl["shrink_only_pct"]:+.1f}% 未低于 {ABL_SHRINK_ONLY_MAX_PCT}%'
    assert abl['props_only_pct'] * abl['shrink_only_pct'] < 0.0, 'P4-C4：两主效应未反号'
    assert abs(abl['interaction_h']) > 1.0, 'P4-C4：交互项过小，与"强耦合"结论不符'
    vi = payload['volume_inversion']
    assert any(v < 0.0 for v in vi['C_vol']), 'P4-C5：体积反演未出现 C<0（一致性核查失效）'
    assert vi['first_negative_h'] is not None and vi['first_negative_h'] <= HOURS_PER_DAY, \
        f'P4-C5：C<0 应在 {HOURS_PER_DAY:.0f} h 及之前出现'
    assert vi['dev_at_18h'] > VOL_INV_MIN_DEV, \
        (f'P4-C5：{VOL_INV_PROBE_H:.0f} h 偏离 {vi["dev_at_18h"] * PCT:.1f}% '
         f'未超 {VOL_INV_MIN_DEV * PCT:.0f}%')
    assert res['cons_res'] < CONS_RES_TOL, '干质量守恒残差超限'
    grid = payload['grid_convergence']
    assert grid['rel_diff'][-1] < GRID_REL_TOL, 'P4-C4：网格相对差未达门限'
    echo('[capability] P4-C1..P4-C5 运行时检查全部通过')


def main() -> int:
    set_all_seeds()
    t_start = time.time()
    props = props_appendix4()
    air = AirState()
    radius = RadiusHistory()
    echo(f'[Q4] 附录 4：rho(2.55)={float(props.rho(C0_KGKG, T0_C)):.1f}, '
         f'cp(2.55)={float(props.cp(C0_KGKG, T0_C)):.1f}, k(2.55)={float(props.k(C0_KGKG, T0_C)):.5f}, '
         f'D(2.55,28C)={float(props.D(C0_KGKG, T0_C)):.4e}')
    t_6h, t_end_h = Q34_TABLE_DT_H * SEC_PER_HOUR, DURATION_PRIOR_H[1] * SEC_PER_HOUR
    echo(f'[Q4] 附件 2：{radius.n_points} 点，R(0)={float(radius.R_of(0.0)) * CM_PER_M:.4f} cm，'
         f'R({Q34_TABLE_DT_H:.0f}h)={float(radius.R_of(t_6h)) * CM_PER_M:.4f} cm，'
         f'R({DURATION_PRIOR_H[1]:.0f}h)={float(radius.R_of(t_end_h)) * CM_PER_M:.4f} cm')

    res = solve_transport(props, air, radius, N_PROD, DT_Q4,
                          DURATION_PRIOR_H[1] * 3.0 * SEC_PER_HOUR,
                          float(Q34_FILE_DT_S), C_th=C_TH)
    assert res['t_dry'] is not None, '阈值未在时窗内达成'
    t_dry_s = float(res['t_dry'])
    t_dry_h = t_dry_s / SEC_PER_HOUR
    echo(f'[Q4] 生产网格 N={N_PROD}, dt={DT_Q4} s：t_dry={t_dry_h:.6f} h ({t_dry_h / HOURS_PER_DAY:.4f} d)，'
         f'步数={res["n_steps"]}，守恒残差={res["cons_res"]:.3e}，耗时={time.time() - t_start:.1f}s')
    echo(f'[Q4] 末步 max_r C={float(np.max(res["C"][-1])):.6f} < {C_TH}，'
         f'前一步 {float(res["Cmax"][IDX_PREV_STEP]):.6f}；'
         f'最大正梯度={res["max_pos_grad"]:.3e}')

    # ---- 表 6：每 6 h 一行 + 烘干结束行；固定半径列只保留全程 <= min R(t) 的三列 ----
    step_6h = int(round(Q34_TABLE_DT_H * SEC_PER_HOUR / Q34_FILE_DT_S))
    C_fix = field_at_radii(res, Q4_TABLE_RADII_CM, 'C')
    tab_rows = []
    for kk in range(1, int(np.floor(t_dry_h / Q34_TABLE_DT_H)) + 1):
        idx = kk * step_6h
        tab_rows.append({'label': f'{kk * Q34_TABLE_DT_H:g}',
                         't_h': kk * Q34_TABLE_DT_H,
                         'R_cm': float(res['R'][idx]) * CM_PER_M,
                         'C_fixed': [float(v) for v in C_fix[idx]],
                         'C_surface': float(res['C'][idx][-1]),
                         'Cbar': float(res['Cbar'][idx])})
    C_end, R_end = snapshot_at_time(res, t_dry_s, 'C')
    T_end, _ = snapshot_at_time(res, t_dry_s, 'T')
    s_q = [min((r / (R_end * CM_PER_M)) ** 2, 1.0) for r in Q4_TABLE_RADII_CM]
    tab_rows.append({'label': f'烘干结束 {t_dry_h:.4f}', 't_h': t_dry_h,
                     'R_cm': R_end * CM_PER_M,
                     'C_fixed': [float(np.interp(v, res['s'], C_end)) for v in s_q],
                     'C_surface': float(C_end[-1]),
                     'Cbar': float(res['m'] @ C_end)})
    for row in tab_rows[:N_ECHO_HEAD] + tab_rows[-N_ECHO_TAIL:]:
        echo(f'[Q4] 表6 t={row["label"]:>24s} h  R={row["R_cm"]:.4f} cm  C=' +
             ' '.join(f'{v:.4f}' for v in row['C_fixed']) +
             f'  表面={row["C_surface"]:.4f}  Cbar={row["Cbar"]:.4f}')

    # ---- result4.xlsx：0.1 cm 步长的物理半径列 + 末列"药材表面"，超出 R(t) 留空 ----
    radii_file = file_radii_cm(R0_M * CM_PER_M)
    Cf = field_at_radii(res, radii_file, 'C')
    surf = surface_value(res, 'C').reshape(-1, 1)
    matrix = np.hstack([Cf, surf])
    header = [HEADER_CORNER] + radii_file + [Q4_SURFACE_LABEL]
    times_file = [int(v) for v in res['t'][1:]]
    n_blank = int(np.count_nonzero(np.isnan(matrix[1:])))
    write_result_xlsx(os.path.join(ROOT, 'result4.xlsx'),
                      [('Sheet1', header, times_file, matrix[1:])])
    echo(f'[Q4] result4.xlsx 已写：1 表 × {len(times_file)} 行 × {len(header)} 列，'
         f'超出 R(t) 留空 {n_blank} 格（占 {n_blank / (len(times_file) * len(radii_file)) * 100:.1f}%）')

    # ---- 四组合消融（真实重算，N=N_ABL）----
    p3 = props_appendix3()
    fixed = FixedRadius()
    abl_defs = [('baseline_h', '基准（问题三）', p3, fixed),
                ('props_only_h', '仅换物性（附录 4 + R=R0）', props, fixed),
                ('shrink_only_h', '仅加收缩（附录 3 + R(t)）', p3, radius),
                ('full_h', '问题四（全）', props, radius)]
    abl = {}
    for key, label, pp, rr in abl_defs:
        t0 = time.time()
        th, _ = t_dry_hours(pp, air, rr, N_ABL, DT_Q4)
        abl[key] = th
        echo(f'[Q4] 消融 {label:<24s} t_dry={th:.4f} h  [{time.time() - t0:.1f}s]')
    base = abl['baseline_h']
    abl['props_only_pct'] = (abl['props_only_h'] - base) / base * 100.0
    abl['shrink_only_pct'] = (abl['shrink_only_h'] - base) / base * 100.0
    abl['full_pct'] = (abl['full_h'] - base) / base * 100.0
    abl['interaction_h'] = (abl['full_h'] - base) - (abl['props_only_h'] - base) \
        - (abl['shrink_only_h'] - base)
    echo(f'[Q4] 消融相对基准：仅换物性 {abl["props_only_pct"]:+.1f}%，'
         f'仅加收缩 {abl["shrink_only_pct"]:+.1f}%，全 {abl["full_pct"]:+.1f}%；'
         f'交互项 {abl["interaction_h"]:+.4f} h')

    # ---- 网格无关性 ----
    grid_N, grid_t = list(N_CONV_Q4), []
    for n in grid_N:
        th, _ = t_dry_hours(props, air, radius, n, DT_Q4)
        grid_t.append(th)
        echo(f'[Q4] 网格自校验 N={n:3d}: t_dry={th:.6f} h')
    rel_diff = [abs(grid_t[i] - grid_t[i - 1]) / grid_t[i] for i in range(1, len(grid_t))]
    echo('[Q4] 相邻相对差 ' + ', '.join(f'{d:.3e}' for d in rel_diff))

    # ---- §8.4 体积反演一致性核查 ----
    vi_times_h = [0.0] + [k * Q34_TABLE_DT_H for k in range(1, N_VI_INTERVALS + 1)]
    vi_times_s = [t * SEC_PER_HOUR for t in vi_times_h]
    R_vi, C_vol, V_ratio = volume_inversion(radius, vi_times_s)
    neg = [t for t, c in zip(vi_times_h, C_vol) if c < 0.0]
    first_neg = min(neg) if neg else None
    t_probe_s = VOL_INV_PROBE_H * SEC_PER_HOUR
    idx18 = int(round(t_probe_s / Q34_FILE_DT_S))
    Cbar_18h_transport = float(res['Cbar'][idx18])
    C_va, kappa, V18 = volume_additive_C(radius, t_probe_s)
    dev18 = abs(C_va - Cbar_18h_transport) / Cbar_18h_transport
    echo(f'[Q4] 体积反演 Cbar_vol(每{Q34_TABLE_DT_H:.0f}h) = ' +
         ' '.join(f'{v:+.4f}' for v in C_vol))
    echo(f'[Q4] 首次出现 C<0 于 t={first_neg:.0f} h（物理不可能）；'
         f'{DURATION_PRIOR_H[1]:.0f} h 处 {C_vol[-1]:+.4f}')
    echo(f'[Q4] 体积可加性对照 kappa=rho_s/rho_w={kappa:.4f}（1/kappa={1 / kappa:.4f}）：'
         f'C({VOL_INV_PROBE_H:.0f}h)={C_va:.4f} vs 输运 {Cbar_18h_transport:.4f}，'
         f'偏离 {dev18 * PCT:.1f}%')

    # ---- 附录 4 密度式的可行性根：C 使 rho_4(C)R^2=(1+C)K 在 R=R_min 时的解 ----
    a4 = FACTS['appendix4_params']
    b0, b1 = float(a4['rho_intercept']), float(a4['rho_slope_on_C'])
    R_min = float(radius.R_of(DURATION_PRIOR_H[1] * SEC_PER_HOUR))
    C_root = (DRY_MASS_LINE_A4 - b0 * R_min ** 2) / (b1 * R_min ** 2 - DRY_MASS_LINE_A4)
    echo(f'[Q4] 末态 R={R_min * CM_PER_M:.4f} cm 对应的密度式根 C={C_root:+.4f}，'
         f'负值说明附件 2 的收缩幅度大于附录 4 密度式所允许的范围')

    payload = {
        'meta': {'problem': 4, 'props': props.name, 'n_cells': N_PROD, 'dt_s': DT_Q4,
                 'shrinkage_source': SHRINKAGE_SOURCE, 'decimals': DECIMALS,
                 'criterion': f'max_r C(r,t) < {C_TH} (pointwise)',
                 'runtime_s': round(time.time() - t_start, 2)},
        't_dry_h': t_dry_h, 't_dry_s': t_dry_s, 't_dry_d': t_dry_h / HOURS_PER_DAY,
        't_dry_reported': round(t_dry_h, DECIMALS),   # 与 RESULTS.md 引用的四位小数同口径
        'C_max_at_end': float(np.max(res['C'][-1])),
        'C_max_prev_step': float(res['Cmax'][IDX_PREV_STEP]),
        'max_pos_grad': res['max_pos_grad'],
        'radius_history': {'n_points': radius.n_points,
                           'R0_cm': float(radius.R_of(0.0)) * CM_PER_M,
                           'R_6h_cm': float(radius.R_of(6 * SEC_PER_HOUR)) * CM_PER_M,
                           'R_min_cm': float(np.min(radius.R)) * CM_PER_M,
                           'R_at_tdry_cm': R_end * CM_PER_M,
                           't_s': to_list(radius.t), 'R_cm': to_list(radius.R * CM_PER_M, 4)},
        'table6': {'columns': [f'{r:g} cm' for r in Q4_TABLE_RADII_CM] + [Q4_SURFACE_LABEL],
                   'radii_cm': list(Q4_TABLE_RADII_CM),
                   'R_cm': [r['R_cm'] for r in tab_rows],
                   'rows': [{'label': r['label'], 't_h': r['t_h'],
                             'R_cm': round(r['R_cm'], DECIMALS),
                             'C_fixed': [round(v, DECIMALS) for v in r['C_fixed']],
                             'C_surface': round(r['C_surface'], DECIMALS),
                             'Cbar': round(r['Cbar'], DECIMALS)} for r in tab_rows]},
        'T_at_tdry': to_list(T_end, DECIMALS),
        'series_t_h': to_list(res['t'] / SEC_PER_HOUR, SERIES_DECIMALS),
        'series_R_cm': to_list(res['R'] * CM_PER_M, SERIES_DECIMALS),
        'series_Cmax': to_list(res['Cmax'], SERIES_DECIMALS),
        'series_Cbar': to_list(res['Cbar'], SERIES_DECIMALS),
        'series_Csurf': to_list(surface_value(res, 'C'), SERIES_DECIMALS),
        'series_Tcenter': to_list(res['T'][:, 0], 4),
        'series_Tsurf': to_list(res['T'][:, -1], 4),
        'file_radii_cm': radii_file,
        'file_last_column': Q4_SURFACE_LABEL,
        'field_C_full': to_list(Cf[::JSON_STRIDE], DECIMALS),
        'field_C_full_stride': JSON_STRIDE,
        'ablation': abl,
        'grid_convergence': {'N': grid_N, 't_dry_h': grid_t, 'rel_diff': rel_diff},
        'volume_inversion': {'times_h': vi_times_h, 'R_cm': to_list(R_vi * CM_PER_M, 4),
                             'C_vol': to_list(C_vol, 4), 'V_ratio': to_list(V_ratio, 4),
                             'first_negative_h': first_neg,
                             'dry_mass_line_kg_per_m': DRY_MASS_LINE_A4,
                             'kappa_rho_s_over_rho_w': kappa,
                             'C_volume_additive_18h': C_va,
                             'Cbar_transport_18h': Cbar_18h_transport,
                             'V_ratio_18h': V18,
                             'dev_at_18h': dev18,
                             'density_feasible_root_at_Rmin': C_root},
        'conservation_residual': res['cons_res'],
        'xlsx_rows': len(times_file),
        'xlsx_cols': len(header),
        'xlsx_blank_cells': n_blank,
        'xlsx_last_header': header[-1],
    }
    validate_capability(payload, res, radius, abl)
    dump_json(os.path.join(FIG_DIR, 'problem_4_results.json'), payload)
    dump_json(os.path.join(FIG_DIR, 'q4_results.json'), payload)
    echo(f'[Q4] DONE t_dry={t_dry_h:.4f} h（{t_dry_h / HOURS_PER_DAY:.4f} d），'
         f'总耗时 {time.time() - t_start:.1f}s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

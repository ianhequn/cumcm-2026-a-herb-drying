"""问题三：烘干时长（逐点阈值首次到达反问题）。

与问题二同一个初边值问题（附录 3 全程、R=R0），唯一新增内容是终止判据
t_dry = min{ t : max_r C(r,t) < 0.15 }，穿越点在相邻 60 s 输出步间线性插值。
产出：t_dry、表 5、result3.xlsx（60 s × 0.1 cm）、网格与时间步无关性表。
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


def validate_capability(payload, res, t_dry_h):
    """P3-C1 / P3-C3 / P3-C4 的可伪证检查（逐点判据 + 收敛性）。"""
    C_profile = res['C'][-1]
    assert float(np.max(C_profile)) < C_TH, f'P3-C1：末步 max_r C 未低于 {C_TH}'
    assert float(res['Cmax'][IDX_PREV_STEP]) >= C_TH, 'P3-C1：前一输出步应仍未达标（首次穿越）'
    assert abs(float(np.max(C_profile)) - float(C_profile[0])) < GRAD_POS_TOL, (
        'P3-C1：径向最大值须等于中心值（判据等价性）')
    assert res['max_pos_grad'] < GRAD_POS_TOL, (
        f'R12/P3-C1：出现正径向梯度 {res["max_pos_grad"]:.3e}，判据等价性不成立')
    assert DURATION_PRIOR_H[0] <= t_dry_h <= DURATION_PRIOR_H[1], (
        f'R10/P2-C4：t_dry={t_dry_h:.4f} h 落在题面 2-3 天先验之外')
    Cbar_at_36h = payload['Cbar_at_36h']
    assert Cbar_at_36h < C_TH, '均值判据对照：Cbar 在制表探针时刻已低于阈值（说明均值判据会低估）'
    grid = payload['grid_convergence']
    diffs = [d for d in grid['rel_diff'] if d is not None]
    assert all(diffs[i] > diffs[i + 1] for i in range(len(diffs) - 1)), (
        'P3-C4：相邻网格相对差未单调下降')
    assert all(d < GRID_REL_TOL for d in diffs[IDX_PREV_STEP:]), 'P3-C4：末两组相对差未达门限'
    assert payload['grid_rel_diff'] < GRID_REL_TOL, \
        f'R13：生产网格相对差 {payload["grid_rel_diff"]:.3e}'
    assert len(set(grid['t_dry_h'])) == len(grid['t_dry_h']), 'P3-C4：网格表数值重复，疑非真实重算'
    assert res['cons_res'] < CONS_RES_TOL, '干质量守恒残差超限'
    echo('[capability] P3-C1/P3-C2/P3-C3/P3-C4 运行时检查全部通过')


def main() -> int:
    set_all_seeds()
    t_start = time.time()
    props = props_appendix3()
    air = AirState()
    radius = FixedRadius()

    res = solve_transport(props, air, radius, N_PROD, DT_Q3,
                          DURATION_PRIOR_H[1] * 3.0 * SEC_PER_HOUR,
                          float(Q34_FILE_DT_S), C_th=C_TH)
    assert res['t_dry'] is not None, '阈值未在时窗内达成'
    t_dry_s = float(res['t_dry'])
    t_dry_h = t_dry_s / SEC_PER_HOUR
    echo(f'[Q3] 生产网格 N={N_PROD}, dt={DT_Q3} s：t_dry={t_dry_h:.6f} h '
         f'({t_dry_h / HOURS_PER_DAY:.4f} d)，步数={res["n_steps"]}，'
         f'守恒残差={res["cons_res"]:.3e}，'
         f'耗时={time.time() - t_start:.1f}s')
    echo(f'[Q3] 末步 max_r C={float(np.max(res["C"][-1])):.6f} < {C_TH}，'
         f'前一步 {float(res["Cmax"][IDX_PREV_STEP]):.6f}；'
         f'径向最大正梯度={res["max_pos_grad"]:.3e}')

    table_radii = [round(k * Q34_TABLE_DR_CM, 4)
                   for k in range(int(R0_M * CM_PER_M / Q34_TABLE_DR_CM) + 1)]
    C_tab_all = field_at_radii(res, table_radii, 'C')
    T_tab_all = field_at_radii(res, table_radii, 'T')
    step_6h = int(round(Q34_TABLE_DT_H * SEC_PER_HOUR / Q34_FILE_DT_S))
    rows_t_h, rows_C, rows_Cbar, rows_T = [], [], [], []
    for k in range(1, int(np.floor(t_dry_h / Q34_TABLE_DT_H)) + 1):
        idx = k * step_6h
        rows_t_h.append(k * Q34_TABLE_DT_H)
        rows_C.append(C_tab_all[idx])
        rows_T.append(T_tab_all[idx])
        rows_Cbar.append(float(res['Cbar'][idx]))
    C_end, _ = snapshot_at_time(res, t_dry_s, 'C')
    T_end, _ = snapshot_at_time(res, t_dry_s, 'T')
    Cbar_end = float(res['m'] @ C_end)
    rows_t_h.append(t_dry_h)
    rows_C.append(np.array([np.interp(min((r / (R0_M * CM_PER_M)) ** 2, 1.0), res['s'], C_end)
                            for r in table_radii]))
    rows_T.append(np.array([np.interp(min((r / (R0_M * CM_PER_M)) ** 2, 1.0), res['s'], T_end)
                            for r in table_radii]))
    rows_Cbar.append(Cbar_end)
    table5 = np.asarray(rows_C)
    echo('[Q3] 表 5 首行 (6 h) C = ' + ' '.join(f'{v:.4f}' for v in table5[0]) +
         f' | Cbar={rows_Cbar[0]:.4f}')
    echo(f'[Q3] 表 5 末行 (烘干结束 {t_dry_h:.4f} h) C = ' +
         ' '.join(f'{v:.4f}' for v in table5[-1]) + f' | Cbar={rows_Cbar[-1]:.4f}')

    idx_mean_probe = int(round(MEAN_CRIT_PROBE_H * SEC_PER_HOUR / Q34_FILE_DT_S))
    Cbar_at_36h = float(res['Cbar'][idx_mean_probe])
    echo(f'[Q3] 判据口径对照：Cbar({MEAN_CRIT_PROBE_H:.0f} h)={Cbar_at_36h:.4f} 已低于阈值 → '
         f'误用均值判据会得约 {MEAN_CRIT_PROBE_H:.0f} h，'
         f'比逐点判据少 {(1 - MEAN_CRIT_PROBE_H / t_dry_h) * PCT:.0f}%')

    radii_file = file_radii_cm(R0_M * CM_PER_M)
    Cf = field_at_radii(res, radii_file, 'C')
    header = [HEADER_CORNER] + radii_file
    times_file = [int(v) for v in res['t'][1:]]
    write_result_xlsx(os.path.join(ROOT, 'result3.xlsx'),
                      [('Sheet1', header, times_file, Cf[1:])])
    echo(f'[Q3] result3.xlsx 已写：1 表 × {len(times_file)} 行 × {len(radii_file)} 列 '
         f'(末行 t={times_file[-1]} s >= t_dry)')

    # ---- 网格无关性（建模报告 §7.3 的 N 集合）与能力清单补充 N 集合 ----
    grid_N = list(N_CONV)
    grid_t = []
    for n in grid_N:
        t0 = time.time()
        th, _ = t_dry_hours(props, air, radius, n, DT_Q3)
        grid_t.append(th)
        echo(f'[Q3] 网格自校验 N={n:3d}: t_dry={th:.6f} h ({time.time() - t0:.1f}s)')
    rel_diff = [None] + [abs(grid_t[i] - grid_t[i - 1]) / grid_t[i] for i in range(1, len(grid_t))]
    echo('[Q3] 相邻相对差 ' + ', '.join('—' if d is None else f'{d:.3e}' for d in rel_diff))

    grid_N_alt = list(GRID_N_ALT)
    grid_t_alt = []
    for n in grid_N_alt:
        th, _ = t_dry_hours(props, air, radius, n, DT_Q3)
        grid_t_alt.append(th)
        echo(f'[Q3] 补充网格 N={n:3d}: t_dry={th:.6f} h')
    rel_alt = [None] + [abs(grid_t_alt[i] - grid_t_alt[i - 1]) / grid_t_alt[i]
                        for i in range(1, len(grid_t_alt))]
    echo('[Q3] 补充集相邻相对差 ' + ', '.join('—' if d is None else f'{d:.3e}' for d in rel_alt))

    # ---- 时间步无关性（N 固定为 N_ABL） ----
    dt_list, dt_t = list(DT_CONV), []
    for dtv in dt_list:
        th, _ = t_dry_hours(props, air, radius, N_ABL, dtv)
        dt_t.append(th)
        echo(f'[Q3] 时间步自校验 N={N_ABL}, dt={dtv} s: t_dry={th:.6f} h')

    payload = {
        'meta': {'problem': 3, 'props': props.name, 'n_cells': N_PROD, 'dt_s': DT_Q3,
                 'decimals': DECIMALS, 'criterion': f'max_r C(r,t) < {C_TH} (pointwise)',
                 'runtime_s': round(time.time() - t_start, 2)},
        't_dry_h': t_dry_h, 't_dry_s': t_dry_s, 't_dry_d': t_dry_h / HOURS_PER_DAY,
        't_dry_reported': round(t_dry_h, DECIMALS),   # 与 RESULTS.md 引用的四位小数同口径
        'C_max_at_end': float(np.max(res['C'][-1])),
        'C_max_prev_step': float(res['Cmax'][IDX_PREV_STEP]),
        'max_pos_grad': res['max_pos_grad'],
        'table_radii_cm': table_radii,
        'table5_times_h': rows_t_h,
        'table5_C': to_list(np.asarray(rows_C), DECIMALS),
        'table5_T': to_list(np.asarray(rows_T), DECIMALS),
        'table5_Cbar': to_list(np.asarray(rows_Cbar), DECIMALS),
        'Cbar_at_36h': Cbar_at_36h,
        'series_t_h': to_list(res['t'] / SEC_PER_HOUR, SERIES_DECIMALS),
        'series_Cmax': to_list(res['Cmax'], SERIES_DECIMALS),
        'series_Cbar': to_list(res['Cbar'], SERIES_DECIMALS),
        'series_Csurf': to_list(surface_value(res, 'C'), SERIES_DECIMALS),
        'series_Tcenter': to_list(res['T'][:, 0], 4),
        'series_Tsurf': to_list(res['T'][:, -1], 4),
        'profile_C_5r': to_list(C_tab_all, DECIMALS),
        'file_radii_cm': radii_file,
        'field_C_full': to_list(Cf[::JSON_STRIDE], DECIMALS),
        'field_C_full_stride': JSON_STRIDE,
        'grid_convergence': {'N': grid_N, 't_dry_h': grid_t, 'rel_diff': rel_diff},
        'grid_convergence_alt': {'N': grid_N_alt, 't_dry_h': grid_t_alt, 'rel_diff': rel_alt},
        'dt_convergence': {'n_cells': N_ABL, 'dt_s': dt_list, 't_dry_h': dt_t},
        'grid_rel_diff': rel_diff[grid_N.index(N_PROD)],
        'conservation_residual': res['cons_res'],
        'xlsx_rows': len(times_file),
        'baseline_ablation_n_cells': N_ABL,
        'baseline_ablation_t_dry_h': grid_t[grid_N.index(N_ABL)],
    }
    validate_capability(payload, res, t_dry_h)
    dump_json(os.path.join(FIG_DIR, 'problem_3_results.json'), payload)
    dump_json(os.path.join(FIG_DIR, 'q3_results.json'), payload)
    echo(f'[Q3] DONE t_dry={t_dry_h:.4f} h（{t_dry_h / HOURS_PER_DAY:.4f} d），'
         f'总耗时 {time.time() - t_start:.1f}s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

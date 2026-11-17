"""问题二：整个烘干过程的变物性热-质耦合场（制表窗口 3 h）。

附录 3 闭合全程统一（题面明示"相关经验公式统一采用附录 3 中的公式"），
不按预热/恒温阶段切换；R=R0；h 与 hm 沿用附录 2。
产出：表 3（温度）、表 4（水分浓度）、result2.xlsx（1 s 步长 × 0.1 cm 步长，全时窗）。
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


def validate_capability(payload, res, tab_T, tab_C, air):
    """P2-C1 / P2-C3 / P2-C4 的可伪证检查。"""
    props = props_appendix3()
    k_at_C0 = float(props.k(C0_KGKG, T0_C))
    k_lo, k_hi = K_A3_AT_C0_RANGE
    assert k_lo <= k_at_C0 <= k_hi, (
        f'P2-C1：k({C0_KGKG})={k_at_C0:.5f} 越界，湿基组合 C/(C+1) 可能被误写为 C')
    T_kelvin_used = T0_C + KELVIN_OFFSET
    assert T_kelvin_used > KELVIN_OFFSET, 'P2-C1：D 中温度须为开尔文'
    # 独立于 params.props_appendix3 再按 FACTS 的系数装一次 D，两路必须逐位一致
    a3 = FACTS['appendix3_params']
    D_from_facts = (a3['D_prefactor_mantissa'] * 10.0 ** a3['D_prefactor_exponent']
                    * np.exp(-a3['D_activation_over_C'] / C0_KGKG)
                    * np.exp(-a3['D_activation_over_T'] / T_kelvin_used))
    assert rel_err(float(props.D(C0_KGKG, T0_C)), D_from_facts) < EXACT_TOL, \
        'P2-C1：D 表达式与附录 3 逐字形式不符'
    assert air.n_points == N_A1_POINTS, \
        f'R8/P2-C4：附件 1 应读入 {N_A1_POINTS} 点，实为 {air.n_points}'
    assert rel_err(air.T_tail_mean, AIR_T_PLATEAU_C) < PLATEAU_REL_TOL, \
        f'P2-C4：末 {AIR_TAIL_POINTS} 点温度均值偏离 {AIR_T_PLATEAU_C}'
    assert rel_err(air.W_tail_mean, AIR_W_PLATEAU) < PLATEAU_W_REL_TOL, \
        f'P2-C4：末 {AIR_TAIL_POINTS} 点含湿均值偏离 {AIR_W_PLATEAU}'
    assert abs(float(air.T_of(AIR_PROBE_LATE_S)) - air.T_tail_mean) < EXACT_TOL, (
        f'R9/P2-C4：t={AIR_PROBE_LATE_S:.0f} s 须返回平台均值'
        f'（线性外推会超过 {AIR_T_LINEXTRAP_MAX_C} °C）')
    assert float(air.T_of(AIR_PROBE_LATE_S)) < AIR_T_LINEXTRAP_MAX_C, \
        'R9/P2-C4：延拓值异常，疑似线性外推'
    assert payload['xlsx_rows'] == Q2_HORIZON_S, (
        f'P2-C3：result2.xlsx 应有 {Q2_HORIZON_S} 数据行')
    n_tab_rows = int(round(Q2_HORIZON_S / SEC_PER_HOUR / Q2_TABLE_DT_H))
    assert tab_T.shape[1] == len(Q1_TABLE_RADII_CM) and tab_C.shape[0] == n_tab_rows, \
        'P2-C2：表 3/4 形状错'
    assert np.all(np.diff(tab_T[-1]) > 0.0), '温度剖面须径向递增'
    assert np.all(np.diff(tab_C[-1]) < 0.0), '含水率剖面须径向递减'
    assert res['cons_res'] < CONS_RES_TOL, '干质量守恒残差超限'
    echo('[capability] P2-C1/P2-C2/P2-C3/P2-C4 运行时检查全部通过')


def main() -> int:
    set_all_seeds()
    t_start = time.time()
    props = props_appendix3()
    air = AirState()
    radius = FixedRadius()
    k_at_C0 = float(props.k(C0_KGKG, T0_C))
    echo(f'[Q2] 附录 3 变物性：k(2.55)={k_at_C0:.5f}, rho(2.55)={float(props.rho(C0_KGKG, T0_C)):.1f}, '
         f'cp(2.55)={float(props.cp(C0_KGKG, T0_C)):.1f}, D(2.55,28C)={float(props.D(C0_KGKG, T0_C)):.4e}')
    echo(f'[Q2] 附件 1：{air.n_points} 点，末 {AIR_TAIL_POINTS} 点均值 '
         f'T={air.T_tail_mean:.4f} °C / '
         f'W={air.W_tail_mean:.6f} kg/kg（延拓自 t={air.tail_t_start:.0f} s）')
    t_end, t_late = float(Q2_HORIZON_S), AIR_PROBE_LATE_S
    echo(f'[Q2] 延拓核验：T_air({t_end:.0f})={float(air.T_of(t_end)):.4f}, '
         f'W_air({t_end:.0f})={float(air.W_of(t_end)):.5f}; '
         f'T_air({t_late:.0f})={float(air.T_of(t_late)):.6f}, '
         f'W_air({t_late:.0f})={float(air.W_of(t_late)):.6f}')

    res = solve_transport(props, air, radius, N_PROD, DT_Q2, float(Q2_HORIZON_S),
                          float(Q2_FILE_DT_S))
    echo(f'[Q2] N={N_PROD}, dt={DT_Q2} s, 步数={res["n_steps"]}, '
         f'守恒残差={res["cons_res"]:.3e}, 耗时={time.time() - t_start:.1f}s')

    table_times_s = [int(round(k * Q2_TABLE_DT_H * SEC_PER_HOUR))
                     for k in range(1, int(Q2_HORIZON_S / (Q2_TABLE_DT_H * SEC_PER_HOUR)) + 1)]
    idx_tab = [int(round(t / Q2_FILE_DT_S)) for t in table_times_s]
    T_all = field_at_radii(res, Q1_TABLE_RADII_CM, 'T')
    C_all = field_at_radii(res, Q1_TABLE_RADII_CM, 'C')
    tab_T, tab_C = T_all[idx_tab], C_all[idx_tab]
    echo('[Q2] 表 3 末行 T(3 h) = ' + ' '.join(f'{v:.4f}' for v in tab_T[-1]))
    echo('[Q2] 表 4 末行 C(3 h) = ' + ' '.join(f'{v:.4f}' for v in tab_C[-1]))
    echo(f'[Q2] 3 h 内外温差 {tab_T[-1][-1] - tab_T[-1][0]:.4f} °C，'
         f'中心距平台温度 {air.T_tail_mean - tab_T[-1][0]:.4f} °C，'
         f'中心含水率完成全程降幅 {(C0_KGKG - tab_C[-1][0]) / (C0_KGKG - C_TH) * 100:.1f}%')

    radii_file = file_radii_cm(R0_M * CM_PER_M)
    Tf = field_at_radii(res, radii_file, 'T')
    Cf = field_at_radii(res, radii_file, 'C')
    header = [HEADER_CORNER] + radii_file
    times_file = [int(v) for v in res['t'][1:]]
    out_path = os.path.join(ROOT, 'result2.xlsx')
    write_result_xlsx(out_path, [(SHEET_T, header, times_file, Tf[1:]),
                                 (SHEET_C, header, times_file, Cf[1:])])
    echo(f'[Q2] result2.xlsx 已写：2 表 × {len(times_file)} 行 × {len(radii_file)} 列')

    # JSON 里给 60 s 粒度的曲线（Q34_FILE_DT_S 即 60 s），xlsx 保留 1 s 全量
    stride = int(round(Q34_FILE_DT_S / Q2_FILE_DT_S))
    payload = {
        'meta': {'problem': 2, 'props': props.name, 'n_cells': N_PROD, 'dt_s': DT_Q2,
                 'horizon_s': Q2_HORIZON_S, 'decimals': DECIMALS,
                 f'k_at_{C0_KGKG}': k_at_C0, 'D_at_C0_T0': float(props.D(C0_KGKG, T0_C)),
                 'runtime_s': round(time.time() - t_start, 2)},
        'field_shape': list(res['T'].shape),
        'file_radii_cm': radii_file,
        'table_times_h': [t / SEC_PER_HOUR for t in table_times_s],
        'table_radii_cm': list(Q1_TABLE_RADII_CM),
        'table3_T': to_list(tab_T, DECIMALS),
        'table4_C': to_list(tab_C, DECIMALS),
        'profile_t_s': to_list(res['t'][::stride]),
        'profile_T_5r': to_list(T_all[::stride], DECIMALS),
        'profile_C_5r': to_list(C_all[::stride], DECIMALS),
        'field_T_full': to_list(Tf[::stride], DECIMALS),
        'field_C_full': to_list(Cf[::stride], DECIMALS),
        'Cbar': to_list(res['Cbar'][::stride], SERIES_DECIMALS),
        'air': {'n_points': air.n_points,
                'T_tail_mean': air.T_tail_mean, 'W_tail_mean': air.W_tail_mean,
                'tail_t_start_s': air.tail_t_start,
                'T_at_20000': float(air.T_of(t_late)), 'W_at_20000': float(air.W_of(t_late)),
                'T_at_10800': float(air.T_of(t_end)), 'W_at_10800': float(air.W_of(t_end)),
                't_s': to_list(air.t), 'T_c': to_list(air.T, DECIMALS),
                'W_kgkg': to_list(air.W, DECIMALS + 1)},
        'conservation_residual': res['cons_res'],
        'xlsx_rows': len(times_file),
    }
    validate_capability(payload, res, tab_T, tab_C, air)
    dump_json(os.path.join(FIG_DIR, 'problem_2_results.json'), payload)
    dump_json(os.path.join(FIG_DIR, 'q2_results.json'), payload)
    echo(f'[Q2] DONE 总耗时 {time.time() - t_start:.1f}s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

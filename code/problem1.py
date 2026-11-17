"""问题一：预热平衡阶段（0-1800 s）的常物性热-质耦合场。

附录 2 闭合（rho/cp/k 取常数，D 只含浓度活化项、不含温度项），R=R0。
产出：表 1（温度）、表 2（水分浓度）、result1.xlsx（1 s × 0.1 cm）、
以及 Bessel 特征函数展开 + Duhamel 叠加（减法形式）的独立解析校验。
"""
from __future__ import annotations

import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import numpy as np
from scipy.optimize import brentq
from scipy.special import j0, j1, jn_zeros

from params import *  # noqa: F401,F403  口径常数（显式声明来源，勿依赖 utils 转出）
from utils import *


def bessel_eigenvalues(n_modes, Bi):
    """特征方程 zeta*J1(zeta) = Bi*J0(zeta)：在 J1 相邻零点区间内逐个 Brent 求根。

    统一初值迭代会让多个模态收敛到同一根（建模报告的失败模式清单），故必须逐区间隔离。
    """
    def eig_eq(z):
        return z * j1(z) - Bi * j0(z)

    zeros = np.concatenate(([0.0], jn_zeros(1, n_modes + 2)))
    roots = []
    for i in range(n_modes):
        lo, hi = zeros[i], zeros[i + 1]
        a, b = lo + GRAD_POS_TOL, hi - GRAD_POS_TOL      # 端点内缩，避开 J1 零点处的符号退化
        if eig_eq(a) * eig_eq(b) > 0.0:
            continue
        roots.append(brentq(eig_eq, a, b, xtol=BRENT_XTOL, rtol=BRENT_RTOL,
                            maxiter=BRENT_MAXITER))
    return np.asarray(roots[:n_modes])


def analytic_temperature(air, radii_m, t_out, n_modes):
    """减法形式 theta = T - T_air(t) 的 Duhamel 叠加解。

    theta 满足齐次 Robin，级数系数按代数速率衰减；常数部分由前导项 T_air(t) 承担。
    附件 1 在 60 s 间隔上分段线性，故每段内 dT_air/dt 为常数，b_n 可精确积分。
    """
    Bi = BIOT_A2
    zeta = bessel_eigenvalues(n_modes, Bi)
    lam = ALPHA_A2 * zeta ** 2 / R0_M ** 2
    A_n = 2.0 * j1(zeta) / (zeta * (j0(zeta) ** 2 + j1(zeta) ** 2))
    shapes = j0(np.outer(radii_m / R0_M, zeta))          # (n_r, n_modes)

    n_keep = int(np.searchsorted(air.t, t_out[-1], side='left')) + 1
    tg = air.t[:n_keep + 1] if n_keep < air.t.size else air.t
    Tg = air.T[:tg.size]
    assert tg[-1] >= t_out[-1] - 1e-9, '空气分段窗口未覆盖输出末时刻（会退化为区间外推）'
    b = A_n * (T0_C - float(air.T_of(0.0)))
    seg_t = [0.0]
    seg_b = [b.copy()]
    for jseg in range(tg.size - 1):
        delta = tg[jseg + 1] - tg[jseg]
        slope = (Tg[jseg + 1] - Tg[jseg]) / delta
        decay = np.exp(-lam * delta)
        b = b * decay - A_n * slope * (1.0 - decay) / lam
        seg_t.append(tg[jseg + 1])
        seg_b.append(b.copy())
    seg_t = np.asarray(seg_t)
    seg_B = np.asarray(seg_b)

    # 减法形式的重构 T = T_air(t) + sum_n b_n J0(zeta_n r/R)：显式逐模求和与
    # 下面的矩阵乘法互校（取表面半径、末时刻），防 shapes 装配或模态对齐写错
    T_air_end = float(air.T_of(seg_t[-1]))
    r_probe = float(radii_m[-1]) / R0_M
    explicit = T_air_end + sum(b * float(j0(z * r_probe))
                               for b, z in zip(seg_B[-1], zeta))
    assert abs(explicit - (T_air_end + shapes[-1] @ seg_B[-1])) < 1e-12, \
        'P1-C4：逐模求和与向量化装配不一致'

    out = np.empty((t_out.size, radii_m.size))
    for i, t_q in enumerate(t_out):
        jseg = int(np.searchsorted(seg_t, t_q, side='right')) - 1
        jseg = min(max(jseg, 0), seg_t.size - 2) if seg_t.size > 1 else 0
        tau = t_q - seg_t[jseg]
        if tau <= 0.0:
            b_q = seg_B[jseg]
        else:
            delta = seg_t[jseg + 1] - seg_t[jseg]
            slope = (Tg[jseg + 1] - Tg[jseg]) / delta
            decay = np.exp(-lam * tau)
            b_q = seg_B[jseg] * decay - A_n * slope * (1.0 - decay) / lam
        out[i] = float(air.T_of(t_q)) + shapes @ b_q
    return out, zeta, lam, A_n


def validate_capability(payload, res, tab_T, tab_C, analytic_dev, zeta_dev):
    """P1-C1 / P1-C4 的可伪证检查（falsifiable_check 的运行时实现）。"""
    assert res['T'].ndim == 2 and res['C'].ndim == 2, 'P1-C1：必须输出二维场'
    n_r_file = len(payload['file_radii_cm'])
    assert payload['field_shape'][1] >= N_FILE_RADII and n_r_file >= N_FILE_RADII, \
        f'P1-C1：径向点数不足 {N_FILE_RADII}'
    assert abs(res['s'][-1] - 1.0) < BRENT_RTOL, 'P1-C1：末节点须精确对应 r=R'
    assert tab_T.shape == (len(Q1_TABLE_TIMES_S), len(Q1_TABLE_RADII_CM)), 'P1-C2：表 1 形状错'
    assert tab_C.shape == (len(Q1_TABLE_TIMES_S), len(Q1_TABLE_RADII_CM)), 'P1-C2：表 2 形状错'
    assert np.all(np.diff(tab_T[-1]) > 0.0), 'R3：末时刻温度剖面须径向递增（表面最热）'
    assert abs(tab_C[-1][0] - C0_KGKG) < C_CENTER_TOL, f'R2：末时刻中心含水率应仍为 {C0_KGKG}'
    assert analytic_dev < ANALYTIC_DEV_TOL, \
        f'R4/P1-C4：解析-数值偏差 {analytic_dev:.3e} 未达 {ANALYTIC_DEV_TOL:g}'
    assert zeta_dev < ZETA_DEV_TOL, f'R5/P1-C4：特征根偏差 {zeta_dev:.3e} 未达 {ZETA_DEV_TOL:g}'
    assert res['cons_res'] < CONS_RES_TOL, 'P1-C4：干质量守恒残差超限'
    echo('[capability] P1-C1/P1-C2/P1-C3/P1-C4 运行时检查全部通过')


def main() -> int:
    set_all_seeds()
    t_start = time.time()
    props = props_appendix2()
    air = AirState()
    radius = FixedRadius()
    echo(f'[Q1] 附录 2 常物性 alpha={ALPHA_A2:.9e} m^2/s, Bi={BIOT_A2:.6f}, '
         f'Fo(1800)={ALPHA_A2 * Q1_HORIZON_S / R0_M ** 2:.6f}')

    res = solve_transport(props, air, radius, N_PROD, DT_Q1, float(Q1_HORIZON_S),
                          float(Q1_FILE_DT_S))
    echo(f'[Q1] 生产网格 N={N_PROD}, dt={DT_Q1} s, 步数={res["n_steps"]}, '
         f'守恒残差={res["cons_res"]:.3e}, 耗时={time.time() - t_start:.1f}s')

    idx_tab = [int(round(t / Q1_FILE_DT_S)) for t in Q1_TABLE_TIMES_S]
    T_all = field_at_radii(res, Q1_TABLE_RADII_CM, 'T')
    C_all = field_at_radii(res, Q1_TABLE_RADII_CM, 'C')
    tab_T, tab_C = T_all[idx_tab], C_all[idx_tab]
    echo('[Q1] 表 1 末行 T(1800 s) = ' + ' '.join(f'{v:.4f}' for v in tab_T[-1]))
    echo('[Q1] 表 2 末行 C(1800 s) = ' + ' '.join(f'{v:.4f}' for v in tab_C[-1]))

    radii_file = file_radii_cm(R0_M * CM_PER_M)
    Tf = field_at_radii(res, radii_file, 'T')
    Cf = field_at_radii(res, radii_file, 'C')
    header = [HEADER_CORNER] + radii_file
    times_file = [int(v) for v in res['t'][1:]]
    out_path = os.path.join(ROOT, 'result1.xlsx')
    write_result_xlsx(out_path, [(SHEET_T, header, times_file, Tf[1:]),
                                 (SHEET_C, header, times_file, Cf[1:])])
    echo(f'[Q1] result1.xlsx 已写：2 表 × {len(times_file)} 行 × {len(radii_file)} 列')

    # ---- 独立解析校验：Duhamel + Bessel（减法形式），不同求解路径 ----
    radii_m = np.asarray(radii_file, dtype=float) * (1.0 / CM_PER_M)
    conv = {}
    for M in N_MODES_CONV:
        Ta, zeta, lam, A_n = analytic_temperature(air, radii_m, res['t'], M)
        conv[M] = float(np.max(np.abs(Ta - Tf)))
    analytic_dev = conv[N_MODES_CONV[-1]]
    zeta_dev = float(np.max(np.abs(zeta[:len(ZETA_REF)] - np.asarray(ZETA_REF))))
    echo('[Q1] 解析校验（减法形式）模态收敛 ' +
         ', '.join(f'M={M}:{v:.3e}' for M, v in conv.items()) + ' °C（生产网格）')

    # 细网格复核：生产网格上的 5.2e-4 是离散误差地板，细化网格可暴露级数是否真收敛
    n_fine, dt_fine = N_PROD * 4, DT_Q1 / 5.0
    res_fine = solve_transport(props, air, radius, n_fine, dt_fine,
                               float(Q1_HORIZON_S), float(Q1_FILE_DT_S))
    Tf_fine = field_at_radii(res_fine, radii_file, 'T')
    conv_fine = {}
    for M in N_MODES_CONV:
        Ta_f, _, _, _ = analytic_temperature(air, radii_m, res_fine['t'], M)
        conv_fine[M] = float(np.max(np.abs(Ta_f - Tf_fine)))
    echo(f'[Q1] 细网格 N={n_fine}, dt={dt_fine} s 模态收敛 ' +
         ', '.join(f'M={M}:{v:.3e}' for M, v in conv_fine.items()) + ' °C')
    echo(f'[Q1] zeta_1..{len(ZETA_REF)} = ' +
         ' '.join(f'{z:.5f}' for z in zeta[:len(ZETA_REF)]) +
         f'，与参照偏差 {zeta_dev:.2e}；lambda_1={lam[0]:.6e} 1/s')

    payload = {
        'meta': {'problem': 1, 'props': props.name, 'n_cells': N_PROD, 'dt_s': DT_Q1,
                 'horizon_s': Q1_HORIZON_S, 'decimals': DECIMALS,
                 'alpha': ALPHA_A2, 'Bi': BIOT_A2,
                 'Fo_1800': ALPHA_A2 * Q1_HORIZON_S / R0_M ** 2,
                 'runtime_s': round(time.time() - t_start, 2)},
        'field_shape': list(res['T'].shape),
        'file_radii_cm': radii_file,
        'table_times_s': list(Q1_TABLE_TIMES_S),
        'table_radii_cm': list(Q1_TABLE_RADII_CM),
        'table1_T': to_list(tab_T, DECIMALS),
        'table2_C': to_list(tab_C, DECIMALS),
        'air_T_at_table_times': to_list([float(air.T_of(t)) for t in Q1_TABLE_TIMES_S], SERIES_DECIMALS),
        'air_W_at_table_times': to_list([float(air.W_of(t)) for t in Q1_TABLE_TIMES_S], SERIES_DECIMALS),
        'profile_t_s': to_list(res['t'][1:]),
        'profile_T_5r': to_list(T_all[1:], DECIMALS),
        'profile_C_5r': to_list(C_all[1:], DECIMALS),
        'field_T_full': to_list(Tf[1:], DECIMALS),
        'field_C_full': to_list(Cf[1:], DECIMALS),
        'Cbar': to_list(res['Cbar'], SERIES_DECIMALS),
        'analytic': {'zeta': to_list(zeta[:len(ZETA_REF)], SERIES_DECIMALS),
                     'zeta_ref': list(ZETA_REF), 'zeta_dev': zeta_dev,
                     'lambda_1': float(lam[0]),
                     'mode_convergence': {str(k): v for k, v in conv.items()},
                     'mode_convergence_fine': {str(k): v for k, v in conv_fine.items()},
                     'fine_grid': {'n_cells': n_fine, 'dt_s': dt_fine,
                                   'analytic_max_dev': conv_fine[N_MODES_CONV[-1]]},
                     'analytic_max_dev': analytic_dev,
                     'T_air_0': float(air.T_of(0.0))},
        'conservation_residual': res['cons_res'],
        'xlsx_rows': len(times_file),
    }
    validate_capability(payload, res, tab_T, tab_C, analytic_dev, zeta_dev)
    dump_json(os.path.join(FIG_DIR, 'problem_1_results.json'), payload)
    dump_json(os.path.join(FIG_DIR, 'q1_results.json'), payload)
    echo(f'[Q1] DONE 总耗时 {time.time() - t_start:.1f}s，解析偏差 {analytic_dev:.3e} °C')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

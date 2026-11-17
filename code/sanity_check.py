"""结构性验证（建模报告的退化检验清单）：四项退化/守恒检验 + 网格无关性汇总。

这些检验不产出题目答案，只用来证明离散化本身没写错：
零输入退化查 Robin 符号，常空气解析退化查首模衰减率，h_m=0 查内部通量反对称，
收缩关闭退化查"问题四代码能复现问题三"这一跨问题闭环。
"""
from __future__ import annotations

import json
import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import numpy as np
from scipy.special import j0, j1, jn_zeros
from scipy.optimize import brentq

from params import *  # noqa: F401,F403  口径常数（显式声明来源，勿依赖 utils 转出）
from utils import *


def test_zero_input():
    """T_air=T0、W_air=C0：解必须恒为初值（检验 Robin 项符号与源项装配）。"""
    air = ConstantAir(T0_C, C0_KGKG)
    res = solve_transport(props_appendix3(), air, FixedRadius(), N_ABL, DT_Q3,
                          DEGEN_HORIZON_S, float(Q34_FILE_DT_S))
    dT = float(np.max(np.abs(res['T'] - T0_C)))
    dC = float(np.max(np.abs(res['C'] - C0_KGKG)))
    ok = dT < DEGEN_TOL and dC < DEGEN_TOL
    echo(f'[sanity] 零输入退化：max|T-{T0_C}|={dT:.3e}, max|C-{C0_KGKG}|={dC:.3e} '
         f'(<{DEGEN_TOL:g}) → {"PASS" if ok else "FAIL"}')
    return {'name': 'zero_input_degeneracy', 'max_dT': dT, 'max_dC': dC,
            'tol': DEGEN_TOL, 'pass': bool(ok)}


def test_constant_air_lambda1():
    """T_air=CONST_AIR_T + 附录 2 常物性：长时衰减率须等于 lambda_1=alpha*zeta_1^2/R0^2。"""
    zeta1 = brentq(lambda z: z * j1(z) - BIOT_A2 * j0(z),
                   GRAD_POS_TOL, jn_zeros(1, 1)[0] - GRAD_POS_TOL,
                   xtol=BRENT_XTOL, rtol=BRENT_RTOL)
    lam1 = ALPHA_A2 * zeta1 ** 2 / R0_M ** 2
    air = ConstantAir(CONST_AIR_T, C0_KGKG)
    res = solve_transport(props_appendix2(), air, FixedRadius(), N_ABL, DT_Q2,
                          LAMBDA1_HORIZON_S, float(Q34_FILE_DT_S))
    t = res['t']
    dev = np.abs(res['T'][:, 0] - CONST_AIR_T)
    w_lo, w_hi = LAMBDA1_FIT_WINDOW_S            # 高阶模在窗左端前已衰减到 1e-6 以下
    win = (t >= w_lo) & (t <= w_hi)
    slope = np.polyfit(t[win], np.log(dev[win]), 1)[0]
    lam_fit = -float(slope)
    rel = abs(lam_fit - lam1) / lam1
    ok = rel < LAMBDA_FIT_REL_TOL
    echo(f'[sanity] 常空气解析退化：拟合衰减率 {lam_fit:.6e} vs lambda_1={lam1:.6e} 1/s，'
         f'相对偏差 {rel:.3e} (<{LAMBDA_FIT_REL_TOL:g}) → {"PASS" if ok else "FAIL"}')
    return {'name': 'constant_air_lambda1', 'zeta_1': float(zeta1), 'lambda_1': float(lam1),
            'lambda_fit': lam_fit, 'rel_dev': rel, 'fit_window_s': [w_lo, w_hi],
            'tol': LAMBDA_FIT_REL_TOL, 'pass': bool(ok)}


def test_pure_diffusion_conservation():
    """h_m=0：Cbar 必须严格守恒（检验内部通量装配的反对称性）。"""
    res = solve_transport(props_appendix3(), AirState(), FixedRadius(), N_ABL, DT_Q3,
                          SEC_PER_HOUR, float(Q34_FILE_DT_S), hm=0.0)
    Cbar = res['Cbar']
    drift = float(np.max(np.abs(Cbar - Cbar[0])) / abs(Cbar[0]))
    ok = drift < CONSERVE_DRIFT_TOL
    echo(f'[sanity] 纯扩散守恒（h_m=0）：Cbar 相对漂移 {drift:.3e} '
         f'(<{CONSERVE_DRIFT_TOL:g}) → {"PASS" if ok else "FAIL"}')
    return {'name': 'pure_diffusion_conservation', 'Cbar_0': float(Cbar[0]),
            'Cbar_end': float(Cbar[-1]), 'rel_drift': drift,
            'tol': CONSERVE_DRIFT_TOL, 'pass': bool(ok)}


def test_shrinkage_off_reproduces_q3():
    """R(t)=R0 + 附录 3：问题四的求解路径必须复现问题三的 t_dry（相对差<1e-6）。"""
    with open(os.path.join(FIG_DIR, 'problem_3_results.json'), encoding='utf-8') as fh:
        q3 = json.load(fh)
    grid = q3['grid_convergence']
    t_q3 = float(grid['t_dry_h'][grid['N'].index(N_ABL)])
    t_off, _ = t_dry_hours(props_appendix3(), AirState(), FixedRadius(), N_ABL, DT_Q4)
    rel = abs(t_off - t_q3) / t_q3
    ok = rel < CROSS_PROBLEM_TOL
    echo(f'[sanity] 收缩关闭退化：问题四路径 {t_off:.6f} h vs 问题三 {t_q3:.6f} h，'
         f'相对差 {rel:.3e} (<{CROSS_PROBLEM_TOL:g}) → {"PASS" if ok else "FAIL"}')
    return {'name': 'shrinkage_off_reproduces_q3', 't_dry_shrinkage_off_h': t_off,
            't_dry_q3_h': t_q3, 'rel_dev': rel, 'n_cells': N_ABL,
            'tol': CROSS_PROBLEM_TOL, 'pass': bool(ok)}


def grid_summary():
    """网格无关性汇总：直接引用问题三/四已落盘的真实重算结果，不重复计算。"""
    out = {}
    for prob, fname in ((3, 'problem_3_results.json'), (4, 'problem_4_results.json')):
        with open(os.path.join(FIG_DIR, fname), encoding='utf-8') as fh:
            grid = json.load(fh)['grid_convergence']
        diffs = [d for d in grid['rel_diff'] if d is not None]
        mono = all(diffs[i] > diffs[i + 1] for i in range(len(diffs) - 1))
        ok = all(d < GRID_REL_TOL for d in diffs[IDX_PREV_STEP:])
        out[f'problem{prob}'] = {'N': grid['N'], 't_dry_h': grid['t_dry_h'],
                                 'rel_diff': diffs, 'monotone_decreasing': bool(mono),
                                 'last_within_tol': bool(ok)}
        echo(f'[sanity] 网格无关性（问题{prob}，引用已落盘重算）N={grid["N"]}，'
             f'相邻相对差 ' + ', '.join(f'{d:.3e}' for d in diffs) +
             f' → 单调下降 {mono}，末组达标 {ok}')
    return out


def main() -> int:
    set_all_seeds()
    t_start = time.time()
    checks = [test_zero_input(), test_constant_air_lambda1(),
              test_pure_diffusion_conservation(), test_shrinkage_off_reproduces_q3()]
    grids = grid_summary()
    n_pass = sum(1 for c in checks if c['pass'])
    payload = {
        'meta': {'section': 'structural_verification', 'source': 'MODELING_REPORT.md',
                 'n_checks': len(checks), 'n_pass': n_pass,
                 'runtime_s': round(time.time() - t_start, 2)},
        'checks': checks,
        'grid_independence': grids,
        'all_pass': bool(n_pass == len(checks)),
    }
    assert payload['all_pass'], f'结构性验证未全通过：{n_pass}/{len(checks)}'
    assert all(g['last_within_tol'] for g in grids.values()), '网格无关性判据未达标'
    dump_json(os.path.join(FIG_DIR, 'sanity_check_results.json'), payload)
    echo(f'[sanity] DONE {n_pass}/{len(checks)} 项结构性验证通过，'
         f'耗时 {time.time() - t_start:.1f}s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

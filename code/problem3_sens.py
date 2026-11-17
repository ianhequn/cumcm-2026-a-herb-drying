"""灵敏度与鲁棒性（建模报告的灵敏度表）：每行一次完整真实重算到阈值。

基准为 N=N_ABL、dt=DT_Q3 下的问题三口径；R(t) 缩放两行以问题四口径为基准
（缩放半径必须配附录 4 物性，否则不是同一模型的扰动）。
禁止由基准结果解析放缩——每行都调用 t_dry_hours 重跑一次。
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


def run_row(label, props, air, radius, base_h, hm=HM, h=H_CONV, driver=''):
    """跑一行灵敏度：返回 (label, t_dry_h, 相对变化, 基准标签)。"""
    t0 = time.time()
    th, _ = t_dry_hours(props, air, radius, N_ABL, DT_Q3, hm=hm, h=h)
    rel = (th - base_h) / base_h * PCT
    echo(f'[Q3-sens] {label:<28s} t_dry={th:.4f} h  ({rel:+.3f}%)  [{time.time() - t0:.1f}s]')
    return {'label': label, 'driver': driver, 't_dry_h': th, 'rel_pct': rel}


def validate_capability(payload):
    """P3-C4 / R14 的可伪证检查：真实重算 + h_m 双向反号 + 排序。"""
    rows = payload['rows']
    assert len(rows) >= SENS_MIN_ROWS, f'P3-C4：灵敏度行数 {len(rows)} < {SENS_MIN_ROWS}'
    vals = [r['t_dry_h'] for r in rows]
    assert len(set(round(v, SERIES_DECIMALS) for v in vals)) == len(vals), (
        'P3-C4：出现重复 t_dry，疑似解析放缩而非真实重算')
    hm_up = next(r for r in rows if r['driver'] == 'hm_up')
    hm_dn = next(r for r in rows if r['driver'] == 'hm_down')
    assert hm_up['rel_pct'] * hm_dn['rel_pct'] < 0.0, (
        f'R14：h_m 双向扰动未反号（{hm_up["rel_pct"]:+.3f}% / {hm_dn["rel_pct"]:+.3f}%）')
    assert hm_up['rel_pct'] < 0.0, 'R14：提高 h_m 必须缩短烘干时长'
    T_up = next(r for r in rows if r['driver'] == 'T_air_up')
    T_dn = next(r for r in rows if r['driver'] == 'T_air_down')
    assert T_up['rel_pct'] < 0.0 < T_dn['rel_pct'], 'R14：T_air 升高须缩短时长'
    W_up = next(r for r in rows if r['driver'] == 'W_air_up')
    assert W_up['rel_pct'] > 0.0, 'R14：空气含湿升高须延长时长（推动力减小）'
    R_up = next(r for r in rows if r['driver'] == 'R_scale_up')
    R_dn = next(r for r in rows if r['driver'] == 'R_scale_down')
    assert R_up['rel_pct'] > 0.0 > R_dn['rel_pct'], 'R14：半径放大须延长时长'
    order = payload['ranking']
    assert order[0] == SENS_TOP_DRIVER and order[-1] == SENS_BOTTOM_DRIVER, \
        f'灵敏度排序异常：{order}'
    assert payload['elasticity_R'] > 1.0, (
        f'R14：dlnt/dlnR={payload["elasticity_R"]:.3f} 应接近扩散标度 '
        f'{DIFFUSION_SCALING_EXPONENT:g}')
    echo('[capability] P3-C4/R14 灵敏度运行时检查全部通过')


def main() -> int:
    set_all_seeds()
    t_start = time.time()
    p3, p4 = props_appendix3(), props_appendix4()
    air0 = AirState()
    fixed = FixedRadius()

    base_h, _ = t_dry_hours(p3, air0, fixed, N_ABL, DT_Q3)
    echo(f'[Q3-sens] 基准（问题三口径 N={N_ABL}, dt={DT_Q3} s）t_dry={base_h:.4f} h')

    rows = [{'label': '基准（问题三）', 'driver': 'baseline',
             't_dry_h': base_h, 'rel_pct': 0.0}]
    rows.append(run_row(f'T_air+{SENS_DT_AIR} °C', p3, AirState(dT=+SENS_DT_AIR), fixed,
                        base_h, driver='T_air_up'))
    rows.append(run_row(f'T_air-{SENS_DT_AIR} °C', p3, AirState(dT=-SENS_DT_AIR), fixed,
                        base_h, driver='T_air_down'))
    rows.append(run_row(f'W_air+{SENS_DW_AIR} kg/kg', p3, AirState(dW=+SENS_DW_AIR), fixed,
                        base_h, driver='W_air_up'))
    rows.append(run_row(f'W_air-{SENS_DW_AIR} kg/kg', p3, AirState(dW=-SENS_DW_AIR), fixed,
                        base_h, driver='W_air_down'))
    rows.append(run_row(f'h_m x{SENS_SCALE_UP}', p3, air0, fixed, base_h,
                        hm=HM * SENS_SCALE_UP, driver='hm_up'))
    rows.append(run_row(f'h_m x{SENS_SCALE_DOWN}', p3, air0, fixed, base_h,
                        hm=HM * SENS_SCALE_DOWN, driver='hm_down'))
    rows.append(run_row(f'h x{SENS_SCALE_UP}', p3, air0, fixed, base_h,
                        h=H_CONV * SENS_SCALE_UP, driver='h_up'))
    rows.append(run_row(f'h x{SENS_SCALE_DOWN}', p3, air0, fixed, base_h,
                        h=H_CONV * SENS_SCALE_DOWN, driver='h_down'))

    base4_h, _ = t_dry_hours(p4, air0, RadiusHistory(), N_ABL, DT_Q4)
    echo(f'[Q3-sens] 半径缩放行的基准（问题四口径）t_dry={base4_h:.4f} h')
    rows.append(run_row(f'R(t) x{SENS_R_SCALE_UP}（基准问题四）', p4, air0,
                        RadiusHistory(scale=SENS_R_SCALE_UP), base4_h, driver='R_scale_up'))
    rows.append(run_row(f'R(t) x{SENS_R_SCALE_DOWN}（基准问题四）', p4, air0,
                        RadiusHistory(scale=SENS_R_SCALE_DOWN), base4_h, driver='R_scale_down'))

    by = {r['driver']: r for r in rows}
    span = {'R(t)': max(abs(by['R_scale_up']['rel_pct']), abs(by['R_scale_down']['rel_pct'])),
            'h_m': max(abs(by['hm_up']['rel_pct']), abs(by['hm_down']['rel_pct'])),
            'T_air': max(abs(by['T_air_up']['rel_pct']), abs(by['T_air_down']['rel_pct'])),
            'W_air': max(abs(by['W_air_up']['rel_pct']), abs(by['W_air_down']['rel_pct'])),
            'h': max(abs(by['h_up']['rel_pct']), abs(by['h_down']['rel_pct']))}
    ranking = sorted(span, key=lambda kk: -span[kk])
    elast_R = ((by['R_scale_up']['t_dry_h'] - by['R_scale_down']['t_dry_h']) / base4_h) / \
              (SENS_R_SCALE_UP - SENS_R_SCALE_DOWN)
    echo('[Q3-sens] 排序 ' + ' > '.join(f'{kk}({span[kk]:.3f}%)' for kk in ranking))
    echo(f'[Q3-sens] 半径弹性 dlnt/dlnR={elast_R:.3f}'
         f'（扩散标度理论值 {DIFFUSION_SCALING_EXPONENT:g}）；'
         f'h_m 双向不对称 {by["hm_up"]["rel_pct"]:+.3f}% vs {by["hm_down"]["rel_pct"]:+.3f}%')

    air_hold = AirState(hold_last_value=True)
    echo(f'[Q3-sens] 鲁棒性对照：空气延拓 {AIR_EXTRAPOLATION}={air0.T_tail_mean:.4f} °C vs '
         f'last_value={air_hold.T_hold:.4f} °C，差 {air_hold.T_hold - air0.T_tail_mean:.4f} °C '
         f'（已被 T_air±{SENS_DT_AIR} °C 行覆盖）')

    payload = {
        'meta': {'problem': 3, 'section': 'sensitivity', 'n_cells': N_ABL, 'dt_s': DT_Q3,
                 'recompute': 'full_simulation_per_row', 'n_rows': len(rows),
                 'runtime_s': round(time.time() - t_start, 2)},
        'baseline_q3_h': base_h, 'baseline_q4_h': base4_h,
        'rows': rows,
        'span_pct': span, 'ranking': ranking, 'elasticity_R': elast_R,
        'air_extrapolation_contrast': {'last24_mean_c': air0.T_tail_mean,
                                       'last_value_c': air_hold.T_hold,
                                       'delta_c': air_hold.T_hold - air0.T_tail_mean},
        'perturbations': {'dT_air_c': SENS_DT_AIR, 'dW_air_kgkg': SENS_DW_AIR,
                          'coef_scale': [SENS_SCALE_UP, SENS_SCALE_DOWN],
                          'R_scale': [SENS_R_SCALE_UP, SENS_R_SCALE_DOWN]},
    }
    validate_capability(payload)
    dump_json(os.path.join(FIG_DIR, 'sensitivity_results.json'), payload)
    echo(f'[Q3-sens] DONE {len(rows)} 行真实重算，总耗时 {time.time() - t_start:.1f}s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

"""结果约束审计：从已落盘的 figures/*_results.json 重新计算 R1-R24。

设计原则：**不信任求解器自报的任何布尔标志**。本脚本只读 JSON 里的数值，
按建模报告的约束清单以 lambda 逐条重算判据，再独立复核可疑数字（报告值与全精度值一致性、
四位小数落盘、跨问题数值不重复）。另附 MC1-MC14 方法签名静态检查——
因为 MODELING_REPORT.md 的 METHOD_CLAIMS_MACHINE 是 JSON 块，
_utils/claim_code_check.py 的 _parse_contract() 只识别 must:/forbid: 管道行，无法解析它。

输出只有结论：每条 PASS/FAIL、计数、最大误差、最多 5 条违规定位。不打印数据体。
"""
from __future__ import annotations

import json
import os
import re
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import numpy as np

from params import (ABL_PROPS_ONLY_MIN_PCT, ABL_SHRINK_ONLY_MAX_PCT, AIR_PROBE_LATE_S,
                    AIR_T_PLATEAU_C, ANALYTIC_DEV_TOL, C0_KGKG, C_CENTER_TOL, C_TH,
                    CM_PER_M, DECIMALS, DESC_WIDTH, DUP_DECIMALS,
                    DURATION_PRIOR_H, EXACT_TOL, FIG_DIR, GRAD_POS_TOL, GRID_REL_TOL,
                    K_A3_AT_C0_RANGE, N_A1_POINTS, N_A2_POINTS, N_FILE_RADII, PCT,
                    PLATEAU_REL_TOL, Q2_HORIZON_S, R_AT_6H_MAX_M, ROOT, SEC_PER_HOUR,
                    SERIES_DECIMALS, VOL_INV_MIN_DEV, VOL_INV_PROBE_H, ZETA_DEV_TOL)

MC_SPEC = [
    ('MC1', ['solve_banded|linalg\\.solve|tridiag', '(?i)flux|通量|a\\s*=\\s*coef'],
     ['exp\\(-\\s*k\\s*\\*\\s*t\\s*\\*\\*\\s*n', 'MR\\s*=\\s*np\\.exp\\(-']),
    ('MC2', ['d\\[-1\\]\\s*\\+=\\s*2\\.?0?\\s*\\*\\s*hm\\s*/\\s*R',
             'd\\[-1\\]\\s*\\+=\\s*2\\.?0?\\s*\\*\\s*h\\s*/\\s*R'],
     ['0\\.5\\s*\\*\\s*ds\\s*/\\s*D.*surface_extrap']),
    ('MC3', ['solve_banded', 'for\\s+_\\s+in\\s+range\\(\\s*nsweep|range\\(\\s*3\\s*\\)'],
     ['explicit_euler|forward_euler']),
    ('MC4', ['7e-9|7\\.0e-9|7e-09', 'exp\\(\\s*-\\s*0\\.89'], []),
    ('MC5', ['650', '128', '2736', '0\\.38', '2\\.4e-3', 'exp\\(\\s*-\\s*0\\.45'], []),
    ('MC6', ['760', '90', '2150', '0\\.20|0\\.2\\b', '4\\.2e-4',
             'exp\\(\\s*-\\s*0\\.30|exp\\(\\s*-\\s*0\\.3\\b'], []),
    ('MC7', ['273\\.15'], ['exp\\(\\s*-\\s*3850\\.?0?\\s*/\\s*T\\s*\\)']),
    ('MC8', ['C\\s*/\\s*\\(\\s*C\\s*\\+\\s*1'], []),
    ('MC9', ['4\\.?0?\\s*/\\s*R\\s*\\*\\*\\s*2', '/\\s*R\\s*\\)\\s*\\*\\*\\s*2'],
     ['front_tracking|moving_mesh_remap']),
    ('MC10', ['附件2|attachment2|a2', 'np\\.interp'], ['R_from_C|radius_from_moisture']),
    ('MC11', ['max\\(.*C|C.*\\.max\\('], ['mean\\(\\s*C\\s*\\)\\s*<\\s*0\\.15']),
    ('MC12', ['j0\\(', 'j1\\(', 'T_?air.*\\+.*b\\s*\\*|Tg\\[j\\+1\\]\\s*\\+'], []),
    ('MC13', ['\\[-24:\\]\\.mean\\(\\)', 'right\\s*='], ['fill_value\\s*=\\s*[\'"]extrapolate']),
    ('MC14', ['ablation|abl_'], []),
]


def load(name):
    with open(os.path.join(FIG_DIR, name), encoding='utf-8') as fh:
        return json.load(fh)


SELF = 'constraint_audit.py'


def load_code():
    """拼接 code/*.py 的可执行行（去掉整行注释，与 claim_code_check 的口径一致）。

    排除本文件：MC_SPEC 里的 forbid 正则以字符串字面量形式存在，
    若把审计器自身也当被审代码，会把"禁止模式的正则文本"误判成"命中禁止模式"。
    """
    out = []
    for fn in sorted(os.listdir(os.path.join(ROOT, 'code'))):
        if not fn.endswith('.py') or fn == SELF:
            continue
        with open(os.path.join(ROOT, 'code', fn), encoding='utf-8') as fh:
            for line in fh:
                if not line.strip().startswith('#'):
                    out.append(line)
    return ''.join(out)


def check(rid, desc, ok, detail):
    return {'id': rid, 'desc': desc, 'pass': bool(ok), 'detail': detail}


def audit_constraints(q1, q2, q3, q4, sens):
    """按约束清单逐条重算。所有取值来自 JSON 数值，不读任何 'pass' 字段。"""
    out = []
    T1 = np.asarray(q1['field_T_full'], dtype=float)
    C1 = np.asarray(q1['field_C_full'], dtype=float)
    lo_h, hi_h = DURATION_PRIOR_H
    k_lo, k_hi = K_A3_AT_C0_RANGE
    out.append(check('R1', f'q1_T_shape: T 列数>={N_FILE_RADII}',
                     T1.shape[1] >= N_FILE_RADII, f'shape={list(T1.shape)}'))
    out.append(check('R2', f'q1_C_center: |C(end,0)-{C0_KGKG}|<{C_CENTER_TOL:g}',
                     abs(C1[-1][0] - C0_KGKG) < C_CENTER_TOL, f'C={C1[-1][0]:.6f}'))
    dT = np.diff(T1[-1])
    out.append(check('R3', 'q1_T_monotone_r: 末行温度径向严格递增',
                     bool(np.all(dT > 0)), f'min_diff={dT.min():.6f}, n_violate={int((dT <= 0).sum())}'))
    out.append(check('R4', f'q1_analytic: 解析偏差<{ANALYTIC_DEV_TOL:g}',
                     q1['analytic']['analytic_max_dev'] < ANALYTIC_DEV_TOL,
                     f'dev={q1["analytic"]["analytic_max_dev"]:.3e}'))
    out.append(check('R5', f'q1_bessel: 特征根偏差<{ZETA_DEV_TOL:g}',
                     q1['analytic']['zeta_dev'] < ZETA_DEV_TOL, f'dev={q1["analytic"]["zeta_dev"]:.3e}'))
    k0 = q2['meta'][f'k_at_{C0_KGKG}']
    out.append(check('R6', f'q2_k_at_C0: {k_lo}<=k(C0)<={k_hi}', k_lo <= k0 <= k_hi, f'k={k0:.5f}'))
    out.append(check('R7', f'q2_rows: result2 数据行=={Q2_HORIZON_S}',
                     q2['xlsx_rows'] == Q2_HORIZON_S, f'rows={q2["xlsx_rows"]}'))
    out.append(check('R8', f'q2_air_points: 附件1=={N_A1_POINTS} 点',
                     q2['air']['n_points'] == N_A1_POINTS, f'n={q2["air"]["n_points"]}'))
    out.append(check('R9', f'q2_air_extrap: |T_air({AIR_PROBE_LATE_S:.0f})-{AIR_T_PLATEAU_C}|'
                           f'<{ANALYTIC_DEV_TOL:g}',
                     abs(q2['air']['T_at_20000'] - AIR_T_PLATEAU_C) < ANALYTIC_DEV_TOL,
                     f'T={q2["air"]["T_at_20000"]:.6f}'))
    out.append(check('R10', f'q3_tdry_range: {lo_h:.0f}<=t_dry<={hi_h:.0f} h',
                     lo_h <= q3['t_dry_h'] <= hi_h, f't_dry={q3["t_dry_h"]:.6f} h'))
    out.append(check('R11', f'q3_criterion: 末步 max_r C<{C_TH}',
                     q3['C_max_at_end'] < C_TH,
                     f'max={q3["C_max_at_end"]:.6f}, prev={q3["C_max_prev_step"]:.6f}'))
    out.append(check('R12', f'q3_monotone: 最大正径向梯度<{GRAD_POS_TOL:g}',
                     q3['max_pos_grad'] < GRAD_POS_TOL, f'max_pos_grad={q3["max_pos_grad"]:.3e}'))
    out.append(check('R13', f'q3_grid_conv: 生产网格相邻相对差<{GRID_REL_TOL:g}',
                     q3['grid_rel_diff'] < GRID_REL_TOL, f'rel={q3["grid_rel_diff"]:.3e}'))
    by = {r['driver']: r['rel_pct'] for r in sens['rows']}
    out.append(check('R14', 'q3_hm_signs: h_m 双向扰动反号',
                     by['hm_up'] * by['hm_down'] < 0.0,
                     f'up={by["hm_up"]:+.3f}%, down={by["hm_down"]:+.3f}%'))
    out.append(check('R15', f'q4_tdry_range: {lo_h:.0f}<=t_dry<={hi_h:.0f} h',
                     lo_h <= q4['t_dry_h'] <= hi_h, f't_dry={q4["t_dry_h"]:.6f} h'))
    out.append(check('R16', 'q4_tdry_lt_q3: t_dry(Q4)<t_dry(Q3)',
                     q4['t_dry_h'] < q3['t_dry_h'],
                     f'{q4["t_dry_h"]:.4f} < {q3["t_dry_h"]:.4f} h'))
    rh = q4['radius_history']
    out.append(check('R17', f'q4_R_points: 附件2=={N_A2_POINTS} 点',
                     rh['n_points'] == N_A2_POINTS, f'n={rh["n_points"]}'))
    dR = np.diff(np.asarray(rh['R_cm'], dtype=float))
    out.append(check('R18', 'q4_R_monotone: R(t) 单调不增',
                     bool(np.all(dR <= EXACT_TOL)), f'max_diff={dR.max():.3e} cm, '
                     f'n_violate={int((dR > EXACT_TOL).sum())}'))
    out.append(check('R19', f'q4_R_6h: R(6h)<{R_AT_6H_MAX_M * CM_PER_M:.2f} cm',
                     rh['R_6h_cm'] / CM_PER_M < R_AT_6H_MAX_M, f'R={rh["R_6h_cm"]:.4f} cm'))
    abl = q4['ablation']
    out.append(check('R20', f'q4_abl_shrink: 仅加收缩<{ABL_SHRINK_ONLY_MAX_PCT:.0f}%',
                     abl['shrink_only_pct'] < ABL_SHRINK_ONLY_MAX_PCT,
                     f'{abl["shrink_only_pct"]:+.2f}%'))
    out.append(check('R21', f'q4_abl_props: 仅换物性>+{ABL_PROPS_ONLY_MIN_PCT:.0f}%',
                     abl['props_only_pct'] > ABL_PROPS_ONLY_MIN_PCT,
                     f'{abl["props_only_pct"]:+.2f}%'))
    vi = q4['volume_inversion']
    out.append(check('R22', 'q4_inv_negative: 体积反演出现 C<0',
                     min(vi['C_vol']) < 0.0,
                     f'min={min(vi["C_vol"]):+.4f} @ {vi["first_negative_h"]:.0f} h'))
    out.append(check('R23', f'q4_inv_diverge: {VOL_INV_PROBE_H:.0f}h 偏离'
                            f'>{VOL_INV_MIN_DEV * PCT:.0f}%',
                     vi['dev_at_18h'] > VOL_INV_MIN_DEV,
                     f'{vi["dev_at_18h"] * PCT:.2f}%'))
    dec = [q1['meta']['decimals'], q2['meta']['decimals'],
           q3['meta']['decimals'], q4['meta']['decimals']]
    out.append(check('R24', f'decimals: 四问统一 {DECIMALS} 位小数',
                     all(d == DECIMALS for d in dec), f'decimals={dec}'))
    return out


def audit_suspicious(q1, q2, q3, q4, sens):
    """可疑数字复核：报告值与全精度一致、四位小数落盘、跨问题数值不重复。"""
    issues = []
    for tag, payload in (('q3', q3), ('q4', q4)):
        if abs(payload['t_dry_reported'] - round(payload['t_dry_h'], DECIMALS)) > EXACT_TOL:
            issues.append(f'{tag}.t_dry_reported 与全精度四舍五入不一致')
    for tag, payload, key in (('q1', q1, 'table1_T'), ('q1', q1, 'table2_C'),
                              ('q2', q2, 'table3_T'), ('q2', q2, 'table4_C'),
                              ('q3', q3, 'table5_C'), ('q3', q3, 'table5_T')):
        arr = np.asarray(payload[key], dtype=float)
        bad = int(np.count_nonzero(np.round(arr, DECIMALS) != arr))
        if bad:
            issues.append(f'{tag}.{key} 有 {bad} 个非四位小数值')
    if abs(q3['t_dry_h'] - q4['t_dry_h']) < GRAD_POS_TOL:
        issues.append('q3/q4 的 t_dry 完全相同，疑似复制')
    grid3 = q3['grid_convergence']['t_dry_h']
    if len(set(round(v, DUP_DECIMALS) for v in grid3)) != len(grid3):
        issues.append('q3 网格表出现重复 t_dry，疑非真实重算')
    vals = [r['t_dry_h'] for r in sens['rows']]
    if len(set(round(v, SERIES_DECIMALS) for v in vals)) != len(vals):
        issues.append('灵敏度表出现重复 t_dry，疑非真实重算')
    abl = q4['ablation']
    if abl['props_only_pct'] * abl['shrink_only_pct'] >= 0.0:
        issues.append('消融两主效应未反号')
    xl3 = q3['xlsx_rows'] * q34_step(q3)
    if xl3 < q3['t_dry_s']:
        issues.append('result3.xlsx 覆盖时长短于 t_dry')
    if q4['xlsx_rows'] * q34_step(q4) < q4['t_dry_s']:
        issues.append('result4.xlsx 覆盖时长短于 t_dry')
    return issues


def q34_step(payload):
    t = payload['series_t_h']
    return (t[1] - t[0]) * SEC_PER_HOUR


def audit_methods(code):
    """MC1-MC14：must 全部命中、forbid 全部不命中。"""
    rows = []
    for mid, must, forbid in MC_SPEC:
        miss = [p for p in must if not re.search(p, code)]
        hit = [p for p in forbid if re.search(p, code)]
        rows.append({'id': mid, 'pass': not miss and not hit,
                     'missing_must': miss, 'hit_forbid': hit})
    return rows


def main() -> int:
    q1, q2 = load('problem_1_results.json'), load('problem_2_results.json')
    q3, q4 = load('problem_3_results.json'), load('problem_4_results.json')
    sens, sanity = load('sensitivity_results.json'), load('sanity_check_results.json')

    cons = audit_constraints(q1, q2, q3, q4, sens)
    n_fail = sum(1 for c in cons if not c['pass'])
    for c in cons:
        print(f'  {c["id"]:<4s} {"PASS" if c["pass"] else "FAIL"}  '
              f'{c["desc"]:<{DESC_WIDTH}s} {c["detail"]}')
    print(f'CONSTRAINTS n={len(cons)} pass={len(cons) - n_fail} fail={n_fail}')

    issues = audit_suspicious(q1, q2, q3, q4, sens)
    for msg in issues[:5]:
        print(f'  SUSPICIOUS {msg}')
    print(f'SUSPICIOUS_NUMBERS n={len(issues)}')

    mc = audit_methods(load_code())
    mc_fail = [m for m in mc if not m['pass']]
    for m in mc_fail[:5]:
        print(f'  {m["id"]} FAIL missing={m["missing_must"]} forbid_hit={m["hit_forbid"]}')
    print(f'METHOD_CLAIMS n={len(mc)} implemented={len(mc) - len(mc_fail)}')

    struct_ok = sanity['all_pass']
    print(f'STRUCTURAL_CHECKS n={sanity["meta"]["n_checks"]} pass={sanity["meta"]["n_pass"]}')

    payload = {
        'constraints': cons, 'n_constraints': len(cons), 'n_fail': n_fail,
        'suspicious_numbers': issues, 'n_suspicious': len(issues),
        'method_claims': mc, 'n_method_claims': len(mc),
        'n_method_implemented': len(mc) - len(mc_fail),
        'structural_checks_pass': struct_ok,
        'verdict': 'PASS' if (n_fail == 0 and not issues and not mc_fail and struct_ok) else 'FAIL',
    }
    with open(os.path.join(FIG_DIR, 'constraint_audit.json'), 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1)
    print(f'CONSTRAINT_AUDIT {payload["verdict"]}')
    return 0 if payload['verdict'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())

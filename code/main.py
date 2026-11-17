"""总编排：按依赖顺序跑完四问 + 灵敏度 + 结构性验证 + 约束审计，并汇总 all_results.json。

依赖关系（决定顺序，不可颠倒）：
  data_check → problem1/2/3 → problem3_sens → problem4 → sanity_check（读 Q3/Q4 的 JSON）
  → constraint_audit（读全部 JSON，重算 R1-R24）
每步在子进程中跑，任一步非零退出即整体失败——不吞异常、不跳过。
用法：python code/main.py [--skip-solved]，后者只在对应 JSON 已存在时跳过重算。
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from params import BYTES_PER_KB, FIG_DIR, LOG_NAME_WIDTH, PCT, ROOT

STEPS = [
    ('data_check.py', [], '数据摄入核对'),
    ('problem1.py', ['problem_1_results.json'], '问题一：预热阶段常物性场 + 解析校验'),
    ('problem2.py', ['problem_2_results.json'], '问题二：全程变物性耦合场'),
    ('problem3.py', ['problem_3_results.json'], '问题三：烘干时长逐点阈值反问题'),
    ('problem3_sens.py', ['sensitivity_results.json'], '灵敏度：逐行真实重算'),
    ('problem4.py', ['problem_4_results.json'], '问题四：附录 4 + 附件 2 收缩'),
    ('sanity_check.py', ['sanity_check_results.json'], '结构性验证四项'),
    ('constraint_audit.py', ['constraint_audit.json'], '约束审计 R1-R24 + MC1-MC14'),
]

AGGREGATE = [('problem1', 'problem_1_results.json'), ('problem2', 'problem_2_results.json'),
             ('problem3', 'problem_3_results.json'), ('problem4', 'problem_4_results.json'),
             ('sensitivity', 'sensitivity_results.json'),
             ('sanity_check', 'sanity_check_results.json'),
             ('constraint_audit', 'constraint_audit.json')]


def run_step(script, outputs, desc, skip_solved):
    path = os.path.join(ROOT, 'code', script)
    done = outputs and all(os.path.isfile(os.path.join(FIG_DIR, o)) for o in outputs)
    if skip_solved and done:
        print(f'[main] SKIP  {script:<{LOG_NAME_WIDTH}s} {desc}（输出已存在）', flush=True)
        return 0.0
    print(f'[main] RUN   {script:<{LOG_NAME_WIDTH}s} {desc}', flush=True)
    t0 = time.time()
    proc = subprocess.run([sys.executable, path], cwd=ROOT)
    dt = time.time() - t0
    if proc.returncode != 0:
        raise SystemExit(f'[main] FAIL  {script} 退出码 {proc.returncode}')
    print(f'[main] OK    {script:<{LOG_NAME_WIDTH}s} {dt:.1f}s', flush=True)
    return dt


def aggregate():
    """汇总各问 JSON 为 figures/all_results.json（保持全精度，不做二次舍入）。"""
    all_res = {}
    for key, fname in AGGREGATE:
        with open(os.path.join(FIG_DIR, fname), encoding='utf-8') as fh:
            all_res[key] = json.load(fh)
    q3, q4 = all_res['problem3'], all_res['problem4']
    all_res['summary'] = {
        'q1_T_center_1800s': all_res['problem1']['table1_T'][-1][0],
        'q1_T_surface_1800s': all_res['problem1']['table1_T'][-1][-1],
        'q1_analytic_max_dev_c': all_res['problem1']['analytic']['analytic_max_dev'],
        'q2_T_center_3h': all_res['problem2']['table3_T'][-1][0],
        'q2_C_center_3h': all_res['problem2']['table4_C'][-1][0],
        'q3_t_dry_h': q3['t_dry_h'], 'q3_t_dry_d': q3['t_dry_d'],
        'q4_t_dry_h': q4['t_dry_h'], 'q4_t_dry_d': q4['t_dry_d'],
        'q4_vs_q3_pct': (q4['t_dry_h'] - q3['t_dry_h']) / q3['t_dry_h'] * PCT,
        'ablation': {k: v for k, v in q4['ablation'].items() if k.endswith(('_h', '_pct'))},
        'sensitivity_ranking': all_res['sensitivity']['ranking'],
        'conservation_residual_max': max(all_res[f'problem{i}']['conservation_residual']
                                         for i in (1, 2, 3, 4)),
        'constraints_pass': all_res['constraint_audit']['n_constraints']
        - all_res['constraint_audit']['n_fail'],
        'constraints_total': all_res['constraint_audit']['n_constraints'],
        'method_claims_implemented': all_res['constraint_audit']['n_method_implemented'],
        'method_claims_total': all_res['constraint_audit']['n_method_claims'],
        'verdict': all_res['constraint_audit']['verdict'],
    }
    out = os.path.join(FIG_DIR, 'all_results.json')
    with open(out, 'w', encoding='utf-8') as fh:
        json.dump(all_res, fh, ensure_ascii=False, indent=1, allow_nan=False)
    size_kb = os.path.getsize(out) / BYTES_PER_KB
    print(f'[main] all_results.json 已写：{len(AGGREGATE)} 段，{size_kb:.0f} KB', flush=True)
    return all_res['summary']


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--skip-solved', action='store_true',
                    help='对应 JSON 已存在时跳过该步（用于只重跑汇总）')
    args = ap.parse_args()
    t_start = time.time()
    for script, outputs, desc in STEPS:
        run_step(script, outputs, desc, args.skip_solved)
    summary = aggregate()
    print(f'[main] 汇总：Q3 t_dry={summary["q3_t_dry_h"]:.4f} h，'
          f'Q4 t_dry={summary["q4_t_dry_h"]:.4f} h（{summary["q4_vs_q3_pct"]:+.1f}%），'
          f'约束 {summary["constraints_pass"]}/{summary["constraints_total"]}，'
          f'方法签名 {summary["method_claims_implemented"]}'
          f'/{summary["method_claims_total"]}，'
          f'守恒残差上界 {summary["conservation_residual_max"]:.2e}')
    print(f'[main] ALL DONE 总耗时 {time.time() - t_start:.1f}s，判定 {summary["verdict"]}')
    return 0 if summary['verdict'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())

"""数据读取验证：附件1（烘房温湿度）、附件2（药材半径）、result*.xlsx 模板。

先跑本脚本确认读取口径，再写求解代码。只打印结论与小规模摘要。
"""
from __future__ import annotations

import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import numpy as np
import pandas as pd

from params import *  # noqa: F401,F403  ROOT / PROFILE / 口径常数

RULE_BAR = '=' * RULE_WIDTH


def read_all_sheets(rel_path: str) -> dict:
    """Excel 一律读全部 sheet（不写 sheet_name 时 pandas 只读首表且不告警）。"""
    path = os.path.join(ROOT, rel_path)
    return pd.read_excel(path, sheet_name=None)


def check_profile(fname: str, total_rows: int) -> None:
    prof = PROFILE.get(fname, {})
    assert 'total_rows' in prof, (
        f"[未建档] {fname} 不在 DATA_PROFILE.json 或建档失败（{prof.get('error', '缺记录')}）")
    assert total_rows == prof['total_rows'], (
        f"[摄入不全] {fname} 实读 {total_rows} 行 != 建档 {prof['total_rows']} 行"
        f"（该文件共 {prof['n_sheets']} 张 sheet）")
    print(f"[data_ingest] {fname} 全量 {total_rows} 行 / {prof['n_sheets']} sheet 已核对通过")


def main() -> int:
    print(RULE_BAR)
    print("附件1：烘房温度与水分浓度")
    sheets1 = read_all_sheets('user_data/附件1.xlsx')
    print(f"  sheet 数 {len(sheets1)}: " + ", ".join(f"{n}({len(d)}行)" for n, d in sheets1.items()))
    a1 = sheets1['Sheet1']
    print(f"  列名 {list(a1.columns)}, dtypes {[str(t) for t in a1.dtypes]}")
    print(f"  缺失值 {int(a1.isnull().sum().sum())}")
    t1 = a1['时间'].to_numpy(float)
    Tair = a1['温度'].to_numpy(float)
    Wair = a1['水分浓度'].to_numpy(float)
    check_profile('附件1.xlsx', len(a1))
    dt1 = np.diff(t1)
    print(f"  t 范围 [{t1[0]:.0f}, {t1[-1]:.0f}] s, 间隔唯一值 {np.unique(dt1)}")
    print(f"  T 范围 [{Tair.min():.4f}, {Tair.max():.4f}] C, 首值 {Tair[0]:.6f}, 末值 {Tair[-1]:.6f}")
    print(f"  W 范围 [{Wair.min():.5f}, {Wair.max():.5f}] kg/kg, 首值 {Wair[0]:.6f}, 末值 {Wair[-1]:.6f}")
    tail = slice(-AIR_TAIL_POINTS, None)
    print(f"  末 {AIR_TAIL_POINTS} 点均值: T={Tair[tail].mean():.6f} C, "
          f"W={Wair[tail].mean():.6f} kg/kg")
    print(f"  末 {AIR_TAIL_POINTS} 点 T 波动 "
          f"+-{HALF * (Tair[tail].max() - Tair[tail].min()):.4f} C")

    print(RULE_BAR)
    print("附件2：药材半径")
    sheets2 = read_all_sheets('user_data/附件2.xlsx')
    print(f"  sheet 数 {len(sheets2)}: " + ", ".join(f"{n}({len(d)}行)" for n, d in sheets2.items()))
    a2 = sheets2['Sheet1']
    print(f"  列名 {list(a2.columns)}, dtypes {[str(t) for t in a2.dtypes]}")
    print(f"  缺失值 {int(a2.isnull().sum().sum())}")
    t2 = a2['时间'].to_numpy(float)
    Rcm = a2['半径'].to_numpy(float)
    check_profile('附件2.xlsx', len(a2))
    print(f"  t 范围 [{t2[0]:.0f}, {t2[-1]:.0f}] s ({t2[-1] / SEC_PER_HOUR:.1f} h), "
          f"间隔唯一值 {np.unique(np.diff(t2))}")
    print(f"  R 范围 [{Rcm.min():.4f}, {Rcm.max():.4f}] cm, R(0)={Rcm[0]:.6f}, R(end)={Rcm[-1]:.6f}")
    dR = np.diff(Rcm)
    print(f"  单调非增: {bool(np.all(dR <= EXACT_TOL))}, 最大正增量 {dR.max():.3e} cm")
    dt_probe = Q34_TABLE_DT_H * SEC_PER_HOUR
    idx6 = int(np.searchsorted(t2, dt_probe))
    print(f"  R({Q34_TABLE_DT_H:.0f}h)={Rcm[idx6]:.6f} cm (t={t2[idx6]:.0f} s), "
          f"末端面积比 (R/R0)^2={(Rcm[-1] / Rcm[0]) ** 2:.6f}")
    per6 = [Rcm[int(np.searchsorted(t2, k * dt_probe))] for k in range(N_VI_INTERVALS + 1)]
    print(f"  每 {Q34_TABLE_DT_H:.0f} h 的 R(cm): " + ", ".join(f"{v:.3f}" for v in per6))

    print(RULE_BAR)
    print("result*.xlsx 模板表头")
    for f in ['result1.xlsx', 'result2.xlsx', 'result3.xlsx', 'result4.xlsx']:
        sheets = read_all_sheets('user_data/' + f)
        check_profile(f, sum(len(d) for d in sheets.values()))
        for name, df in sheets.items():
            print(f"  {f}::{name} 表头 {list(df.columns)}")

    print(RULE_BAR)
    print(f"DATA_CHECK 结论: 附件1 {len(a1)} 点 / 附件2 {len(a2)} 点 / 模板表头已核对，读取口径无误")
    assert len(a1) == N_A1_POINTS and len(a2) == N_A2_POINTS, '实读点数与行数账本不符'
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

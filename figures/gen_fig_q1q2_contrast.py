"""问题一与问题二在同一时刻的口径对照（两面板分组柱：温度 / 含水率）。

本图讲什么：t=1800 s 同一时刻、同一 5 个输出半径上，附录 2 常物性闭合（问题一）
与附录 3 变物性闭合（问题二）的读数并排对照——两题结论不同源于物性闭合不同，
不是相互矛盾。
(a) 温度对照；(b) 含水率对照。柱顶数值为该位置读数，虚线为初始含水率参考。
数据来源：附录 2 侧取 figures/problem_1_results.json → table1_T / table2_C 的
末行（对应 table_times_s 的 1800 s）；附录 3 侧取 problem_2_results.json →
field_T_full / field_C_full 在 profile_t_s 中 t=1800 s 精确命中行，按 5 个输出
半径抽取（径向网格 problem_1_results.json → file_radii_cm）。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS, _lighten

setup_style()

d1 = fb.load('problem_1_results')
d2 = fb.load('problem_2_results')

t_tab = fb.arr(d1['table_times_s'])
i1 = int(np.argmin(np.abs(t_tab - 1800.0)))
T1 = fb.mat(d1['table1_T'])[i1]
C1 = fb.mat(d1['table2_C'])[i1]

r_cm = fb.arr(d1['file_radii_cm'])
r_out = np.linspace(0.0, float(r_cm[-1]), T1.size)
cols = [int(np.argmin(np.abs(r_cm - rr))) for rr in r_out]
t2 = fb.arr(d2['profile_t_s'])
i2 = int(np.argmin(np.abs(t2 - 1800.0)))
assert abs(float(t2[i2]) - 1800.0) < 1e-6, '问题二时间轴未精确命中 1800 s'
T2 = fb.mat(d2['field_T_full'])[i2][cols]
C2 = fb.mat(d2['field_C_full'])[i2][cols]

labels = [f'{rr:.1f}' for rr in r_out]
x = np.arange(len(labels))
w = 0.36
CA, CB = PALETTE[2], PALETTE[0]

fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.1), layout='constrained')
fb.set_paper_placement(fig)

for ax, va, vb, ylab, tag, fmt in (
        (axes[0], T1, T2, '温度 $T$ (°C)', '(a) 温度对照', '%.1f'),
        (axes[1], C1, C2, '干基含水率 $C$ (kg/kg)', '(b) 含水率对照', '%.2f')):
    b1 = ax.bar(x - w / 2, va, w, color=_lighten(CA, 0.45), edgecolor=CA,
                linewidth=1.0, label='附录 2 常物性')
    b2 = ax.bar(x + w / 2, vb, w, color=_lighten(CB, 0.45), edgecolor=CB,
                linewidth=1.0, label='附录 3 变物性')
    ax.bar_label(b1, fmt=fmt, padding=1.5, fontsize=8.0, color=CA, rotation=90)
    ax.bar_label(b2, fmt=fmt, padding=1.5, fontsize=8.0, color=CB, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel('径向坐标 $r$ (cm)')
    ax.set_ylabel(ylab)
    top = max(float(np.nanmax(va)), float(np.nanmax(vb)))
    ax.set_ylim(0.0, top * 1.30)
    fb.declutter_axes(ax, grid=False)
    fb.panel(ax, tag)

fb.consolidate_shared_legends(fig, list(axes), where='top', ncol=2)

save_fig(fig, fb.out('fig_q1q2_contrast'))

"""附件 2 药材半径收缩历程（单面板双轴）。

本图讲什么：药材外半径 R(t) 随干燥时间的单调收缩过程，以及按 6 h 分箱统计
的收缩增量分布——用于说明收缩集中在前 6 h、后期几何冻结。
左轴：外半径 R (cm) 折线；右轴：每 6 h 区间的半径减少量 (cm) 柱。
数据来源：figures/problem_4_results.json → radius_history.t_s / radius_history.R_cm
（n=145 全量，按时间升序排序后使用；R0 与末值取自同一数组）。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS, _lighten

setup_style()

rh = fb.load('problem_4_results')['radius_history']
t_raw = fb.arr(rh['t_s']) / 3600.0
R_raw = fb.arr(rh['R_cm'])
order = np.argsort(t_raw)
t_h, R = t_raw[order], R_raw[order]

BIN_H = 6.0
edges = np.arange(0.0, float(t_h[-1]) + BIN_H, BIN_H)
centers, drops = [], []
for lo, hi in zip(edges[:-1], edges[1:]):
    R_lo = float(np.interp(lo, t_h, R))
    R_hi = float(np.interp(min(hi, t_h[-1]), t_h, R))
    centers.append(0.5 * (lo + min(hi, t_h[-1])))
    drops.append(R_lo - R_hi)
centers = np.array(centers)
drops = np.array(drops)

C_R, C_D = PALETTE[2], PALETTE[0]
total_drop = float(R[0] - R[-1])
first_share = drops[0] / total_drop * 100.0 if total_drop > 0 else np.nan

fig, ax = plt.subplots(figsize=(6.4, 3.9))
fb.set_paper_placement(fig)
ax2 = ax.twinx()

ax2.bar(centers, drops, width=BIN_H * 0.72, color=_lighten(C_D, 0.45),
        edgecolor=C_D, linewidth=0.9, zorder=1, label='每 6 h 收缩量')
ax.plot(t_h, R, color=C_R, linewidth=1.6, zorder=3, label='外半径 $R(t)$')
ax.scatter([t_h[0], t_h[-1]], [R[0], R[-1]], s=22, color=C_R, zorder=4,
           edgecolors='white', linewidths=0.8)

ax.annotate(f'{R[0]:.3f} cm', xy=(t_h[0], R[0]), xytext=(5, -2),
            textcoords='offset points', ha='left', va='top',
            fontsize=8.5, color=C_R)
ax.annotate(f'{R[-1]:.3f} cm', xy=(t_h[-1], R[-1]), xytext=(-3, 6),
            textcoords='offset points', ha='right', va='bottom',
            fontsize=8.5, color=C_R)
ax2.annotate(f'首段占 {first_share:.0f}%', xy=(centers[0], drops[0]),
             xytext=(6, 2), textcoords='offset points', ha='left', va='bottom',
             fontsize=8.5, color=C_D)

ax.set_xlabel('干燥时间 $t$ (h)')
ax.set_ylabel('外半径 $R$ (cm)', color=C_R)
ax2.set_ylabel('每 6 h 半径减少量 (cm)', color=C_D)
ax.tick_params(axis='y', colors=C_R)
ax2.tick_params(axis='y', colors=C_D)
ax.set_xlim(-1.5, float(t_h[-1]) + 1.5)
ax.set_ylim(float(R.min()) - 0.08, float(R.max()) + 0.06)
ax2.set_ylim(0.0, float(drops.max()) * 1.45)
ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.spines['left'].set_color(C_R)
ax2.spines['right'].set_color(C_D)
ax.set_zorder(ax2.get_zorder() + 1)
ax.patch.set_visible(False)

handles = ax.get_lines()[:1] + list(ax2.containers)
fb.auto_legend(ax, handles=handles,
               labels=[h.get_label() for h in handles], loc='upper right')

save_fig(fig, fb.out('fig_shrinkage_curve'))

"""收缩与物性的加性分解瀑布图（单面板）。

本图讲什么：从问题三基准工况到问题四完整工况，干燥时间 t_dry 的变化如何分解为
「仅开收缩」「仅换附录 4 物性」与两者的交互项三部分。色带层显示每一步的增量区间，
阶梯线与圆点显示累积值，圆点上方标出各步累积 t_dry。两个单因子效应符号相反且量级
相近，交互项吸收其非线性叠加残差，四项之和精确回到完整工况。
数据来源：figures/problem_4_results.json → ablation.baseline_h /
.shrink_only_h / .props_only_h / .full_h / .interaction_h（增量在本脚本内由
这些绝对值相减得到，未额外引入数字）。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS, _lighten

setup_style()

ab = fb.load('problem_4_results')['ablation']
base = float(ab['baseline_h'])
d_shrink = float(ab['shrink_only_h']) - base
d_props = float(ab['props_only_h']) - base
d_inter = float(ab['interaction_h'])
full = float(ab['full_h'])
assert abs(base + d_shrink + d_props + d_inter - full) < 1e-6, '加性分解不闭合'

labels = ['基准', '+收缩', '+变物性', '+交互项']
deltas = [d_shrink, d_props, d_inter]
cum = [base]
for dv in deltas:
    cum.append(cum[-1] + dv)
cum = np.array(cum)
assert abs(cum[-1] - full) < 1e-6
x = np.arange(len(labels))
layer_cols = [COLORS['down'] if dv < 0 else COLORS['up'] for dv in deltas]

fig, ax = plt.subplots(figsize=(6.4, 3.9))
fb.set_paper_placement(fig)

ax.fill_between([x[0] - 0.5, x[-1] + 0.5], 0.0, cum[0],
                color=_lighten(PALETTE[2], 0.88), zorder=0)
for i, col in enumerate(layer_cols, start=1):
    lo, hi = min(cum[i - 1], cum[i]), max(cum[i - 1], cum[i])
    ax.fill_between([x[i] - 0.5, x[-1] + 0.5], lo, hi,
                    color=_lighten(col, 0.72), zorder=i)
    ax.plot([x[i] - 0.5, x[-1] + 0.5], [cum[i], cum[i]], color=col,
            linewidth=0.8, linestyle='--', dashes=(4, 2), zorder=i)

ax.step(x, cum, where='mid', color=PALETTE[2], linewidth=2.2, zorder=10)
for i in range(len(labels)):
    col = PALETTE[2] if i == 0 else layer_cols[min(i - 1, len(layer_cols) - 1)]
    ax.scatter(x[i], cum[i], s=52, color=col, zorder=11, edgecolors='white',
               linewidths=1.6)
    ax.annotate(f'{cum[i]:.2f}', xy=(x[i], cum[i]), xytext=(0, 8),
                textcoords='offset points', ha='center', va='bottom',
                fontsize=8.5, color=col, zorder=12,
                bbox=dict(boxstyle='square,pad=0.14', facecolor='white',
                          edgecolor='none'))

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel('干燥完成时刻 $t_{\\mathrm{dry}}$ (h)')
ax.set_xlim(x[0] - 0.5, x[-1] + 0.5)
ax.set_ylim(0.0, float(cum.max()) * 1.16)
fb.declutter_axes(ax, grid=False)

save_fig(fig, fb.out('fig_q4_ablation'))

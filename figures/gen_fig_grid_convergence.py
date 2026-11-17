"""网格收敛性检验（单面板，问题三/问题四双曲线 + Richardson 外推）。

本图讲什么：干燥完成时刻 t_dry 随径向控制体数 N 的收敛行为。两条曲线分别为
问题三（固定半径）与问题四（移动边界）；水平虚线为按最后三档做 Richardson 外推
得到的网格无关极限；阴影为最细两档之间的相对差带（即残余离散误差的量级）。
右侧竖直短线标出各题最终采用的 N。
数据来源：figures/problem_3_results.json → grid_convergence.N / .t_dry_h /
.rel_diff 与 dt_convergence.n_cells；figures/problem_4_results.json →
grid_convergence.N / .t_dry_h。Richardson 外推在本脚本内按二阶假定现算。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS, _lighten

setup_style()

d3 = fb.load('problem_3_results')
d4 = fb.load('problem_4_results')
g3, g4 = d3['grid_convergence'], d4['grid_convergence']
N3, T3 = fb.arr(g3['N']), fb.arr(g3['t_dry_h'])
N4, T4 = fb.arr(g4['N']), fb.arr(g4['t_dry_h'])


def richardson(N, T):
    """相邻两档等比细化下的二阶 Richardson 外推极限。"""
    r = float(N[-1]) / float(N[-2])
    return float(T[-1] + (T[-1] - T[-2]) / (r ** 2 - 1.0))


R3, R4 = richardson(N3, T3), richardson(N4, T4)

fig, ax = plt.subplots(figsize=(6.4, 3.9))
fb.set_paper_placement(fig)

for N, T, R, col, name in ((N3, T3, R3, PALETTE[2], '问题三 固定半径'),
                           (N4, T4, R4, PALETTE[0], '问题四 移动边界')):
    band = abs(float(T[-1] - T[-2]))
    fb.uncertainty_band(ax, N, np.full_like(N, R - band),
                        np.full_like(N, R + band), color=col, alpha=0.13)
    ax.plot(N, T, color=col, linewidth=1.5, marker='o', markersize=4.2,
            markerfacecolor='white', markeredgewidth=1.2, label=name)
    ax.plot([float(N[0]) * 0.85, float(N[-1])], [R, R], color=col,
            linewidth=0.9, linestyle='--', dashes=(5, 2.5))
    ax.annotate(f'{R:.3f} h', xy=(float(N[-1]), R), xytext=(7, 0),
                textcoords='offset points', ha='left', va='center',
                fontsize=8.5, color=col,
                bbox=dict(boxstyle='square,pad=0.14', facecolor='white',
                          edgecolor='none'))

N_used = int(d3['dt_convergence']['n_cells'])
ax.axvline(N_used, color=COLORS['ref_line'], linewidth=0.9, linestyle=':',
           zorder=1)
ax.annotate(f'$N$={N_used}', xy=(N_used, float(np.nanmin(T4))),
            xytext=(4, 6), textcoords='offset points', ha='left', va='bottom',
            fontsize=8.5, color=COLORS['text'])

ax.set_xlabel('径向控制体数 $N$')
ax.set_ylabel('干燥完成时刻 $t_{\\mathrm{dry}}$ (h)')
ax.set_xlim(float(min(N3.min(), N4.min())) * 0.85,
            float(max(N3.max(), N4.max())) * 1.14)
fb.declutter_axes(ax, grid=False)
fb.auto_legend(ax, loc='center right')

save_fig(fig, fb.out('fig_grid_convergence'))

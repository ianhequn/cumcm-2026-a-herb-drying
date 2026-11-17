"""问题二温度场 T(r,t) 云图（单面板，附录 3 变物性，0-3 h）。

本图讲什么：物性随含水率与温度变化后的温度场分布，横轴时间（h）、纵轴半径，
填充色为温度并叠加等值线；与问题一同色标口径便于跨题对读。
数据来源：figures/problem_2_results.json → field_T_full (181×21 全量)、
profile_t_s（换算为小时）、problem_1_results.json → file_radii_cm（同一径向网格）。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS

setup_style()

d2 = fb.load('problem_2_results')
t_h = fb.arr(d2['profile_t_s']) / 3600.0
r_cm = fb.arr(fb.load('problem_1_results')['file_radii_cm'])
T = fb.mat(d2['field_T_full'])

fig, ax = plt.subplots(figsize=(6.4, 3.8), layout='constrained')
fb.set_paper_placement(fig)

levels = np.linspace(float(np.nanmin(T)), float(np.nanmax(T)), 15)
cf = ax.contourf(t_h, r_cm, T.T, levels=levels, cmap=fb.CMAP_T)
cs = ax.contour(t_h, r_cm, T.T, levels=levels[1::3], colors=[COLORS['dark']],
                linewidths=0.45, alpha=0.55)
fb.label_contours(ax, cs, '%.0f')

cb = fig.colorbar(cf, ax=ax, pad=0.02, aspect=26)
cb.set_label('温度 $T$ (°C)')
cb.ax.tick_params(length=2.0)

t_q1 = 1800.0 / 3600.0
ax.axvline(t_q1, color=COLORS['text'], linewidth=0.9, linestyle='--',
           dashes=(4, 2))
fb.field_note(ax, '$t$=1800 s', (t_q1, float(r_cm[-1])), xytext=(4, -6),
              ha='left')

ax.set_xlabel('干燥时间 $t$ (h)')
ax.set_ylabel('径向坐标 $r$ (cm)')
ax.set_xlim(float(t_h[0]), float(t_h[-1]))
ax.set_ylim(0.0, float(r_cm[-1]))

save_fig(fig, fb.out('fig_q2_field_T'))

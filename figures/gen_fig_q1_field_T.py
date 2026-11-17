"""问题一温度场 T(r,t) 云图（单面板，附录 2 常物性，0-1800 s）。

本图讲什么：短时程内圆柱药材内部温度的二维分布——横轴时间、纵轴半径，
填充色为温度，叠加等值线读数；虚线标出中心与表面两条边界位置。
数据来源：figures/problem_1_results.json → field_T_full (1800×21 全量)
配合 profile_t_s (时间轴) 与 file_radii_cm (半径轴)。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS

setup_style()

d = fb.load('problem_1_results')
t_s = fb.arr(d['profile_t_s'])
r_cm = fb.arr(d['file_radii_cm'])
T = fb.mat(d['field_T_full'])

fig, ax = plt.subplots(figsize=(6.4, 3.8), layout='constrained')
fb.set_paper_placement(fig)

levels = np.linspace(float(np.nanmin(T)), float(np.nanmax(T)), 15)
cf = ax.contourf(t_s, r_cm, T.T, levels=levels, cmap=fb.CMAP_T)
cs = ax.contour(t_s, r_cm, T.T, levels=levels[1::3], colors=[COLORS['dark']],
                linewidths=0.45, alpha=0.55)
fb.label_contours(ax, cs, '%.0f')

cb = fig.colorbar(cf, ax=ax, pad=0.02, aspect=26)
cb.set_label('温度 $T$ (°C)')
cb.ax.tick_params(length=2.0)

ax.axhline(float(r_cm[-1]), color=COLORS['text'], linewidth=0.8,
           linestyle='--', dashes=(4, 2))
fb.field_note(ax, f'表面 $r$={r_cm[-1]:.1f} cm',
              (float(t_s[-1]), float(r_cm[-1])), xytext=(-4, -6), ha='right')

ax.set_xlabel('干燥时间 $t$ (s)')
ax.set_ylabel('径向坐标 $r$ (cm)')
ax.set_xlim(float(t_s[0]), float(t_s[-1]))
ax.set_ylim(0.0, float(r_cm[-1]))

save_fig(fig, fb.out('fig_q1_field_T'))

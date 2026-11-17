"""问题二含水率场 C(r,t) 云图（单面板，附录 3 变物性，0-3 h）。

本图讲什么：变物性闭合下 3 h 内的含水率二维分布。横轴时间（h）、纵轴半径，
填充色为干基含水率并叠加等值线；虚线标出末时刻湿分影响前沿（由场数据反算：
偏离初值超过 1% 的最内侧半径）。
数据来源：figures/problem_2_results.json → field_C_full (181×21 全量)、
profile_t_s；径向网格取 problem_1_results.json → file_radii_cm。
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
C = fb.mat(d2['field_C_full'])

C0 = float(C[0, 0])
moved = np.where(np.abs(C[-1] - C0) > 0.01 * C0)[0]
r_front = float(r_cm[moved[0]]) if moved.size else float(r_cm[-1])
depth_mm = (float(r_cm[-1]) - r_front) * 10.0

fig, ax = plt.subplots(figsize=(6.4, 3.8), layout='constrained')
fb.set_paper_placement(fig)

levels = np.linspace(float(np.nanmin(C)), float(np.nanmax(C)), 15)
cf = ax.contourf(t_h, r_cm, C.T, levels=levels, cmap=fb.CMAP_C)
cs = ax.contour(t_h, r_cm, C.T, levels=levels[1::3], colors=[COLORS['dark']],
                linewidths=0.45, alpha=0.55)
fb.label_contours(ax, cs, '%.2f')

cb = fig.colorbar(cf, ax=ax, pad=0.02, aspect=26)
cb.set_label('干基含水率 $C$ (kg/kg)')
cb.ax.tick_params(length=2.0)

ax.axhline(r_front, color=COLORS['text'], linewidth=0.9, linestyle='--',
           dashes=(4, 2))
# 内区含水率仍为初值、无等值线，把读数放这里可避开自动摆放的等值线标签
fb.field_note(ax, f'影响深度 {depth_mm:.1f} mm',
              (float(t_h[-1]) * 0.52, r_front * 0.42), xytext=(0, 0),
              va='center')

ax.set_xlabel('干燥时间 $t$ (h)')
ax.set_ylabel('径向坐标 $r$ (cm)')
ax.set_xlim(float(t_h[0]), float(t_h[-1]))
ax.set_ylim(0.0, float(r_cm[-1]))

save_fig(fig, fb.out('fig_q2_field_C'))

"""问题一含水率场 C(r,t) 云图（单面板，附录 2 常物性，0-1800 s）。

本图讲什么：与温度场同口径的含水率二维分布。横轴时间、纵轴半径，填充色为
干基含水率；叠加等值线，并用虚线标出末时刻失水影响深度（由场数据反算：
末时刻含水率偏离初值超过 1% 的最内侧半径）。
数据来源：figures/problem_1_results.json → field_C_full (1800×21 全量)、
profile_t_s、file_radii_cm。
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
C = fb.mat(d['field_C_full'])

C0 = float(C[0, 0])
end = C[-1, :]
moved = np.where(np.abs(end - C0) > 0.01 * C0)[0]
r_front = float(r_cm[moved[0]]) if moved.size else float(r_cm[-1])
depth_mm = (float(r_cm[-1]) - r_front) * 10.0

fig, ax = plt.subplots(figsize=(6.4, 3.8), layout='constrained')
fb.set_paper_placement(fig)

levels = np.linspace(float(np.nanmin(C)), float(np.nanmax(C)), 15)
cf = ax.contourf(t_s, r_cm, C.T, levels=levels, cmap=fb.CMAP_C)
cs = ax.contour(t_s, r_cm, C.T, levels=levels[1::3], colors=[COLORS['dark']],
                linewidths=0.45, alpha=0.55)
fb.label_contours(ax, cs, '%.2f')

cb = fig.colorbar(cf, ax=ax, pad=0.02, aspect=26)
cb.set_label('干基含水率 $C$ (kg/kg)')
cb.ax.tick_params(length=2.0)

ax.axhline(r_front, color=COLORS['text'], linewidth=0.9, linestyle='--',
           dashes=(4, 2))
fb.field_note(ax, f'影响深度 {depth_mm:.1f} mm',
              (float(t_s[-1]) * 0.62, r_front), xytext=(0, -7))

ax.set_xlabel('干燥时间 $t$ (s)')
ax.set_ylabel('径向坐标 $r$ (cm)')
ax.set_xlim(float(t_s[0]), float(t_s[-1]))
ax.set_ylim(0.0, float(r_cm[-1]))

save_fig(fig, fb.out('fig_q1_field_C'))

"""问题一热/湿尺度分离图（单面板双轴）。

本图讲什么：同一时程内温度与含水率在中心、表面两处的推进速度差异。
左轴：温度时程（中心与表面几乎重合，说明热扰动早已穿透）；
右轴：含水率时程（表面显著下降、中心仍为初值）——湿分渗透深度远小于半径。
虚线为湿分影响前沿深度（由 field_C_full 逐时刻反算：含水率偏离初值超过 1%
的最内侧半径到表面的距离）。
数据来源：figures/problem_1_results.json → profile_T_5r / profile_C_5r
（1800×5，第 1 列为中心 r=0、末列为表面 r=R）与 field_C_full、file_radii_cm。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS, _lighten

setup_style()

d = fb.load('problem_1_results')
t_s = fb.arr(d['profile_t_s'])
r_cm = fb.arr(d['file_radii_cm'])
T5 = fb.mat(d['profile_T_5r'])
C5 = fb.mat(d['profile_C_5r'])
C = fb.mat(d['field_C_full'])

C0 = float(C[0, 0])
front_mm = np.empty(t_s.size)
for i in range(t_s.size):
    moved = np.where(np.abs(C[i] - C0) > 0.01 * C0)[0]
    r_in = float(r_cm[moved[0]]) if moved.size else float(r_cm[-1])
    front_mm[i] = (float(r_cm[-1]) - r_in) * 10.0

C_T, C_C, C_F = PALETTE[2], PALETTE[4], PALETTE[0]

fig, ax = plt.subplots(figsize=(6.4, 3.9))
fb.set_paper_placement(fig)
ax2 = ax.twinx()

ax.plot(t_s, T5[:, 0], color=C_T, linewidth=1.5, label='中心温度')
ax.plot(t_s, T5[:, -1], color=_lighten(C_T, 0.45), linewidth=1.5,
        linestyle='--', dashes=(5, 2.5), label='表面温度')
ax2.plot(t_s, C5[:, 0], color=C_C, linewidth=1.5, label='中心含水率')
ax2.plot(t_s, C5[:, -1], color=_lighten(C_C, 0.45), linewidth=1.5,
         linestyle='--', dashes=(5, 2.5), label='表面含水率')
ax3 = ax.twinx()
ax3.spines['right'].set_position(('axes', 1.13))
ax3.plot(t_s, front_mm, color=C_F, linewidth=1.3, label='湿分前沿深度')
ax3.set_ylabel('湿分前沿深度 (mm)', color=C_F)
ax3.tick_params(axis='y', colors=C_F)
ax3.spines['right'].set_color(C_F)
ax3.spines['top'].set_visible(False)
ax3.set_ylim(0.0, float(front_mm.max()) * 1.6)
ax3.annotate(f'{front_mm[-1]:.1f} mm', xy=(float(t_s[-1]), float(front_mm[-1])),
             xytext=(-4, 5), textcoords='offset points', ha='right',
             va='bottom', fontsize=8.5, color=C_F)

ax.set_xlabel('干燥时间 $t$ (s)')
ax.set_ylabel('温度 $T$ (°C)', color=C_T)
ax2.set_ylabel('干基含水率 $C$ (kg/kg)', color=C_C)
ax.tick_params(axis='y', colors=C_T)
ax2.tick_params(axis='y', colors=C_C)
ax.set_xlim(float(t_s[0]), float(t_s[-1]))
# 中心含水率是一条贴近量程上沿的常值线，抬高右轴上限把它压下来，给图例腾出净空
C_lo = float(np.nanmin(C5)) - 0.02 * float(np.nanmax(C5))
ax2.set_ylim(C_lo, float(np.nanmax(C5)) * 1.16)
ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.spines['left'].set_color(C_T)
ax2.spines['right'].set_color(C_C)

handles = ax.get_lines()[:2] + ax2.get_lines()[:2] + ax3.get_lines()[:1]
leg = fb.auto_legend(ax, handles=handles,
                     labels=[h.get_label() for h in handles],
                     loc='center right', ncol=1)
if leg is not None:
    leg.get_frame().set_alpha(1.0)
    leg.get_frame().set_facecolor('white')
    leg.set_zorder(20)

save_fig(fig, fb.out('fig_q1_penetration'))

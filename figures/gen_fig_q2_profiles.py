"""问题二径向剖面族（两面板：温度 / 含水率，6 个输出时刻）。

本图讲什么：把表 3、表 4 的读数还原为连续径向剖面，看变物性闭合下的剖面演化。
(a) 温度剖面：3 h 内已接近热风温度，径向温差很小；
(b) 含水率剖面：表层持续下坡而芯部仍在初值附近，湿分前沿缓慢内移。
散点为表 3/表 4 的 5 个输出半径读数，实线为同一时刻的 21 点全剖面。
数据来源：figures/problem_2_results.json → table_times_h、table3_T、table4_C
（散点）与 field_T_full / field_C_full、profile_t_s（曲线）；
径向网格取 problem_1_results.json → file_radii_cm。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS, _lighten

setup_style()

d2 = fb.load('problem_2_results')
t_h = fb.arr(d2['profile_t_s']) / 3600.0
r_cm = fb.arr(fb.load('problem_1_results')['file_radii_cm'])
T = fb.mat(d2['field_T_full'])
C = fb.mat(d2['field_C_full'])
t_tab = fb.arr(d2['table_times_h'])
T_tab = fb.mat(d2['table3_T'])
C_tab = fb.mat(d2['table4_C'])
r_tab = np.linspace(0.0, float(r_cm[-1]), T_tab.shape[1])

idx = [int(np.argmin(np.abs(t_h - tt))) for tt in t_tab]
shades = np.linspace(0.62, 0.0, len(idx))

fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.1), layout='constrained')
fb.set_paper_placement(fig)

for ax, field, tab, lab in ((axes[0], T, T_tab, '温度 $T$ (°C)'),
                            (axes[1], C, C_tab, '干基含水率 $C$ (kg/kg)')):
    for k, (i, sh) in enumerate(zip(idx, shades)):
        ax.plot(r_cm, field[i, :], color=_lighten(PALETTE[4], sh),
                linewidth=1.3, label=f'{t_tab[k]:.1f} h')
        ax.scatter(r_tab, tab[k, :], s=13, color=_lighten(PALETTE[5], sh),
                   zorder=4, edgecolors='white', linewidths=0.5)
    ax.set_xlabel('径向坐标 $r$ (cm)')
    ax.set_ylabel(lab)
    ax.set_xlim(0.0, float(r_cm[-1]))
    fb.declutter_axes(ax, grid=False)

fb.panel(axes[0], '(a) 温度剖面')
fb.panel(axes[1], '(b) 含水率剖面')

fb.shared_legend(fig, list(axes), where='right', ncol=1)

save_fig(fig, fb.out('fig_q2_profiles'))

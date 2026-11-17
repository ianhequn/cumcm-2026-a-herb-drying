"""问题一径向剖面族（两面板：温度 / 含水率，7 个输出时刻）。

本图讲什么：把表 1、表 2 的离散读数还原成连续径向剖面，看两个场的形状差异。
(a) 温度剖面：整根药材几乎同步升温，剖面近乎水平；
(b) 含水率剖面：仅表层下降，芯部保持初值——两者尺度分离的直接证据。
散点为表 1/表 2 的 5 个输出半径处读数，实线为同一时刻的 21 点全剖面。
数据来源：figures/problem_1_results.json → table_times_s、table1_T、table2_C
（散点）与 field_T_full / field_C_full、profile_t_s、file_radii_cm（曲线）。
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
T = fb.mat(d['field_T_full'])
C = fb.mat(d['field_C_full'])
t_tab = fb.arr(d['table_times_s'])
T_tab = fb.mat(d['table1_T'])
C_tab = fb.mat(d['table2_C'])
r_tab = np.linspace(0.0, float(r_cm[-1]), T_tab.shape[1])

idx = [int(np.argmin(np.abs(t_s - tt))) for tt in t_tab]
shades = np.linspace(0.62, 0.0, len(idx))

fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.1), layout='constrained')
fb.set_paper_placement(fig)

for ax, field, tab, lab in ((axes[0], T, T_tab, '温度 $T$ (°C)'),
                            (axes[1], C, C_tab, '干基含水率 $C$ (kg/kg)')):
    for k, (i, sh) in enumerate(zip(idx, shades)):
        ax.plot(r_cm, field[i, :], color=_lighten(PALETTE[2], sh),
                linewidth=1.3, label=f'{t_tab[k]:.0f} s')
        ax.scatter(r_tab, tab[k, :], s=13, color=_lighten(PALETTE[5], sh),
                   zorder=4, edgecolors='white', linewidths=0.5)
    ax.set_xlabel('径向坐标 $r$ (cm)')
    ax.set_ylabel(lab)
    ax.set_xlim(0.0, float(r_cm[-1]))
    fb.declutter_axes(ax, grid=False)

fb.panel(axes[0], '(a) 温度剖面')
fb.panel(axes[1], '(b) 含水率剖面')

fb.shared_legend(fig, list(axes), where='right', ncol=1)

save_fig(fig, fb.out('fig_q1_profiles'))

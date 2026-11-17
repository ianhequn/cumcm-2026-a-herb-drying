"""附件 1 热风工况的双轴时程图（单面板）。

本图讲什么：干燥箱热风温度与含湿量在 0-4 h 内的升温—稳定过程，以及
求解器对 t > 14400 s 的延拓口径（取尾段 24 点均值作常值延拓）。
左轴：热风温度 (°C)；右轴：热风含湿量 (kg/kg)。阴影为尾段 24 点统计窗口，
虚线段为延拓区间。
数据来源：figures/problem_2_results.json → air.t_s / air.T_c / air.W_kgkg
（n=241 全量，无抽样），延拓值取 air.T_tail_mean / air.W_tail_mean。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS, _lighten

setup_style()

air = fb.load('problem_2_results')['air']
t_h = fb.arr(air['t_s']) / 3600.0
T = fb.arr(air['T_c'])
W = fb.arr(air['W_kgkg'])
tail_h = float(air['tail_t_start_s']) / 3600.0
T_tail = float(air['T_tail_mean'])
W_tail = float(air['W_tail_mean'])

t_end = float(t_h[-1])
t_ext = np.array([t_end, t_end + 2.0])

C_T, C_W = PALETTE[2], PALETTE[3]

fig, ax = plt.subplots(figsize=(6.4, 3.9))
fb.set_paper_placement(fig)
ax2 = ax.twinx()

ax.axvspan(tail_h, t_end, color=_lighten(PALETTE[6], 0.80), zorder=0)

ax.plot(t_h, T, color=C_T, linewidth=1.5, label='热风温度')
ax.plot(t_ext, [T_tail, T_tail], color=C_T, linewidth=1.2, linestyle='--',
        dashes=(4, 2), label='温度常值延拓')
ax2.plot(t_h, W, color=C_W, linewidth=1.5, label='热风含湿量')
ax2.plot(t_ext, [W_tail, W_tail], color=C_W, linewidth=1.2, linestyle='--',
         dashes=(4, 2), label='含湿量常值延拓')

ax.annotate(f'{T_tail:.2f} °C', xy=(t_ext[1], T_tail),
            xytext=(-4, -11), textcoords='offset points',
            ha='right', va='top', fontsize=8.5, color=C_T)
ax2.annotate(f'{W_tail * 1000:.2f} g/kg', xy=(t_ext[1], W_tail),
             xytext=(-4, 7), textcoords='offset points',
             ha='right', va='bottom', fontsize=8.5, color=C_W)
n_tail = int(np.sum(t_h >= tail_h))
ax.annotate(f'尾段 {n_tail} 点窗口',
            xy=(0.5 * (tail_h + t_end), float(np.nanmin(T))),
            xytext=(0, 4), textcoords='offset points',
            ha='center', va='bottom', fontsize=8.5, color=COLORS['text'])

ax.set_xlabel('干燥时间 $t$ (h)')
ax.set_ylabel('热风温度 $T_{\\mathrm{a}}$ (°C)', color=C_T)
ax2.set_ylabel('热风含湿量 $W_{\\mathrm{a}}$ (kg/kg)', color=C_W)
ax.tick_params(axis='y', colors=C_T)
ax2.tick_params(axis='y', colors=C_W)
ax.set_xlim(0.0, float(t_ext[1]))
ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.spines['left'].set_color(C_T)
ax2.spines['right'].set_color(C_W)

handles = ax.get_lines()[:2] + ax2.get_lines()[:2]
fb.auto_legend(ax, handles=handles, labels=[h.get_label() for h in handles],
               loc='lower right')

save_fig(fig, fb.out('fig_air_condition'))

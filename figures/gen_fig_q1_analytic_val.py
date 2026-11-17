"""问题一解析解校验（两面板：曲线重叠 / 残差包络）。

本图讲什么：用 Bessel 特征函数展开 + Duhamel 叠加（减法形式 theta=T-T_air(t)）
独立求解同一定解问题，与有限体积数值解逐点对比。
(a) 中心与表面两处的解析（虚线）与数值（实线）温度时程重叠；
(b) 全部 21 个文件半径上的解析-数值偏差包络（阴影 = 各半径偏差的最小-最大范围，
    实线 = 绝对偏差的径向均值），标出全场最大偏差。
数据来源：数值解取 figures/problem_1_results.json → field_T_full (1800×21 全量)；
解析解由 code/problem1.py::analytic_temperature 以 code/utils.py::AirState
（附件 1 原始 241 点）在同一时间/半径网格上重算，模态数取 analytic.mode_convergence
中已收敛的最大档。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, 'code'))

import matplotlib.pyplot as plt
import numpy as np

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS, _lighten
from problem1 import analytic_temperature
from utils import AirState

setup_style()

d = fb.load('problem_1_results')
t_s = fb.arr(d['profile_t_s'])
r_cm = fb.arr(d['file_radii_cm'])
T_num = fb.mat(d['field_T_full'])
n_modes = max(int(k) for k in d['analytic']['mode_convergence'])

T_ana, zeta, lam, _A = analytic_temperature(AirState(), r_cm / 100.0, t_s, n_modes)
dev = T_ana - T_num
abs_dev = np.abs(dev)
max_dev = float(abs_dev.max())

fig, axes = plt.subplots(2, 1, figsize=(6.4, 4.8), sharex=True,
                         layout='constrained', height_ratios=[1.55, 1.0])
fb.set_paper_placement(fig)
ax1, ax2 = axes

for j, name, col in ((0, '中心', PALETTE[2]), (-1, '表面', PALETTE[3])):
    ax1.plot(t_s, T_num[:, j], color=col, linewidth=1.5,
             label=f'{name} 数值解')
    ax1.plot(t_s, T_ana[:, j], color=_lighten(col, 0.05), linewidth=1.1,
             linestyle='--', dashes=(5, 2.5), label=f'{name} 解析解')
ax1.set_ylabel('温度 $T$ (°C)')
fb.declutter_axes(ax1, grid=False)
fb.panel(ax1, '(a) 时程重叠')
fb.auto_legend(ax1, loc='lower right', ncol=2)

fb.uncertainty_band(ax2, t_s, dev.min(axis=1), dev.max(axis=1),
                    color=PALETTE[4], alpha=0.22,
                    label='各半径偏差范围')
ax2.plot(t_s, abs_dev.mean(axis=1), color=PALETTE[5], linewidth=1.3,
         label='绝对偏差径向均值')
ax2.axhline(0.0, color=COLORS['ref_line'], linewidth=0.7, linestyle='--',
            dashes=(4, 2))
i_max = int(np.unravel_index(abs_dev.argmax(), abs_dev.shape)[0])
ax2.annotate(f'{max_dev:.2e} °C', xy=(t_s[i_max], dev[i_max].max()),
             xytext=(-6, 6), textcoords='offset points', ha='right',
             va='bottom', fontsize=8.5, color=PALETTE[5])
ax2.set_xlabel('干燥时间 $t$ (s)')
ax2.set_ylabel('解析$-$数值偏差 (°C)')
ax2.set_xlim(float(t_s[0]), float(t_s[-1]))
fb.declutter_axes(ax2, grid=False)
fb.panel(ax2, f'(b) 残差包络（{n_modes} 模态）')
fb.auto_legend(ax2, loc='lower left', ncol=2)

save_fig(fig, fb.out('fig_q1_analytic_val'))

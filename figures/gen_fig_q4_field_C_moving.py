"""问题四移动边界含水率场 C(r,t) 云图（单面板）。

本图讲什么：同时开启体积收缩与附录 4 变物性后，含水率场如何在一个不断内缩的
计算域上演化。填充色为干基含水率，仅绘制 r ≤ R(t) 的物理存在区域，R(t) 之外
留白；粗实线为收缩边界 R(t) 本身（由 series_R_cm 给出），可见半径在前若干小时
迅速由 R0 收缩并随后趋于 R_min 平台；另一条粗线为 C=0.15 干燥判据等值线，竖虚线
为逐点判据给出的 t_dry。与问题三同类图对照即可看出移动边界把扩散路径缩短。
数据来源：figures/problem_4_results.json → field_C_full（307×21，外部区域在求解
侧已置 None，本脚本转为 nan 后自然留白；时间轴按 field_C_full_stride 从
series_t_h 抽出，全部输出行未再抽样）、series_R_cm / series_t_h、t_dry_h、
radius_history.R0_cm / .R_min_cm；径向网格取同文件 file_radii_cm。
判据阈值取 code/params.py::C_TH。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, 'code'))

import matplotlib.pyplot as plt
import numpy as np

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS
import params as P

setup_style()

d4 = fb.load('problem_4_results')
C = fb.mat(d4['field_C_full'])
stride = int(d4['field_C_full_stride'])
t_h = fb.arr(d4['series_t_h'])[::stride][:C.shape[0]]
r_cm = fb.arr(d4['file_radii_cm'])
t_R = fb.arr(d4['series_t_h'])
R_cm = fb.arr(d4['series_R_cm'])
t_dry = float(d4['t_dry_h'])
C_th = float(P.C_TH)
R0 = float(d4['radius_history']['R0_cm'])
R_min = float(d4['radius_history']['R_min_cm'])

fig, ax = plt.subplots(figsize=(6.4, 3.8), layout='constrained')
fb.set_paper_placement(fig)

levels = np.linspace(float(np.nanmin(C)), float(np.nanmax(C)), 15)
cf = ax.contourf(t_h, r_cm, C.T, levels=levels, cmap=fb.CMAP_C)
cs = ax.contour(t_h, r_cm, C.T, levels=levels[1::3], colors=[COLORS['dark']],
                linewidths=0.45, alpha=0.55)
fb.label_contours(ax, cs, '%.2f')
ax.contour(t_h, r_cm, C.T, levels=[C_th], colors=[COLORS['dark']],
           linewidths=1.6)

ax.plot(t_R, R_cm, color=PALETTE[5], linewidth=2.0, zorder=6)
ax.axvline(t_dry, color=COLORS['dark'], linewidth=1.0, linestyle='--',
           dashes=(4, 2), zorder=6)

cb = fig.colorbar(cf, ax=ax, pad=0.02, aspect=26)
cb.set_label('干基含水率 $C$ (kg/kg)')
cb.ax.tick_params(length=2.0)

fb.field_note(ax, f'$R(t)$  {R0:.2f}→{R_min:.2f} cm',
              (t_dry * 0.42, R_min), xytext=(0, 9), va='bottom')
fb.field_note(ax, f'$t_{{\\mathrm{{dry}}}}$={t_dry:.2f} h',
              (t_dry, R0 * 0.86), xytext=(-5, 0), ha='right', va='center')

ax.set_xlabel('干燥时间 $t$ (h)')
ax.set_ylabel('径向坐标 $r$ (cm)')
ax.set_xlim(float(t_h[0]), float(t_h[-1]))
ax.set_ylim(0.0, float(r_cm[-1]))

save_fig(fig, fb.out('fig_q4_field_C_moving'))

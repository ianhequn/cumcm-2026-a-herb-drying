"""问题三全程含水率场 C(r,t) 云图（单面板，含 0.15 判据等值线）。

本图讲什么：固定半径、附录 4 物性下 57 h 全程的含水率二维分布。横轴时间（h）、
纵轴半径，填充色为干基含水率；粗实线为 C=0.15 干燥判据等值线，竖虚线为逐点
判据给出的干燥完成时刻 t_dry（最后一个越过 0.15 的位置在中心）。
数据来源：figures/problem_3_results.json → field_C_full（345×21，按
field_C_full_stride 从 series_t_h 抽出对应时间轴，全部输出行，未再抽样）、
t_dry_h、C_max_at_end；径向网格取 problem_1_results.json → file_radii_cm。
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

d3 = fb.load('problem_3_results')
C = fb.mat(d3['field_C_full'])
stride = int(d3['field_C_full_stride'])
t_h = fb.arr(d3['series_t_h'])[::stride][:C.shape[0]]
r_cm = fb.arr(fb.load('problem_1_results')['file_radii_cm'])
t_dry = float(d3['t_dry_h'])
C_th = float(P.C_TH)

fig, ax = plt.subplots(figsize=(6.4, 3.8), layout='constrained')
fb.set_paper_placement(fig)

levels = np.linspace(float(np.nanmin(C)), float(np.nanmax(C)), 15)
cf = ax.contourf(t_h, r_cm, C.T, levels=levels, cmap=fb.CMAP_C)
cs = ax.contour(t_h, r_cm, C.T, levels=levels[1::3], colors=[COLORS['dark']],
                linewidths=0.45, alpha=0.55)
fb.label_contours(ax, cs, '%.2f')
ax.contour(t_h, r_cm, C.T, levels=[C_th], colors=[COLORS['dark']],
           linewidths=1.6)

cb = fig.colorbar(cf, ax=ax, pad=0.02, aspect=26)
cb.set_label('干基含水率 $C$ (kg/kg)')
cb.ax.tick_params(length=2.0)

ax.axvline(t_dry, color=COLORS['dark'], linewidth=1.0, linestyle='--',
           dashes=(4, 2))
fb.field_note(ax, f'$t_{{\\mathrm{{dry}}}}$={t_dry:.2f} h',
              (t_dry, float(r_cm[-1]) * 0.55), xytext=(-5, 0), ha='right',
              va='center')
fb.field_note(ax, f'$C$={C_th:.2f} 判据线',
              (t_dry * 0.60, float(r_cm[-1]) * 0.16), xytext=(0, -4))

ax.set_xlabel('干燥时间 $t$ (h)')
ax.set_ylabel('径向坐标 $r$ (cm)')
ax.set_xlim(float(t_h[0]), float(t_h[-1]))
ax.set_ylim(0.0, float(r_cm[-1]))

save_fig(fig, fb.out('fig_q3_field_C'))

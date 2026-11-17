"""问题三每 6 h 径向剖面族（单面板，含 0.15 判据线）。

本图讲什么：表 5 的每 6 h 输出时刻加末时刻，含水率径向剖面如何整体下移并最终
整根低于判据。曲线由浅到深对应时间推进；水平虚线为 C=0.15 干燥判据；末时刻剖面
加粗，其最大值即逐点判据实际卡住的位置。
数据来源：figures/problem_3_results.json → table5_times_h、table5_C（10 行×5 个
输出半径，全部使用）与 t_dry_h、C_max_at_end；判据阈值取 code/params.py::C_TH。
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
import params as P

setup_style()

d3 = fb.load('problem_3_results')
t_tab = fb.arr(d3['table5_times_h'])
C_tab = fb.mat(d3['table5_C'])
r_max = float(fb.arr(fb.load('problem_1_results')['file_radii_cm'])[-1])
r_out = np.linspace(0.0, r_max, C_tab.shape[1])
C_th = float(P.C_TH)
t_dry = float(d3['t_dry_h'])
C_end = float(d3['C_max_at_end'])

fig, ax = plt.subplots(figsize=(6.4, 3.9), layout='constrained')
fb.set_paper_placement(fig)

# 十条曲线只靠同色系明度难分（尤其灰度打印），故明度跨度收窄 + 线型三循环
shades = np.linspace(0.52, 0.0, t_tab.size)
dash_cycle = ['solid', (0, (5, 2)), (0, (1, 1.6))]
for k in range(t_tab.size):
    last = k == t_tab.size - 1
    ax.plot(r_out, C_tab[k], color=_lighten(PALETTE[4], shades[k]),
            linewidth=2.0 if last else 1.2,
            linestyle='solid' if last else dash_cycle[k % 3],
            marker='o' if last else None, markersize=3.4,
            label=f'{t_tab[k]:.1f} h')

ax.axhline(C_th, color=COLORS['dark'], linewidth=1.1, linestyle='--',
           dashes=(5, 2.5), zorder=5)

ax.set_xlabel('径向坐标 $r$ (cm)')
ax.set_ylabel('干基含水率 $C$ (kg/kg)')
ax.set_xlim(0.0, r_max)
ax.set_ylim(-0.16, float(np.nanmax(C_tab)) * 1.06)
ax.annotate(f'判据 {C_th:.2f}', xy=(r_max * 0.5, C_th), xytext=(0, -12),
            textcoords='offset points', ha='center', va='top', fontsize=8.5,
            color=COLORS['dark'])
ax.annotate(f'{t_dry:.2f} h 末刻 {C_end:.4f}', xy=(r_max * 0.5, 0.0),
            xytext=(0, -13), textcoords='offset points', ha='center', va='top',
            fontsize=8.5, color=PALETTE[4])
fb.declutter_axes(ax, grid=False)
fb.auto_legend(ax, loc='lower left', ncol=3)

save_fig(fig, fb.out('fig_q3_profiles_6h'))

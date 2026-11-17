"""问题三干燥动力学（单面板双轴）。

本图讲什么：中心、表面、截面平均三条含水率时程与瞬时干燥速率的关系。
左轴：三条含水率时程 + C=0.15 判据线 + t_dry 竖线；
右轴：平均含水率的下降速率 −dC̄/dt（由 series_Cbar 对 series_t_h 中心差分求得），
用于显示速率峰值出现在早期、后期进入慢速降速段。
数据来源：figures/problem_3_results.json → series_t_h / series_Cmax /
series_Cbar / series_Csurf（3447 点全量，未抽样）、t_dry_h；
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
from _figbase import setup_style, save_fig, PALETTE, COLORS, _lighten
import params as P

setup_style()

d3 = fb.load('problem_3_results')
t_h = fb.arr(d3['series_t_h'])
C_max = fb.arr(d3['series_Cmax'])
C_bar = fb.arr(d3['series_Cbar'])
C_surf = fb.arr(d3['series_Csurf'])
t_dry = float(d3['t_dry_h'])
C_th = float(P.C_TH)

rate = -np.gradient(C_bar, t_h)

C_M, C_B, C_S, C_R = PALETTE[5], PALETTE[2], PALETTE[4], PALETTE[3]

fig, ax = plt.subplots(figsize=(6.4, 3.9))
fb.set_paper_placement(fig)
ax2 = ax.twinx()

ax2.plot(t_h, rate, color=_lighten(C_R, 0.25), linewidth=1.1,
         label='平均含水率下降速率')
ax.plot(t_h, C_max, color=C_M, linewidth=1.5, label='最大含水率')
ax.plot(t_h, C_bar, color=C_B, linewidth=1.4, linestyle='--', dashes=(5, 2.5),
        label='截面平均含水率')
ax.plot(t_h, C_surf, color=C_S, linewidth=1.3, label='表面含水率')

ax.axhline(C_th, color=COLORS['dark'], linewidth=1.0, linestyle=':',
           zorder=5)
ax.axvline(t_dry, color=COLORS['dark'], linewidth=1.0, linestyle='--',
           dashes=(4, 2), zorder=5)
ax.annotate(f'{t_dry:.2f} h', xy=(t_dry, float(np.nanmax(C_max)) * 0.72),
            xytext=(-5, 0), textcoords='offset points', ha='right',
            va='center', fontsize=8.5, color=COLORS['dark'])
_BOX = dict(boxstyle='square,pad=0.16', facecolor='white', edgecolor='none')
ax.annotate(f'判据 {C_th:.2f}', xy=(t_dry * 0.72, C_th), xytext=(0, 5),
            textcoords='offset points', ha='center', va='bottom', fontsize=8.5,
            color=COLORS['dark'], bbox=_BOX, zorder=12)
i_pk = int(np.nanargmax(rate))
ax2.annotate(f'峰值 {rate[i_pk]:.3f} /h', xy=(float(t_h[i_pk]), float(rate[i_pk])),
             xytext=(26, 10), textcoords='offset points', ha='left',
             va='bottom', fontsize=8.5, color=C_R, bbox=_BOX, zorder=12,
             arrowprops=dict(arrowstyle='-', color=C_R, linewidth=0.8,
                             shrinkA=1.0, shrinkB=2.0))

ax.set_xlabel('干燥时间 $t$ (h)')
ax.set_ylabel('干基含水率 $C$ (kg/kg)')
ax2.set_ylabel('$-\\mathrm{d}\\bar{C}/\\mathrm{d}t$ (1/h)', color=C_R)
ax2.tick_params(axis='y', colors=C_R)
ax.set_xlim(0.0, float(t_h[-1]))
ax.set_ylim(0.0, float(np.nanmax(C_max)) * 1.05)
ax2.set_ylim(0.0, float(np.nanmax(rate)) * 1.55)
ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_color(C_R)
ax.set_zorder(ax2.get_zorder() + 1)
ax.patch.set_visible(False)

handles = ax.get_lines()[:3] + ax2.get_lines()[:1]
fb.auto_legend(ax, handles=handles, labels=[h.get_label() for h in handles],
               loc='upper right')

save_fig(fig, fb.out('fig_q3_drying_kinetics'))

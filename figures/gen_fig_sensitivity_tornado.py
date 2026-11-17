"""单因子敏感性龙卷风图（单面板）。

本图讲什么：五个驱动量各自双向扰动后，干燥完成时刻 t_dry 相对基准的百分比变化。
每个驱动量占一行、两根反向水平条，条长为 rel_pct，颜色区分该扰动是延长还是缩短
t_dry；行序按该驱动量的最大绝对影响（span_pct）从大到小排列，故自上而下即敏感度
排序。R(t) 行以问题四为基准（其余各行以问题三为基准），故其纵轴标签已注明。
数据来源：figures/sensitivity_results.json → rows（11 行中 1 行基准 + 10 行扰动，
全部使用）的 .label / .driver / .rel_pct，行序取 ranking，扰动幅度取
perturbations.dT_air_c / .dW_air_kgkg / .coef_scale / .R_scale。
每行 t_dry 均由 meta.recompute = full_simulation_per_row 重算得到，非线性化估计。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

import _figbase as fb
from _figbase import setup_style, save_fig, PALETTE, COLORS, _lighten

setup_style()

d = fb.load('sensitivity_results')
rows = [r for r in d['rows'] if r['driver'] != 'baseline']
pert = d['perturbations']
order = list(d['ranking'])

by_driver = {}
for r in rows:
    key = r['driver'].rsplit('_', 1)[0]
    by_driver.setdefault(key, []).append(r)

# driver 名到 ranking 名、以及纵轴上显示的扰动幅度
KEY = {'R(t)': 'R_scale', 'h_m': 'hm', 'T_air': 'T_air', 'W_air': 'W_air',
       'h': 'h'}
mag = {
    'R(t)': f"R(t) ±{(pert['R_scale'][0] - 1.0) * 100:.0f}%（问题四基准）",
    'h_m': f"h_m ±{(pert['coef_scale'][0] - 1.0) * 100:.0f}%",
    'T_air': f"T_air ±{pert['dT_air_c']:.1f} °C",
    'W_air': f"W_air ±{pert['dW_air_kgkg']:.3f} kg/kg",
    'h': f"h ±{(pert['coef_scale'][0] - 1.0) * 100:.0f}%",
}

y = np.arange(len(order))[::-1]
fig, ax = plt.subplots(figsize=(6.4, 3.5))
fb.set_paper_placement(fig)

C_UP, C_DN = COLORS['up'], COLORS['down']
for yi, name in zip(y, order):
    pair = by_driver[KEY[name]]
    for r in pair:
        v = float(r['rel_pct'])
        col = C_UP if v > 0 else C_DN
        ax.barh(yi, v, height=0.62, color=_lighten(col, 0.30),
                edgecolor=col, linewidth=0.9, zorder=3)
        ax.annotate(f'{v:+.2f}', xy=(v, yi),
                    xytext=(4 if v > 0 else -4, 0),
                    textcoords='offset points', va='center',
                    ha='left' if v > 0 else 'right', fontsize=8.5, color=col,
                    zorder=5)

ax.axvline(0.0, color=COLORS['dark'], linewidth=1.0, zorder=4)
ax.set_yticks(y)
ax.set_yticklabels([mag[n] for n in order])
ax.set_xlabel('$t_{\\mathrm{dry}}$ 相对基准变化 (%)')
lim = max(abs(float(r['rel_pct'])) for r in rows) * 1.34
ax.set_xlim(-lim, lim)
ax.set_ylim(y.min() - 0.62, y.max() + 0.62)
fb.declutter_axes(ax, grid=False)

fb.auto_legend(ax, handles=[
    Patch(facecolor=_lighten(C_UP, 0.30), edgecolor=C_UP, label='延长干燥'),
    Patch(facecolor=_lighten(C_DN, 0.30), edgecolor=C_DN, label='缩短干燥'),
], labels=['延长干燥', '缩短干燥'], loc='lower right')

save_fig(fig, fb.out('fig_sensitivity_tornado'))

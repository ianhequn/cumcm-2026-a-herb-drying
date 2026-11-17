"""附录 3 物性闭合随状态量的变化（四面板：rho / cp / k / D）。

本图讲什么：问题二采用的变物性闭合式在实际经历的状态区间内变化幅度有多大。
(a) 密度 rho(C)；(b) 比热 cp(C)；(c) 导热系数 k(C)——三者只依赖含水率；
(d) 水分扩散系数 D(C,T)——同时依赖含水率与温度，跨多个数量级，故用对数纵轴，
    并按问题二实际温度区间画三条等温线。
灰色竖线为初始含水率与干燥阈值两条判据位置。
数据来源：闭合式取 code/params.py::props_appendix3（逐字调用，不复写系数）；
状态区间由 figures/problem_2_results.json → field_C_full / field_T_full 的实际
最小-最大值给出（不外推到未经历的区间）。
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

d2 = fb.load('problem_2_results')
C_field = fb.mat(d2['field_C_full'])
T_field = fb.mat(d2['field_T_full'])
props = P.props_appendix3()

C_lo, C_hi = float(np.nanmin(C_field)), float(np.nanmax(C_field))
T_lo, T_hi = float(np.nanmin(T_field)), float(np.nanmax(T_field))
C_grid = np.linspace(C_lo, C_hi, 240)
T_levels = [T_lo, 0.5 * (T_lo + T_hi), T_hi]

fig, axes = plt.subplots(2, 2, figsize=(5.7, 4.9), layout='constrained')
fb.set_paper_placement(fig)
(axA, axB), (axC, axD) = axes

scalar_specs = (
    (axA, props.rho, '密度 rho\n(kg/m^3)', '(a) 密度', PALETTE[2]),
    (axB, props.cp, '比热 cp\n(J/(kg·K))', '(b) 比热', PALETTE[4]),
    (axC, props.k, '导热系数 k\n(W/(m·K))', '(c) 导热系数', PALETTE[3]),
)
for ax, fn, lab, tag, col in scalar_specs:
    y = np.asarray(fn(C_grid, T_levels[-1]), dtype=float)
    ax.plot(C_grid, y, color=col, linewidth=1.6)
    ax.fill_between(C_grid, float(y.min()), y, color=_lighten(col, 0.78),
                    zorder=0)
    ax.set_ylabel(lab)
    fb.panel(ax, tag)

for lev, sh in zip(T_levels, np.linspace(0.55, 0.0, len(T_levels))):
    axD.plot(C_grid, np.asarray(props.D(C_grid, lev), dtype=float),
             color=_lighten(PALETTE[5], sh), linewidth=1.5,
             label=f'{lev:.0f} °C')
axD.set_yscale('log')
axD.minorticks_off()
# 中文字体缺 U+2212，mathtext 的 10^-9 会丢指数 → 刻度标签改为纯文本 1e-9 写法
D_all = np.concatenate([np.asarray(props.D(C_grid, lev), dtype=float)
                        for lev in T_levels])
ticks = np.logspace(np.log10(D_all.min()), np.log10(D_all.max()), 4)
axD.set_yticks(ticks)
# 量级提到轴标签里，刻度只留有效数字，否则刻度过宽会把 (d) 的轴标签挤进 (c) 的画区
axD.set_yticklabels([f'{v / 1e-9:.1f}' for v in ticks])
axD.set_ylabel('扩散系数 D\n(1e-9 m^2/s)')
fb.panel(axD, '(d) 扩散系数')
fb.auto_legend(axD, loc='lower right')

for ax in (axA, axB, axC, axD):
    ax.set_xlim(C_lo, C_hi)
    ax.axvline(P.C_TH, color=COLORS['ref_line'], linewidth=0.8,
               linestyle='--', dashes=(4, 2), zorder=1)
    fb.declutter_axes(ax, grid=False)
# 四格共用同一横轴量程，只在底行标注，避免上行 x 标签压住下行 panel 标题
for ax in (axA, axB):
    ax.tick_params(axis='x', labelbottom=False)
for ax in (axC, axD):
    ax.set_xlabel('干基含水率 C (kg/kg)')

axA.annotate(f'C={P.C_TH:.2f}', xy=(P.C_TH, float(props.rho(P.C_TH, T_hi))),
             xytext=(5, 4), textcoords='offset points', ha='left', va='bottom',
             fontsize=8.5, color=COLORS['text'])

save_fig(fig, fb.out('fig_q2_props_evolution'))

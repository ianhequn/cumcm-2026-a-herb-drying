"""数据图公用引导模块（本文件不是出图脚本，不产 PDF）。

职责：
  1. 把工作区根加入 sys.path，使 `from _utils.plot_utils import ...` 在
     `python figures/gen_fig_xxx.py` 与 `python -m` 两种调用下都可用；
  2. 统一调用 setup_style()（读取 CLAUDE.md 里的自定义 8 色配色板）；
  3. 提供 JSON 读取、面板标签、场域网格构造等共用小工具；
  4. 提供由 PALETTE 派生的连续 colormap（不引入调色板外的新色系）。
数据来源：figures/*.json（由 code/ 下的求解脚本产出），本模块不含任何硬编码数值。
"""
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from _utils.plot_utils import (  # noqa: E402,F401
    setup_style, save_fig, set_paper_placement, auto_legend,
    consolidate_shared_legends, shared_legend, uncertainty_band,
    dynamic_limits, declutter_axes, smart_labels, PALETTE, PALETTE_LIGHT,
    COLORS, _lighten,
)

setup_style()

FIG_DIR = 'figures'


def load(name):
    """读取 figures/<name>.json。"""
    with open(os.path.join(_HERE, name + '.json'), encoding='utf-8') as fh:
        return json.load(fh)


def out(name):
    """出图路径（矢量 PDF，供 LaTeX \\includegraphics 用）。"""
    return os.path.join(FIG_DIR, name + '.pdf')


def arr(values):
    """列表→float 数组，JSON 的 null 转 nan（问题四域外已被求解器置 null）。"""
    return np.array([np.nan if v is None else v for v in values], dtype=float)


def mat(rows):
    """二维列表→float 矩阵，null 转 nan。"""
    return np.array([[np.nan if v is None else v for v in row] for row in rows],
                    dtype=float)


def panel(ax, tag):
    """左上角面板标签，如 panel(ax, '(a) 温度场')。"""
    ax.set_title(tag, loc='left', fontsize=9.5, color=COLORS['text'], pad=4.0)
    return ax


def label_contours(ax, cs, fmt):
    """给等值线打读数标签：深色字 + 不透明白底，保证在任意云图背景上都够对比度。"""
    texts = ax.clabel(cs, fmt=fmt, fontsize=8.5, inline=True, inline_spacing=3)
    for txt in texts:
        txt.set_color(COLORS['dark'])
        txt.set_bbox(dict(boxstyle='square,pad=0.10', facecolor='white',
                          edgecolor='none'))
    return texts


def field_note(ax, text, xy, **kwargs):
    """云图上的短锚点标注：深色字 + 不透明白底。"""
    opts = dict(ha='center', va='top', fontsize=8.5, color=COLORS['dark'],
                textcoords='offset points',
                bbox=dict(boxstyle='square,pad=0.16', facecolor='white',
                          edgecolor='none'))
    opts.update(kwargs)
    return ax.annotate(text, xy=xy, **opts)


def _ramp(name, stops):
    return LinearSegmentedColormap.from_list(name, stops, N=256)


# 由自定义配色板派生的顺序 colormap：浅→深单调，灰度打印仍可读。
CMAP_T = _ramp('mh_T', [_lighten(PALETTE[1], 0.72), _lighten(PALETTE[0], 0.30),
                        PALETTE[3], _lighten(PALETTE[5], 0.10)])
CMAP_C = _ramp('mh_C', [_lighten(PALETTE[1], 0.78), _lighten(PALETTE[4], 0.45),
                        PALETTE[4], _lighten(PALETTE[5], 0.05)])

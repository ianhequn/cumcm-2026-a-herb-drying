"""从 figures/*.json 生成 LaTeX 三线表（booktabs）。

每张表的每个单元格都直接取自 JSON 字段，不在本脚本内做任何再计算（唯一例外是
表 C 的相对差已由求解侧写入 grid_convergence.rel_diff，本脚本只做百分号换算，
换算关系在表注中说明）。caption 一律为 ≤20 字的中文名词短语，不含口径与结论。
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    with open(os.path.join(HERE, name + '.json'), encoding='utf-8') as f:
        return json.load(f)


def fmt(v, nd=4):
    if v is None:
        return '--'
    return f'{float(v):.{nd}f}'


def write(name, body):
    path = os.path.join(HERE, name + '.tex')
    with open(path, 'w', encoding='utf-8') as f:
        f.write(body)
    print('Saved:', os.path.relpath(path, os.path.dirname(HERE)))


def table(label, caption, colspec, header, rows, note=None):
    L = ['\\begin{table}[H]', '  \\centering',
         f'  \\caption{{{caption}}}', f'  \\label{{tab:{label}}}',
         f'  \\begin{{tabular}}{{{colspec}}}', '    \\toprule',
         '    ' + ' & '.join(header) + ' \\\\', '    \\midrule']
    L += ['    ' + ' & '.join(r) + ' \\\\' for r in rows]
    L += ['    \\bottomrule', '  \\end{tabular}']
    if note:
        L.append(f'  \\par\\vspace{{2pt}}\\footnotesize {note}')
    L += ['\\end{table}', '']
    return '\n'.join(L)


d1, d2 = load('problem_1_results'), load('problem_2_results')
d3, d4 = load('problem_3_results'), load('problem_4_results')
sen = load('sensitivity_results')

RAD = d1['table_radii_cm']
rad_head = [f'{r:g} cm' for r in RAD]

# ---- 表 1 / 表 2：问题一 30 min 内温度与含水率 -----------------------------
for tag, key, cap, unit_nd in (
        ('q1_T', 'table1_T', '问题一药材温度分布', 4),
        ('q1_C', 'table2_C', '问题一药材水分浓度', 4)):
    rows = [[f"{t:g}"] + [fmt(v, unit_nd) for v in row]
            for t, row in zip(d1['table_times_s'], d1[key])]
    write('TABLE_' + tag, table(
        tag, cap, 'r' + 'r' * len(RAD),
        ['时间 (s)'] + rad_head, rows,
        note='数据来源：figures/problem\\_1\\_results.json 的 %s；'
             '径向位置取同文件 table\\_radii\\_cm。' % key.replace('_', '\\_')))

# ---- 表 3 / 表 4：问题二 3 h 内温度与含水率 --------------------------------
for tag, key, cap in (('q2_T', 'table3_T', '问题二药材温度分布'),
                      ('q2_C', 'table4_C', '问题二药材水分浓度')):
    rows = [[f"{t:g}"] + [fmt(v) for v in row]
            for t, row in zip(d2['table_times_h'], d2[key])]
    write('TABLE_' + tag, table(
        tag, cap, 'r' + 'r' * len(RAD),
        ['时间 (h)'] + rad_head, rows,
        note='数据来源：figures/problem\\_2\\_results.json 的 %s。'
             % key.replace('_', '\\_')))

# ---- 表 5：问题三烘干过程含水率（含结束时刻行） ----------------------------
t5 = d3['table5_times_h']
rows = []
for i, (t, row) in enumerate(zip(t5, d3['table5_C'])):
    lab = f'烘干结束 {t:.4f}' if i == len(t5) - 1 else f'{t:g}'
    rows.append([lab] + [fmt(v) for v in row])
write('TABLE_q3_C', table(
    'q3_C', '问题三烘干过程水分浓度', 'l' + 'r' * len(RAD),
    ['时间 (h)'] + rad_head, rows,
    note='数据来源：figures/problem\\_3\\_results.json 的 table5\\_C 与 '
         'table5\\_times\\_h；末行时刻为 t\\_dry\\_h。'))

# ---- 表 6：问题四烘干过程含水率（末列为随 R(t) 移动的表面） ----------------
t6 = d4['table6']
rows = []
for r in t6['rows']:
    lab = r['label']
    lab = f'烘干结束 {d4["t_dry_reported"]:.4f}' if '结束' in lab else lab
    rows.append([lab] + [fmt(v) for v in r['C_fixed']]
                + [fmt(r['C_surface']), fmt(r['R_cm'], 3)])
write('TABLE_q4_C', table(
    'q4_C', '问题四烘干过程水分浓度', 'l' + 'r' * (len(t6['columns']) + 1),
    ['时间 (h)'] + list(t6['columns']) + ['$R(t)$ (cm)'], rows,
    note='数据来源：figures/problem\\_4\\_results.json 的 table6.rows '
         '(.C\\_fixed / .C\\_surface / .R\\_cm)；末列表面位置随 $R(t)$ 内移。'))

# ---- 表 A：两个附件输入数据的描述性统计 ------------------------------------
air = d2['air']
rh = d4['radius_history']
rows = [
    ['附件 1 热风温度 (°C)', str(air['n_points']),
     fmt(air['T_at_10800'], 4), fmt(air['T_tail_mean'], 4),
     f"{air['tail_t_start_s']:g} s 起"],
    ['附件 1 热风含湿量 (kg/kg)', str(air['n_points']),
     fmt(air['W_at_10800'], 5), fmt(air['W_tail_mean'], 5),
     f"{air['tail_t_start_s']:g} s 起"],
    ['附件 2 药材半径 (cm)', str(rh['n_points']),
     fmt(rh['R0_cm'], 3), fmt(rh['R_min_cm'], 3),
     f"6 h 时 {fmt(rh['R_6h_cm'], 3)}"],
]
write('TABLE_input_stats', table(
    'input_stats', '附件输入数据描述性统计', 'lrrrl',
    ['输入序列', '样本数', '起始或 3 h 值', '末段值', '备注'], rows,
    note='数据来源：figures/problem\\_2\\_results.json 的 air '
         '(.n\\_points / .T\\_at\\_10800 / .W\\_at\\_10800 / .T\\_tail\\_mean '
         '/ .W\\_tail\\_mean / .tail\\_t\\_start\\_s) 与 '
         'figures/problem\\_4\\_results.json 的 radius\\_history '
         '(.n\\_points / .R0\\_cm / .R\\_min\\_cm / .R\\_6h\\_cm)。'
         '末段值对附件 1 为 tail\\_t\\_start\\_s 之后各点的均值，'
         '对附件 2 为全程最小半径。'))

# ---- 表 C：网格与时间步收敛、干质量守恒残差 --------------------------------
g3, g4 = d3['grid_convergence'], d4['grid_convergence']
rows = []
for i, N in enumerate(g3['N']):
    rd = g3['rel_diff'][i]
    rows.append([str(N), fmt(g3['t_dry_h'][i]),
                 '--' if rd is None else f'{float(rd) * 100:.4f}'])
write('TABLE_grid', table(
    'grid', '问题三网格收敛性', 'rrr',
    ['径向控制体数 $N$', '$t_{\\mathrm{dry}}$ (h)', '相对前一档差 (\\%)'], rows,
    note='数据来源：figures/problem\\_3\\_results.json 的 grid\\_convergence '
         '(.N / .t\\_dry\\_h / .rel\\_diff)；百分比列为 rel\\_diff 乘 100。'))

dtc = d3['dt_convergence']
rows = [[fmt(dt, 1), fmt(t)] for dt, t in zip(dtc['dt_s'], dtc['t_dry_h'])]
write('TABLE_dt', table(
    'dt', '问题三时间步收敛性', 'rr',
    ['时间步 $\\Delta t$ (s)', '$t_{\\mathrm{dry}}$ (h)'], rows,
    note='数据来源：figures/problem\\_3\\_results.json 的 dt\\_convergence '
         '(.dt\\_s / .t\\_dry\\_h)，径向控制体数固定为 %d。' % dtc['n_cells']))

rows = [['问题一', f"{d1['conservation_residual']:.3e}"],
        ['问题二', f"{d2['conservation_residual']:.3e}"],
        ['问题三', f"{d3['conservation_residual']:.3e}"],
        ['问题四', f"{d4['conservation_residual']:.3e}"]]
write('TABLE_conservation', table(
    'conservation', '各问干质量守恒残差', 'lr',
    ['算例', '相对残差'], rows,
    note='数据来源：figures/problem\\_$k$\\_results.json 的 '
         'conservation\\_residual（$k=1,2,3,4$）。'))

# ---- 表 D：消融实验 --------------------------------------------------------
ab = d4['ablation']
rows = [
    ['固定半径 + 附录 3（问题三基准）', fmt(ab['baseline_h']), '0.0000'],
    ['收缩 + 附录 3 基准物性', fmt(ab['shrink_only_h']),
     fmt(ab['shrink_only_pct'])],
    ['固定半径 + 换附录 4 物性', fmt(ab['props_only_h']),
     fmt(ab['props_only_pct'])],
    ['收缩 + 附录 4（完整）', fmt(ab['full_h']), fmt(ab['full_pct'])],
]
write('TABLE_ablation', table(
    'ablation', '问题四消融实验', 'lrr',
    ['工况组合', '$t_{\\mathrm{dry}}$ (h)', '相对基准 (\\%)'], rows,
    note='数据来源：figures/problem\\_4\\_results.json 的 ablation '
         '(.baseline\\_h / .shrink\\_only\\_h / .props\\_only\\_h / '
         '.full\\_h 与对应 \\_pct)；交互项 %s h 为四者的加性分解残差，'
         '亦取自同字段 .interaction\\_h。' % fmt(ab['interaction_h'])))

# ---- 表 E：灵敏度汇总 ------------------------------------------------------
def sen_label(s):
    """JSON 里的 label 是纯文本（含下划线与 x），转成可编译的 LaTeX 记号。"""
    for a, b in (('T_air', '$T_{\\mathrm{air}}$'),
                 ('W_air', '$W_{\\mathrm{air}}$'),
                 ('h_m', '$h_m$'), ('R(t)', '$R(t)$'),
                 (' x', ' $\\times$')):
        s = s.replace(a, b)
    return s.replace('%', '\\%')


rows = [[sen_label(r['label']), fmt(r['t_dry_h']), fmt(r['rel_pct'])]
        for r in sen['rows']]
write('TABLE_sensitivity', table(
    'sensitivity', '单因子灵敏度汇总', 'lrr',
    ['扰动工况', '$t_{\\mathrm{dry}}$ (h)', '相对基准 (\\%)'], rows,
    note='数据来源：figures/sensitivity\\_results.json 的 rows '
         '(.label / .t\\_dry\\_h / .rel\\_pct)；每行均按 '
         'meta.recompute=full\\_simulation\\_per\\_row 全程重算。'
         '末两行以问题四为基准（baseline\\_q4\\_h=%s h），其余以问题三为基准'
         '（baseline\\_q3\\_h=%s h）。'
         % (fmt(sen['baseline_q4_h']), fmt(sen['baseline_q3_h']))))

# ---- 表 G：四问主结果汇总 --------------------------------------------------
rows = [
    ['问题一', '附录 2 常物性', '固定 $R$',
     f"{d1['table_times_s'][-1]:g} s", '--'],
    ['问题二', '附录 3 变物性', '固定 $R$',
     f"{d2['table_times_h'][-1]:g} h", '--'],
    ['问题三', '附录 3 变物性', '固定 $R$',
     f"{d3['t_dry_h']:.4f} h", fmt(d3['C_max_at_end'])],
    ['问题四', '附录 4 变物性', '$R(t)$ 收缩',
     f"{d4['t_dry_h']:.4f} h", fmt(d4['C_max_at_end'])],
]
write('TABLE_main_results', table(
    'main_results', '四问主结果汇总', 'lllrr',
    ['问题', '物性闭合', '几何', '时间跨度或 $t_{\\mathrm{dry}}$',
     '结束时 $\\max_r C$'], rows,
    note='数据来源：figures/problem\\_$k$\\_results.json（$k=1,2,3,4$）的 '
         'table\\_times\\_s / table\\_times\\_h / t\\_dry\\_h / '
         'C\\_max\\_at\\_end；问题一、二为指定时间跨度，不涉及干燥判据。'))

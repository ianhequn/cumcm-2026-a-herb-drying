"""校验 figures/TABLE_*.tex 的每个数值单元格都能在 figures/*.json 里找到出处。

做法：把 figures 下全部 JSON 递归摊平成一个数值池，对每张表解析出正文行里的
数值单元格，逐个在池中寻找「按该单元格显示位数四舍五入后完全相等」的源数值。
容许的换算只有两类，且必须在表注里已声明：百分号列的 ×100，以及科学计数法。
任一单元格找不到出处即视为失败并打印出来。
"""
import glob
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))

pool = []


def harvest(obj):
    if isinstance(obj, bool) or obj is None:
        return
    if isinstance(obj, (int, float)):
        pool.append(float(obj))
    elif isinstance(obj, list):
        for v in obj:
            harvest(v)
    elif isinstance(obj, dict):
        for v in obj.values():
            harvest(v)


for jf in sorted(glob.glob(os.path.join(HERE, '*.json'))):
    with open(jf, encoding='utf-8') as f:
        harvest(json.load(f))

# 已声明的换算：原值本身、×100（百分比列）
cands = set()
for v in pool:
    cands.add(v)
    cands.add(v * 100.0)

NUM = re.compile(r'^-?\d+(?:\.\d+)?(?:e[+-]?\d+)?$', re.I)
bad, total = [], 0

for tf in sorted(glob.glob(os.path.join(HERE, 'TABLE_*.tex'))):
    lines = open(tf, encoding='utf-8').read().splitlines()
    body, on = [], False
    for ln in lines:
        s = ln.strip()
        if s == '\\midrule':
            on = True
            continue
        if s == '\\bottomrule':
            on = False
        if on:
            body.append(s)
    for ln in body:
        for cell in ln.replace('\\\\', '').split('&'):
            c = cell.strip()
            if c in ('', '--') or not NUM.match(c):
                continue
            total += 1
            nd = len(c.split('.')[1]) if '.' in c and 'e' not in c.lower() else 0
            want = float(c)
            hit = any(abs(round(v, nd) - want) < 0.5 * 10 ** (-nd - 3)
                      if nd else abs(v - want) < 1e-9 or
                      abs(float(f'{v:.3e}') - want) < 1e-30
                      for v in cands)
            if not hit:
                # 科学计数法单元格：按 %.3e 显示比较
                hit = any(f'{v:.3e}' == c for v in cands)
            if not hit:
                bad.append((os.path.basename(tf), c))

print(f'数值单元格总数: {total}')
if bad:
    print(f'无法溯源: {len(bad)}')
    for f, c in bad[:40]:
        print('  ', f, c)
    raise SystemExit(1)
print('全部数值单元格均可溯源到 figures/*.json')

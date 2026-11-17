"""公共求解框架：干质量拉格朗日坐标 s=(r/R(t))^2 上的节点中心有限体积 + 隐式欧拉。

四问共用本模块：问题一/二/三取 R(t)=R0，问题四把 R(t) 换成附件 2 的插值。
全部物理常数来自 params.py，本文件不出现裸物性字面量。
"""
from __future__ import annotations

import json
import os
import random
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import numpy as np
import pandas as pd
from scipy.linalg import solve_banded

from params import *  # noqa: F401,F403  口径常数与附录 2/3/4 物性闭合


def set_all_seeds(seed: int = SEED) -> None:
    """本题为确定性 PDE 求解，无随机采样；固定种子仅为可复现性留痕。"""
    random.seed(seed)
    np.random.seed(seed)


def read_all_sheets(path: str) -> dict:
    """Excel 一律读全部 sheet（不写 sheet_name 时 pandas 只读首表且不告警）。"""
    return pd.read_excel(path, sheet_name=None)


def check_ingest(fname: str, total_rows: int) -> None:
    """与 DATA_PROFILE.json 的行数账本对账，防"只读到首表"式摄入不全。"""
    with open(PROFILE_PATH, encoding='utf-8') as fh:
        profile = json.load(fh)['files']
    prof = profile.get(fname, {})
    assert 'total_rows' in prof, f"[未建档] {fname} 不在 DATA_PROFILE.json"
    assert total_rows == prof['total_rows'], (
        f"[摄入不全] {fname} 实读 {total_rows} 行 != 建档 {prof['total_rows']} 行")


class AirState:
    """烘房边界驱动：附件 1 全部点线性插值 + 末段平台均值常数延拓（假设 6）。

    超出附件 1 覆盖范围时用 np.interp 的 right= 参数返回平台均值，绝不线性外推
    （线性外推会让温度无界上升，见建模报告的延拓讨论）。
    """

    def __init__(self, dT: float = 0.0, dW: float = 0.0, hold_last_value: bool = False):
        sheets = read_all_sheets(A1_PATH)
        frame = sheets['Sheet1']
        check_ingest(os.path.basename(A1_PATH), sum(len(d) for d in sheets.values()))
        self.t = frame['时间'].to_numpy(float)
        self.T = frame['温度'].to_numpy(float) + dT
        self.W = frame['水分浓度'].to_numpy(float) + dW
        self.n_points = int(self.t.size)
        # 平台口径固定为末 24 点（建模报告假设 6）；下面三行按此口径逐字取尾段
        assert AIR_TAIL_POINTS == 24, '平台均值口径应为末 24 点'
        self.T_tail_mean = float(self.T[-24:].mean())
        self.W_tail_mean = float(self.W[-24:].mean())
        self.tail_t_start = float(self.t[-AIR_TAIL_POINTS])
        self.T_hold = float(self.T[-1]) if hold_last_value else self.T_tail_mean
        self.W_hold = float(self.W[-1]) if hold_last_value else self.W_tail_mean

    def T_of(self, t):
        return np.interp(t, self.t, self.T, left=self.T[0], right=self.T_hold)

    def W_of(self, t):
        return np.interp(t, self.t, self.W, left=self.W[0], right=self.W_hold)


class ConstantAir:
    """结构性验证用的常值烘房条件（零输入退化、常空气解析退化）。"""

    def __init__(self, T_const: float, W_const: float):
        self.T_const, self.W_const = float(T_const), float(W_const)
        self.n_points = 0
        self.T_tail_mean, self.W_tail_mean = self.T_const, self.W_const

    def T_of(self, t):
        return np.full_like(np.asarray(t, dtype=float), self.T_const) if np.ndim(t) else self.T_const

    def W_of(self, t):
        return np.full_like(np.asarray(t, dtype=float), self.W_const) if np.ndim(t) else self.W_const


class RadiusHistory:
    """药材半径 R(t)：由附件 2 全部实测点线性插值（禁止自拟收缩律）。

    t 超出附件 2 覆盖范围时用末点值保持（right=），末点即先验上界时刻的实测半径。
    """

    def __init__(self, scale: float = 1.0):
        sheets = read_all_sheets(A2_PATH)
        frame = sheets['Sheet1']
        check_ingest(os.path.basename(A2_PATH), sum(len(d) for d in sheets.values()))
        self.t = frame['时间'].to_numpy(float)
        self.R = frame['半径'].to_numpy(float) * (1.0 / CM_PER_M) * float(scale)
        self.n_points = int(self.t.size)

    def R_of(self, t):
        return np.interp(t, self.t, self.R, left=self.R[0], right=self.R[-1])


class FixedRadius:
    """固定几何 R(t)=R0（问题一/二/三，以及消融实验的"关闭收缩"组）。"""

    def __init__(self, R_const: float = R0_M):
        self.R_const = float(R_const)
        self.n_points = 0

    def R_of(self, t):
        return self.R_const


def s_grid(n_cells: int):
    """s in [0,1] 上的均匀节点与节点中心控制体测度 m_i（端点为半格）。"""
    s = np.linspace(0.0, 1.0, n_cells + 1)
    ds = 1.0 / n_cells
    m = np.full(n_cells + 1, ds)
    m[0] = 0.5 * ds
    m[-1] = 0.5 * ds
    return s, ds, m


def jacobian_coef(R: float) -> float:
    """s 坐标下的时变几何系数 4/R(t)^2（收缩的全部效应集中在此，逐步重算）。"""
    return 4.0 / R ** 2                          # S_COORD_JACOBIAN，逐字保留分子形式


def s_of_radius(r, R):
    """干质量拉格朗日坐标 s = (r / R)**2，裁剪到 [0,1]（r>R 由调用方先剔除）。"""
    return np.clip((r / R) ** 2, 0.0, 1.0)


def _banded_diffusion(diff_nodal, s, ds, m, cap, dt, R):
    """装配 (M/dt + A) 的三对角带状矩阵与 M/dt*x^n 部分，A 为 s 坐标径向通量算子。

    界面通量 flux_{i+1/2} = s_{i+1/2} * D_{i+1/2} * (x_{i+1}-x_i)/ds，
    界面扩散系数取算术平均，故内部通量装配严格反对称 → 逐步守恒。
    """
    n = s.size
    flux_coef = jacobian_coef(R) * 0.5 * (diff_nodal[:-1] + diff_nodal[1:]) \
        * (s[:-1] + 0.5 * ds) / ds          # 各界面 i+1/2 的通量系数
    lower = np.zeros(n)
    upper = np.zeros(n)
    d = cap * m / dt
    d[:-1] += flux_coef
    d[1:] += flux_coef
    upper[1:] = -flux_coef                   # 存 A[i, i+1]，带状矩阵上对角
    lower[:-1] = -flux_coef                  # 存 A[i+1, i]，带状矩阵下对角
    return d, lower, upper


def _solve_tridiag(d, lower, upper, rhs):
    ab = np.zeros((3, d.size))
    ab[0, 1:] = upper[1:]
    ab[1, :] = d
    ab[2, :-1] = lower[:-1]
    return solve_banded((1, 1), ab, rhs, overwrite_ab=True, overwrite_b=True,
                        check_finite=False)


def step_moisture(C_old, C_it, T_it, s, ds, m, dt, R, W_air, props, hm=HM):
    """水分方程一次隐式装配求解；Robin 在末节点精确施加（节点中心格式）。"""
    D_nodal = props.D(C_it, T_it)
    cap = np.ones_like(C_old)
    d, lower, upper = _banded_diffusion(D_nodal, s, ds, m, cap, dt, R)
    rhs = m * C_old / dt
    d[-1] += 2.0 * hm / R
    rhs[-1] += 2.0 * hm / R * W_air
    return _solve_tridiag(d, lower, upper, rhs)


def step_energy(T_old, C_new, T_it, s, ds, m, dt, R, T_air, props, h=H_CONV):
    """能量方程一次隐式装配求解；Robin 在末节点精确施加。"""
    k_nodal = props.k(C_new, T_it)
    cap = props.rho(C_new, T_it) * props.cp(C_new, T_it)
    d, lower, upper = _banded_diffusion(k_nodal, s, ds, m, cap, dt, R)
    rhs = cap * m * T_old / dt
    d[-1] += 2.0 * h / R
    rhs[-1] += 2.0 * h / R * T_air
    return _solve_tridiag(d, lower, upper, rhs)


def solve_transport(props, air, radius, n_cells, dt, t_end, out_dt,
                    hm=HM, h=H_CONV, C_th=None, nsweep=NSWEEP,
                    T_init=T0_C, C_init=C0_KGKG, store_fields=True):
    """隐式向后欧拉 + 每步 nsweep 次顺序 Picard 迭代（先水分、再用新 C 装配能量）。

    返回 dict：s 网格、输出时刻、T/C 场（在 s 节点上）、当时 R(t)、Cbar、Cmax、
    守恒残差、径向最大正梯度，以及（给定 C_th 时）阈值首次到达时刻 t_dry。
    """
    s, ds, m = s_grid(n_cells)
    C = np.full(s.size, float(C_init))
    T = np.full(s.size, float(T_init))
    n_sub = int(round(out_dt / dt))
    assert abs(n_sub * dt - out_dt) < 1e-12 and n_sub >= 1, '输出间隔须为 dt 的整数倍'

    t_now = 0.0
    out_t = [0.0]                              # 输出时刻按 k*out_dt 精确给出，不累加 dt
    out_R = [float(radius.R_of(0.0))]
    out_C = [C.copy()]
    out_T = [T.copy()]
    Cbar = [float(m @ C)]
    Cmax = [float(np.max(C))]
    cons_res = 0.0
    max_pos_grad = float(np.max(np.diff(C)))
    t_dry = None
    n_steps = 0
    max_out = int(np.ceil(t_end / out_dt))

    for k_out in range(max_out):
        t_ref = k_out * out_dt
        for j_sub in range(n_sub):
            t_new = t_ref + (j_sub + 1) * dt   # 相对本输出窗起点计时，避免长时程累加漂移
            R_new = float(radius.R_of(t_new))
            T_air = float(air.T_of(t_new))
            W_air = float(air.W_of(t_new))
            C_old, T_old = C, T
            C_it, T_it = C, T
            for _ in range(nsweep):
                C_it = step_moisture(C_old, C_it, T_it, s, ds, m, dt, R_new, W_air, props, hm=hm)
                T_it = step_energy(T_old, C_it, T_it, s, ds, m, dt, R_new, T_air, props, h=h)
            Cbar_old = float(m @ C_old)
            C, T = C_it, T_it
            Cbar_new = float(m @ C)
            surface_flux = -(2.0 * hm / R_new) * (C[-1] - W_air)
            cons_res = max(cons_res, abs((Cbar_new - Cbar_old) / dt - surface_flux)
                           / max(abs(surface_flux), 1e-30))
            t_now = t_new
            n_steps += 1
        t_now = (k_out + 1) * out_dt
        out_t.append(t_now)
        out_R.append(R_new)
        Cbar.append(float(m @ C))
        Cmax.append(float(np.max(C)))
        max_pos_grad = max(max_pos_grad, float(np.max(np.diff(C))))
        if store_fields:
            out_C.append(C.copy())
            out_T.append(T.copy())
        if C_th is not None and Cmax[-1] < C_th:
            prev = Cmax[IDX_PREV_STEP]
            t_dry = out_t[-1] - out_dt * (C_th - Cmax[-1]) / (prev - Cmax[-1])
            break

    return {
        's': s, 'ds': ds, 'm': m,
        't': np.asarray(out_t), 'R': np.asarray(out_R),
        'C': np.asarray(out_C) if store_fields else None,
        'T': np.asarray(out_T) if store_fields else None,
        'C_final': C, 'T_final': T,
        'Cbar': np.asarray(Cbar), 'Cmax': np.asarray(Cmax),
        'cons_res': float(cons_res), 'max_pos_grad': float(max_pos_grad),
        't_dry': t_dry, 'n_steps': n_steps, 'n_cells': n_cells, 'dt': dt,
    }


def t_dry_hours(props, air, radius, n_cells, dt, hm=HM, h=H_CONV,
                out_dt=float(Q34_FILE_DT_S), t_end_h=DURATION_PRIOR_H[1] * 3.0):
    """只求阈值到达时刻（不存场），用于网格无关性/灵敏度/消融的真实重算。"""
    res = solve_transport(props, air, radius, n_cells, dt, t_end_h * SEC_PER_HOUR,
                          out_dt, hm=hm, h=h, C_th=C_TH, store_fields=False)
    assert res['t_dry'] is not None, '阈值未在给定时窗内达成'
    return res['t_dry'] / SEC_PER_HOUR, res


def field_at_radii(res, radii_cm, key='C'):
    """把 s 网格上的场取到物理半径列：s_q = clip((r/R(t))^2, 0, 1) 后线性插值。

    生产网格（N 为制表半径公倍数）下五个制表半径恰落在节点上，插值退化为取值。
    """
    field = res[key]
    s = res['s']
    radii_m = np.asarray(radii_cm, dtype=float) * (1.0 / CM_PER_M)
    out = np.full((field.shape[0], radii_m.size), np.nan)
    if float(np.ptp(res['R'])) == 0.0:        # 固定几何：s_q 与时间无关，按列向量化
        R_c = float(res['R'][0])
        for j, r_q in enumerate(radii_m):
            if r_q > R_c + 1e-15:
                continue
            s_q = float(s_of_radius(r_q, R_c))
            j_hi = int(np.searchsorted(s, s_q))
            j_hi = min(max(j_hi, 1), s.size - 1)
            w = (s_q - s[j_hi - 1]) / (s[j_hi] - s[j_hi - 1])
            out[:, j] = (1.0 - w) * field[:, j_hi - 1] + w * field[:, j_hi]
        return out
    for i in range(field.shape[0]):
        R_i = res['R'][i]
        for j, r_q in enumerate(radii_m):
            if r_q > R_i + 1e-15:
                continue                      # 超出当时药材，留空，不做域外插值
            s_q = float(s_of_radius(r_q, R_i))
            out[i, j] = np.interp(s_q, s, field[i])
    return out


def surface_value(res, key='C'):
    return res[key][:, -1]


def snapshot_at_time(res, t_query, key='C'):
    """在两个输出步之间对场做线性插值（用于烘干结束这一非整步时刻）。"""
    t = res['t']
    idx = int(np.searchsorted(t, t_query))
    idx = min(max(idx, 1), t.size - 1)
    w = (t_query - t[idx - 1]) / (t[idx] - t[idx - 1])
    field = res[key]
    return (1.0 - w) * field[idx - 1] + w * field[idx], (1.0 - w) * res['R'][idx - 1] + w * res['R'][idx]


NUM_FMT = '0.' + '0' * DECIMALS


def r4(x):
    """四位小数（题面：所有结果保留四位小数）。NaN 保持为空。"""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    return float(np.round(float(x), DECIMALS))


def _fill_sheet(ws, header, times, matrix):
    ws.append(list(header))
    for t_val, row in zip(times, matrix):
        ws.append([t_val] + [r4(v) for v in row])
    for row in ws.iter_rows(min_row=2, min_col=2):
        for cell in row:
            cell.number_format = NUM_FMT


def write_result_xlsx(path, sheets):
    """写 result*.xlsx：sheets = [(sheet_name, header, times, matrix), ...]。

    header 首格为模板的"时间\\到药材中心的距离"，其后为距离标签（或"药材表面"）。
    """
    from openpyxl import Workbook
    wb = Workbook()
    wb.remove(wb.active)
    for name, header, times, matrix in sheets:
        ws = wb.create_sheet(title=name)
        _fill_sheet(ws, header, times, matrix)
    wb.save(path)
    return path


def file_radii_cm(r_max_cm):
    """0.1 cm 步长的输出列（不超过 r_max_cm），四位小数标签。"""
    n = int(np.floor(r_max_cm / FILE_DR_CM + 1e-9))
    return [round(i * FILE_DR_CM, 4) for i in range(n + 1)]


def dump_json(path, payload):
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=1, allow_nan=False)
    return path


def to_list(a, nd=None):
    """ndarray → 可 JSON 序列化的列表；nd 给定时按 nd 位小数落盘（表格用）。"""
    arr = np.asarray(a, dtype=float)
    if nd is not None:
        arr = np.round(arr, nd)
    out = arr.tolist()
    return [None if (isinstance(v, float) and np.isnan(v)) else v for v in out] \
        if arr.ndim == 1 else [[None if (isinstance(v, float) and np.isnan(v)) else v for v in row]
                               for row in out]


def rel_err(a, b):
    return abs(a - b) / max(abs(b), 1e-30)


def echo(msg):
    print(msg, flush=True)




"""参数口径表：全部取值来自 PROBLEM_FACTS.json，import 时逐项断言。

本文件是全工程唯一允许出现数值字面量的位置（附录 2/3/4 的经验公式按题面逐字实现），
其余 code/*.py 必须 `from params import *` 使用具名常量，禁止裸数字。
文件末尾的口径一致性断言块把每个字面量与 PROBLEM_FACTS.json 对回去，
防止长上下文下的记忆漂移。
"""
from __future__ import annotations

import json
import os
from collections import namedtuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(_HERE)
FACTS_PATH = os.path.join(ROOT, 'PROBLEM_FACTS.json')
PROFILE_PATH = os.path.join(ROOT, 'DATA_PROFILE.json')
with open(FACTS_PATH, encoding='utf-8') as _fh:
    FACTS = json.load(_fh)

# ---------------- 几何与初值（题面问题 1）----------------
R0_M = 0.02
LENGTH_M = 0.25
T0_C = 28.0
C0_KGKG = 2.55
KELVIN_OFFSET = 273.15

# ---------------- 边界层系数（附录 2，四问共用）----------------
H_CONV = 25.0
HM = 8e-7

# ---------------- 附录 2：常物性（问题 1）----------------
RHO_A2 = 820.0
CP_A2 = 2600.0
K_A2 = 0.36

# ---------------- 附录 3/4 的活化项（共用指数分母系数）----------------
D_ACT_OVER_T = 3850.0

# ---------------- 数值与输出口径 ----------------
C_TH = 0.15
DECIMALS = 4
DURATION_PRIOR_H = (48.0, 72.0)
C_FLOOR = 1e-12
NSWEEP = 3
N_PROD = 160
N_CONV = (48, 80, 112, 160, 224)
N_ABL = 112
DT_Q1 = 0.25
DT_Q2 = 0.5
DT_Q3 = 1.0
DT_Q4 = 1.0
DT_CONV = (4.0, 2.0, 1.0, 0.5)
AIR_TAIL_POINTS = 24
S_COORD_JACOBIAN = 4.0
Q1_TABLE_TIMES_S = (100, 300, 600, 900, 1200, 1500, 1800)
Q1_TABLE_RADII_CM = (0.0, 0.5, 1.0, 1.5, 2.0)
Q1_HORIZON_S = 1800
Q2_HORIZON_S = 10800
Q2_TABLE_DT_H = 0.5
Q34_TABLE_DT_H = 6.0
Q34_TABLE_DR_CM = 0.5
FILE_DR_CM = 0.1
Q1_FILE_DT_S = 1
Q2_FILE_DT_S = 1
Q34_FILE_DT_S = 60
Q4_TABLE_RADII_CM = (0.0, 0.5, 1.0)
N_CONV_Q4 = (48, 112, 160)

# ---------------- 灵敏度扰动幅度（建模报告 §10 表 E）----------------
SENS_DT_AIR = 0.5
SENS_DW_AIR = 0.004
SENS_SCALE_UP = 1.2
SENS_SCALE_DOWN = 0.8
SENS_R_SCALE_UP = 1.05
SENS_R_SCALE_DOWN = 0.95
SENS_MIN_ROWS = 5                    # 灵敏度表最少行数（含基准）
SENS_TOP_DRIVER = 'R(t)'             # 预期最敏感因子
SENS_BOTTOM_DRIVER = 'h'             # 预期最不敏感因子
DIFFUSION_SCALING_EXPONENT = 2.0     # t ~ R^2/D 的理论标度指数

# ---------------- 建模选择开关（建模报告 §2 假设的参数化名，逐一落地为具名常量）----------------
GEOMETRY_MODE = '1d_radial'                   # 备选 '2d_axisym'：加两端 Robin，会加快失水
MASS_MODEL = 'fickian_dry_basis'              # 干基 Fick 扩散，C 为干基含水率
SHRINKAGE_KINEMATICS = 'attachment2_interp'   # 备选 'volume_additive'：见 §8.4 的一致性核查
SHRINKAGE_SOURCE = 'attachment2_interp'
LATENT_HEAT_COUPLING = 'off'                  # 忽略汽化潜热表面降温；灵敏度显示 h 影响 <0.04%
AIR_EXTRAPOLATION = 'last24_mean_hold'        # 备选 'last_value_hold'：差 0.1555 °C

# ---------------- 验收门限（建模报告 §9/§11）----------------
GRID_REL_TOL = 5e-4
CONS_RES_TOL = 1e-8
DEGEN_TOL = 1e-10
CONSERVE_DRIFT_TOL = 1e-12
CROSS_PROBLEM_TOL = 1e-6
ABL_PROPS_ONLY_MIN_PCT = 80.0
ABL_SHRINK_ONLY_MAX_PCT = -40.0
VOL_INV_MIN_DEV = 0.50
CONST_AIR_T = 50.0

# 体积可加性对照（§8.4）用的纯水密度；仅用于一致性核查，不进入主输运求解
RHO_WATER = 1000.0

Q4_SURFACE_LABEL = FACTS['output_spec']['q4_last_column_label']
HEADER_CORNER = '时间\\到药材中心的距离'
SHEET_T = '温度'
SHEET_C = '水分浓度'
A1_PATH = os.path.join(ROOT, 'user_data', '附件1.xlsx')
A2_PATH = os.path.join(ROOT, 'user_data', '附件2.xlsx')
FIG_DIR = os.path.join(ROOT, 'figures')
TMP_DIR = os.path.join(ROOT, '_tmp')
CM_PER_M = 100.0
SEC_PER_HOUR = 3600.0
PCT = 100.0            # 分数→百分数；与 CM_PER_M 数值相同但量纲无关，不可互用
HOURS_PER_DAY = 24.0
BYTES_PER_KB = 1024.0
HALF = 0.5             # 半幅（峰峰值→±幅度）；与 DT_Q2 数值相同但语义无关，不可互用
SEED = 20260912
RULE_WIDTH = 64        # 控制台分隔线宽度
LOG_NAME_WIDTH = 24    # [main] 步骤名对齐宽度
IDX_PREV_STEP = -2     # 前一输出步的负索引（阈值首次穿越的"未达标"一侧）
DESC_WIDTH = 40        # 审计表判据描述列宽
DUP_DECIMALS = 9       # 判"数值重复（疑非真实重算）"时的比较位数
N_ECHO_HEAD = 3        # 控制台回显的表头行数
N_ECHO_TAIL = 2        # 控制台回显的表尾行数
SERIES_DECIMALS = 6    # figures/*.json 里曲线序列的保留位数（表格与 xlsx 仍用 DECIMALS）
JSON_STRIDE = 10       # 全场落 JSON 的时间下采样步长（仅控体积，xlsx 为全量）

# ---------------- 数据行数账本（取自 DATA_PROFILE.json，禁止凭记忆写点数）----------------
with open(PROFILE_PATH, encoding='utf-8') as _ph:
    PROFILE = json.load(_ph)['files']
N_A1_POINTS = PROFILE['附件1.xlsx']['total_rows']
N_A2_POINTS = PROFILE['附件2.xlsx']['total_rows']

# ---------------- 报告检查点（建模报告 §6/§9/§11 明示的核对值，由运行时断言对回实测）----------------
N_FILE_RADII = int(round(R0_M * CM_PER_M / FILE_DR_CM)) + 1   # 0.0..2.0 cm 共 21 列
AIR_T_PLATEAU_C = 50.0095            # 附件 1 末 24 点温度均值
AIR_W_PLATEAU = 0.04999              # 附件 1 末 24 点含湿均值
AIR_T_LINEXTRAP_MAX_C = 50.2         # 线性外推会越过此值；常数延拓必须低于它
AIR_PROBE_LATE_S = 20000.0           # 超出附件 1 覆盖范围的延拓探针时刻
K_A3_AT_C0_RANGE = (0.47, 0.49)      # k_3(2.55) 的合法区间：越界说明湿基组合被误写为 C
R_AT_6H_MAX_M = 0.0138               # 附件 2 的 R(6 h) 上界
MEAN_CRIT_PROBE_H = 36.0             # 均值判据对照时刻
VOL_INV_PROBE_H = 18.0               # 体积反演偏离核查时刻
N_VI_INTERVALS = int(round(DURATION_PRIOR_H[1] / Q34_TABLE_DT_H))   # 72 h / 6 h = 12 段
ZETA_REF = (1.41847, 4.16647, 7.20848, 10.30827, 13.42716, 16.55454)
N_MODES_CONV = (6, 12, 20, 40)       # 解析级数模态数收敛序列
GRID_N_ALT = (40, 60, 80, 120, 160)  # 补充网格集（与 N_CONV 互为独立佐证）
LAMBDA1_HORIZON_S = 6000.0           # 常空气退化检验时窗
LAMBDA1_FIT_WINDOW_S = (2000.0, 6000.0)   # 高阶模衰减到 1e-6 以下后的拟合窗
DEGEN_HORIZON_S = 600.0              # 零输入退化检验时窗（恒定解，无需长时窗）
BRENT_XTOL = 1e-14
BRENT_RTOL = 1e-15
BRENT_MAXITER = 200

# ---------------- 判据容差（建模报告 §9/§11 的数值门限）----------------
ANALYTIC_DEV_TOL = 1e-3      # 解析-数值最大偏差
ZETA_DEV_TOL = 1e-4          # 特征根偏差
C_CENTER_TOL = 1e-3          # 1800 s 中心含水率仍为 C0
PLATEAU_REL_TOL = 1e-4       # 平台温度均值相对偏差
PLATEAU_W_REL_TOL = 1e-3     # 平台含湿均值相对偏差
GRAD_POS_TOL = 1e-9          # 径向正梯度（单调性）判据
LAMBDA_FIT_REL_TOL = 1e-3    # 拟合衰减率对 lambda_1 的相对偏差
EXACT_TOL = 1e-12            # 浮点严格相等判据（延拓一致、单调不增、报告值回代）

# ---------------- 物性闭合：附录 2/3/4 逐字实现 ----------------
Props = namedtuple('Props', 'name rho cp k D')


def _kelvin(T_celsius):
    return np.asarray(T_celsius, dtype=float) + KELVIN_OFFSET


def _arrhenius(T_celsius):
    """exp(-3850 / T_kelvin)，T 以开尔文计（附录 3/4 明示单位为 K）。"""
    return np.exp(-D_ACT_OVER_T / _kelvin(T_celsius))


def _safe_C(C):
    return np.maximum(np.asarray(C, dtype=float), C_FLOOR)


def props_appendix2():
    """附录 2：rho/cp/k 为常数，D = 7e-9 * exp(-0.89/C)（不含温度项）。"""
    def rho(C, T):
        return np.full_like(np.asarray(C, dtype=float), RHO_A2)

    def cp(C, T):
        return np.full_like(np.asarray(C, dtype=float), CP_A2)

    def k(C, T):
        return np.full_like(np.asarray(C, dtype=float), K_A2)

    def D(C, T):
        return 7e-9 * np.exp(-0.89 / _safe_C(C))

    return Props('appendix2', rho, cp, k, D)


def props_appendix3():
    """附录 3：rho=650+128C, cp=1450+2736*C/(C+1), k=0.21+0.38*C/(C+1),
    D = 2.4e-3 * exp(-0.45/C) * exp(-3850/T_kelvin)。"""
    def rho(C, T):
        return 650.0 + 128.0 * np.asarray(C, dtype=float)

    def cp(C, T):
        Cs = np.asarray(C, dtype=float)
        return 1450.0 + 2736.0 * (Cs / (Cs + 1.0))

    def k(C, T):
        Cs = np.asarray(C, dtype=float)
        return 0.21 + 0.38 * (Cs / (Cs + 1.0))

    def D(C, T):
        return 2.4e-3 * np.exp(-0.45 / _safe_C(C)) * _arrhenius(T)

    return Props('appendix3', rho, cp, k, D)


def props_appendix4():
    """附录 4：rho=760+90C, cp=1850+2150*C/(C+1), k=0.12+0.20*C/(C+1),
    D = 4.2e-4 * exp(-0.30/C) * exp(-3850/T_kelvin)。"""
    def rho(C, T):
        return 760.0 + 90.0 * np.asarray(C, dtype=float)

    def cp(C, T):
        Cs = np.asarray(C, dtype=float)
        return 1850.0 + 2150.0 * (Cs / (Cs + 1.0))

    def k(C, T):
        Cs = np.asarray(C, dtype=float)
        return 0.12 + 0.20 * (Cs / (Cs + 1.0))

    def D(C, T):
        return 4.2e-4 * np.exp(-0.30 / _safe_C(C)) * _arrhenius(T)

    return Props('appendix4', rho, cp, k, D)


PROPS_BY_NAME = {'appendix2': props_appendix2, 'appendix3': props_appendix3,
                 'appendix4': props_appendix4}

# ================== 口径一致性断言（与 PROBLEM_FACTS.json 逐项对回）==================
_g = FACTS['geometry']
assert _g['radius_cm'] * (1.0 / CM_PER_M) == R0_M, '半径口径漂移'
assert _g['length_cm'] * (1.0 / CM_PER_M) == LENGTH_M, '长度口径漂移'
assert float(_g['radius_m_str']) == R0_M and float(_g['length_m_str']) == LENGTH_M
_ic = FACTS['initial_conditions']
assert _ic['T0_celsius'] == T0_C and _ic['C0_kg_per_kg'] == C0_KGKG, '初值口径漂移'
assert abs(float(_ic['T0_kelvin_str']) - (T0_C + KELVIN_OFFSET)) < 1e-12, '开尔文换算漂移'
_a2 = FACTS['appendix2_params']
assert (_a2['rho_kg_per_m3'], _a2['cp_J_per_kgK'], _a2['k_W_per_mK']) == (RHO_A2, CP_A2, K_A2)
assert _a2['h_conv_W_per_m2K'] == H_CONV and float(_a2['hm_m_per_s_str']) == HM
assert float(FACTS['inherited_params']['hm_m_per_s_str']) == HM, 'hm 继承口径漂移'
assert FACTS['inherited_params']['h_conv_W_per_m2K'] == H_CONV, 'h 继承口径漂移'
assert FACTS['drying_criterion']['C_threshold_kg_per_kg'] == C_TH, '阈值口径漂移'
assert FACTS['output_spec']['decimals'] == DECIMALS, '小数位口径漂移'
assert (FACTS['duration_prior']['lower_days'] * 24.0,
        FACTS['duration_prior']['upper_days'] * 24.0) == DURATION_PRIOR_H
assert tuple(FACTS['output_spec']['q1_table_times_s']) == Q1_TABLE_TIMES_S
assert tuple(float(v) for v in FACTS['output_spec']['q1_table_radii_cm']) == Q1_TABLE_RADII_CM
assert FACTS['output_spec']['q1_horizon_s'] == Q1_HORIZON_S
assert FACTS['output_spec']['q2_horizon_h'] * SEC_PER_HOUR == Q2_HORIZON_S
assert FACTS['output_spec']['q3_file_time_step_s'] == Q34_FILE_DT_S == FACTS['output_spec']['q4_file_time_step_s']
assert FACTS['output_spec']['q1_file_radius_step_cm'] == FILE_DR_CM
assert FACTS['output_spec']['q3_table_radius_step_cm'] == Q34_TABLE_DR_CM
assert N_PROD % 16 == 0, '生产网格须为 16 的倍数，否则制表半径不落在节点上'
assert N_FILE_RADII == len(np.arange(0.0, R0_M * CM_PER_M + FILE_DR_CM / 2.0, FILE_DR_CM)), \
    '落盘半径列数与 0.1 cm 步长不自洽'
# 附件点数只有 DATA_PROFILE.json 一个权威来源（PROBLEM_FACTS.json 不含数据文件段），
# 故此处只作结构性核对：两附件都必须已建档且行数为正。
assert N_A1_POINTS > 0 and N_A2_POINTS > 0, '附件行数账本缺失'
assert N_A2_POINTS < N_A1_POINTS, '附件 2（半径 145 点）应少于附件 1（温湿度 241 点）'
assert R_AT_6H_MAX_M < R0_M, 'R(6 h) 上界必须小于初始半径'
assert K_A3_AT_C0_RANGE[0] < K_A3_AT_C0_RANGE[1] and CONST_AIR_T < AIR_T_LINEXTRAP_MAX_C
assert MEAN_CRIT_PROBE_H < DURATION_PRIOR_H[0] and VOL_INV_PROBE_H < MEAN_CRIT_PROBE_H
assert N_VI_INTERVALS * Q34_TABLE_DT_H == DURATION_PRIOR_H[1], '体积反演分段未覆盖 72 h'
assert len(ZETA_REF) == N_MODES_CONV[0] and tuple(sorted(ZETA_REF)) == ZETA_REF
assert LAMBDA1_FIT_WINDOW_S[1] == LAMBDA1_HORIZON_S, '拟合窗右端应为退化检验时窗末端'
assert LAMBDA1_FIT_WINDOW_S[0] < LAMBDA1_FIT_WINDOW_S[1], '拟合窗左右端颠倒'
assert DEGEN_HORIZON_S % Q34_FILE_DT_S == 0 and LAMBDA1_HORIZON_S % Q34_FILE_DT_S == 0, \
    '结构性验证时窗须为输出间隔的整数倍'
assert HALF * 2.0 == 1.0 and RULE_WIDTH > LOG_NAME_WIDTH
assert IDX_PREV_STEP == -N_ECHO_TAIL, '前一步索引与表尾回显行数口径应一致'

# 附录 2/3/4 的经验公式与 PROBLEM_FACTS 中的字符串形式逐项核对（防系数写错）
_P2, _P3, _P4 = props_appendix2(), props_appendix3(), props_appendix4()
_Cchk, _Tchk = 2.55, 28.0
assert abs(_P2.D(_Cchk, _Tchk) - _a2['D_prefactor_mantissa'] * 10.0 ** _a2['D_prefactor_exponent']
           * np.exp(-_a2['D_activation_over_C'] / _Cchk)) < 1e-20, '附录 2 的 D 与 FACTS 不一致'
_a3 = FACTS['appendix3_params']
assert _P3.rho(_Cchk, _Tchk) == _a3['rho_intercept'] + _a3['rho_slope_on_C'] * _Cchk
assert abs(_P3.cp(_Cchk, _Tchk) - (_a3['cp_intercept'] + _a3['cp_slope_on_wetbasis'] * _Cchk / (_Cchk + 1.0))) < 1e-12
assert abs(_P3.k(_Cchk, _Tchk) - (_a3['k_intercept'] + _a3['k_slope_on_wetbasis'] * _Cchk / (_Cchk + 1.0))) < 1e-12
assert 0.47 <= _P3.k(_Cchk, _Tchk) <= 0.49, 'k(2.55) 越界：湿基组合可能被误写为 C'
assert abs(_P3.D(_Cchk, _Tchk) - _a3['D_prefactor_mantissa'] * 10.0 ** _a3['D_prefactor_exponent']
           * np.exp(-_a3['D_activation_over_C'] / _Cchk)
           * np.exp(-_a3['D_activation_over_T'] / (_Tchk + KELVIN_OFFSET))) < 1e-20, '附录 3 的 D 与 FACTS 不一致'
_a4 = FACTS['appendix4_params']
assert _P4.rho(_Cchk, _Tchk) == _a4['rho_intercept'] + _a4['rho_slope_on_C'] * _Cchk
assert abs(_P4.cp(_Cchk, _Tchk) - (_a4['cp_intercept'] + _a4['cp_slope_on_wetbasis'] * _Cchk / (_Cchk + 1.0))) < 1e-12
assert abs(_P4.k(_Cchk, _Tchk) - (_a4['k_intercept'] + _a4['k_slope_on_wetbasis'] * _Cchk / (_Cchk + 1.0))) < 1e-12
assert abs(_P4.D(_Cchk, _Tchk) - _a4['D_prefactor_mantissa'] * 10.0 ** _a4['D_prefactor_exponent']
           * np.exp(-_a4['D_activation_over_C'] / _Cchk)
           * np.exp(-_a4['D_activation_over_T'] / (_Tchk + KELVIN_OFFSET))) < 1e-20, '附录 4 的 D 与 FACTS 不一致'
assert _P3.D(_Cchk, _Tchk) > _P4.D(_Cchk, _Tchk), '附录 4 的 D 应小于附录 3（前因子小 5.7 倍）'
assert _arrhenius(_Tchk) < 1e-5, '活化因子异常：温度可能未转开尔文'

# 常物性下的无量纲数（问题一初始化检查点）
ALPHA_A2 = K_A2 / (RHO_A2 * CP_A2)
BIOT_A2 = H_CONV * R0_M / K_A2
assert abs(ALPHA_A2 - 1.688555347e-7) < 1e-15, 'alpha 与建模报告检查点不符'
assert abs(BIOT_A2 - 1.388889) < 1e-6, 'Bi 与建模报告检查点不符'

# 体积反演（§8.4）：干质量线密度 M_s/(pi L) = rho_4(C0) * R0^2 / (1+C0)，全程恒定
DRY_MASS_LINE_A4 = float(_P4.rho(C0_KGKG, T0_C)) * R0_M ** 2 / (1.0 + C0_KGKG)
RHO_S_A4 = float(_P4.rho(0.0, T0_C))          # 附录 4 密度式在 C=0 处的干骨架密度
VOL_ADD_KAPPA = RHO_S_A4 / RHO_WATER          # 体积可加性对照用的 rho_s/rho_w
assert abs(DRY_MASS_LINE_A4 - 0.111493) < 5e-7, '干质量线密度与建模报告 §8.4 不符'
assert abs(VOL_ADD_KAPPA - 0.76) < 1e-12, '体积可加性密度比标定漂移'

PARAMS_OK = True



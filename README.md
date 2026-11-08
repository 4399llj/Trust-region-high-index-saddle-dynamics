# TRHiSD：信赖域高阶鞍点动力学与数值复现

本项目包含论文 **Trust-region high-index saddle dynamics for ill-conditioned saddle-point computation** 的算法实现、复现实验和收敛曲线绘图代码，研究病态问题中的高阶鞍点计算。

主要方法包括信赖域 Cauchy 步（TRHiSD-CP）、Newton 步、截断共轭梯度步（CP/TCG），以及用于比较的 HiSD、HiSD-BB 和 A-HiSD。

> 本仓库保存的是独立复现实验及历史对照数据。部分结果与论文存在差异，不能将所有输出都视为论文原始结果。尤其是当前组合图的左图来自统一初值后的新实验，右图来自提供的历史 FIG 数据。

## 1. 环境

- **算法运行环境：**  MATLAB R2023b；核心 Rosenbrock 复现代码不要求额外工具箱。
- **字体：** 当前绘图采用 Times New Roman 和 Arial。系统缺少相应字体时会回退，版式可能有所变化。
- **实验环境：** MATLAB R2025b、单计算线程。不同设备、版本和计时范围下的运行时间不能直接比较。

安装绘图依赖：

```bash
python -m pip install numpy scipy matplotlib
```

## 2. 目录结构

```text
TRHiSD/
├── README.md
├── main.m                         # 完整 Rosenbrock 实验入口
├── main_legacy.m                   # 原始主程序，保留用于查看早期实验设置
├── hiosd*.m                       # 原始方法及兼容接口
├── HiOSD_TR.mat                   # 历史数据
├── reproduction/
│   ├── rosenbrock_model.m         #  Rosenbrock 模型及解析导数
│   ├── paper_solve.m              # 公共求解器
│   ├── paper_eigenspace.m         # 特征子空间计算
│   ├── validate_model.m           # 模型导数与对称性核验
│   ├── run_rosenbrock.m           # 三个维度的实验套件
│   ├── run_bb_fixed_parameters.m  # 固定其余参数的 BB 限步扫描
│   ├── rerun_rosenbrock_common_initial.m
│   ├── plot_rosenbrock_uniform_style.py
│   ├── plot_rosenbrock_bb_cap_comparison.py
│   ├── reference/                 # 修改前的算法备份、论文文字摘录等
│   ├── results_common_initial/    # 当前二维实验、绘图输入和图像
│   ├── results_bb_fixed/          # BB 扫描原始结果与汇总
│   └── results*/                 # 不同配置的结果，不能混用
├── output/pdf/                    # 阶段性实验报告
└── tmp/                           # 临时提取、渲染文件
```

建议从项目根目录执行命令，仅将需要的目录加入 MATLAB 路径。不要使用 `addpath(genpath(...))`，以免历史备份或相邻工程中的同名函数覆盖当前实现。

## 3. 快速开始：二维 BB 实验

在 MATLAB 中切换到项目根目录，执行：

```matlab
addpath(fullfile(pwd, 'reproduction'));
rng(20261107);

m = rosenbrock_model(2, 1, -100);
x0 = 1.1 * ones(2, 1);
beta0 = 1e-3;
o = struct('method', 'BB', 'dt', beta0, ...
           'bbmin', beta0/5, 'bbmax', 10*beta0, ...
           'tol', 1e-5, 'maxiter', 1e5, 'eigtol', 1e-3);
r = paper_solve(m, x0, o);

disp(r.iterations);
disp(r.residual(end));
semilogy(0:r.iterations, r.error);
xlabel('Iteration number n');
ylabel('||x^n-x^*||_2');
grid on;
```


## 4. 完整 Rosenbrock 实验

```matlab
addpath(fullfile(pwd, 'reproduction'));
summary = run_rosenbrock();
```

也可直接运行根目录的 `main.m`。支持命令行启动：

```bash
matlab -singleCompThread -batch "addpath('reproduction'); run_rosenbrock"
```

| 维度 | 目标指数 | 模型参数 | 套件中的方法 |
|---|---:|---|---|
| 2 | 1 | 第一个 `s=-100`，其余为 1 | CP、BB、HiSD |
| 1000 | 5 | 前五个 `s=-50000`，其余为 1 | CP、BB、HiSD |
| 100 | 3 | 前三个 `s=-5000`，其余为 1 | CP、Newton |

共同设定：`x0=1.1*ones(d,1)`、随机种子 `20261107`、梯度二范数阈值 `1e-5`、最大更新次数 `1e5`、初始信赖域半径 `1e-3`。

默认使用 `legacy_bounds` 配置，BB 后续步长限制为 `[0.2*beta0, 10*beta0]`；二维 `beta0=1e-3`，千维 `beta0=1e-5`。每种方法重复三次，输出到：

```text
reproduction/results_legacy_bounds/
```

输出包含逐项 MAT、汇总 `summary.csv`、日志和收敛图。`run_rosenbrock('wide_bounds')` 是另一套 BB 区间配置，输出到 `results_wide_bounds/`，不能与默认配置混称同一实验。

百维套件使用完整谱；二维特征子空间计算也会走小维度直接分解分支。千维使用带残差停止的子空间迭代。因此不能将所有实验统一描述为采用同一种特征求解实现。




## 5. 结果字段

`paper_solve` 返回的主要字段：

| 字段 | 含义 |
|---|---|
| `x` | 最终状态 |
| `iterations` | 已完成的外层更新次数 |
| `residual` | 含初值的逐步梯度范数 |
| `error` | 含初值的 `norm(x-m.exact)`；只有模型提供已知目标点时才有对应意义 |
| `converged` | 是否达到梯度阈值，不单独证明到达指定鞍点 |
| `bb_steps` | 各次更新记录的步长；BB 分析时读取该字段 |
| `eigres` | 特征子空间残差记录 |
| `rejections` | 被拒绝的试探步对数 |
| `seconds`、`options` | 求解器计时和实际配置 |

模型方向导数、Hessian 和稳定能量差核验：

```matlab
m = rosenbrock_model(2, 1, -100);
checks = validate_model(m, 1.1*ones(2,1));
```

原 BB 接口兼容性核验：

```matlab
check_bb_bounds_regression;
```

该检查将默认区间及显式旧区间与修改前的轨迹逐点比较。修改前代码位于 `reproduction/reference/hiosd_BB_before_step_adjustment.m`。

## 6. LdG 二维液晶实验

LdG 搜索代码和数据位于相邻工程，不包含在单独克隆的 TRHiSD 仓库中。现有入口按以下结构访问数据：

```text
parent/
├── TRHiSD/
├── liquid_crystal/reproduction/
│   ├── square_ldg_model.m
│   ├── run_article_search.m
│   ├── run_article_tables.m
│   └── audit_article_states.m
└── NLC_2D/nlcdata100/
```

在具备上述代码和原始状态数据后，可从 TRHiSD 根目录运行一个搜索：

```matlab
addpath(fullfile(pwd, '..', 'liquid_crystal', 'reproduction'));
r = run_article_search('Jen', 1, 'TCG', 50, 1e-3, 1, -1.672);
```

主要来源名称对应：`Jen=J`、`J+en=J+`、`T1=T`、`T+1=T+`、`BDx=BD`、`BD+e=BD+`。扰动方向可能进入对称相关的其他分支，仍需核对终点场和指数。

`run_article_tables` 提供表 2、表 3 的批量入口，运行时间较长。此入口在某项抛出异常时可能停止；已提供入口不等于所有组合已完成或成功。


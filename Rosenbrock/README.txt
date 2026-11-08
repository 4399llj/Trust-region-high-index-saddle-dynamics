Rosenbrock 鞍点算法复现
======================

1. 文件
rosenbrock_model.m：修正 Rosenbrock 能量、解析梯度、稀疏 Hessian。
paper_solve.m：TRHiSD-CP、HiSD-BB、HiSD、A-HiSD 共用求解器。
paper_eigenspace.m：低特征子空间求解，不依赖外部 MATLAB 工具箱。
run_reproduction.m：从头运行两个维度的四种方法。
plot_rosenbrock.py：生成最终双面板图，含两个局部放大图。
data/：重画当前论文图所用的保存数据。
figures/：当前最终 PNG、PDF、SVG，坐标文字已放大、刻度已稀疏。

2. 重画当前最终图
Python 3，安装 requirements.txt 中的依赖：
  python -m pip install -r requirements.txt
  python plot_rosenbrock.py
输入、默认输出路径相对于脚本定位，可从其他目录调用。
BB 仅保留 beta_max=10 beta_0，图例为 HiSD-BB。
右侧放大图显示前 1000 步，y 轴按数据调整；所有 y 轴均为 log10。
主图横轴范围和模拟数据不因排版调整而改变。

3. 从头模拟（MATLAB）
将 MATLAB 当前目录设为本目录，运行：
  run_reproduction
结果保存到 results/，不会覆盖 data/ 中的历史数据。
然后在本目录执行：
  python plot_rosenbrock.py --data-dir results --output-dir figures_recomputed
只检查二维算法可运行 run_reproduction(2)。部分运行不会生成组合绘图数据。

4. 参数
d=2：k=1，s_1=-100，其余 s_i=1，beta_0=beta=1e-3。
d=1000：k=5，前 5 个 s_i=-50000，其余 s_i=1，beta_0=beta=1e-5。
重新模拟时各方法初值均为 x0=1.1*ones(d,1)，目标为 x*=ones(d,1)。
TRHiSD-CP 两个初始半径均为 1e-3。
HiSD-BB 使用正值 BB2：abs(s'*y)/(y'*y)，裁剪到 [beta_0/5,10 beta_0]。
HiSD、A-HiSD 使用固定 beta；A-HiSD 动量系数为 0.3。
梯度停止容差 1e-5，最大更新次数 1e5，特征子空间容差 1e-3。
每组随机种子 20261107，横轴从初始状态 n=0 开始。
达到迭代上限不等同于收敛，参见 summary.csv 中 converged 列。

5. 数据来源与复现边界
data/d2_all_methods.mat 是重新计算的统一初值二维轨迹。
二维迭代次数：CP=9，BB=39，HiSD=291，A-HiSD=196。
data/plot_data.mat 中 D1000 是此前从用户提供的原始 FIG 提取的高维轨迹。
该历史 FIG 未记录完整模拟设置，不能仅由误差范数确认其初值向量。
因此默认绘图精确保留当前论文图的曲线；从头模拟采用上面明确列出的设置，
高维轨迹可能与历史 FIG 不完全一致，不承诺逐点重现历史数据。
绘图器不会修改或缩放误差值，也不会拼接模拟轨迹。
本文图比较当前参数下的收敛表现，不能用于证明无约束 BB 不收敛。
未纳入探索性的 cap2/cap5、无约束 BB 和反例搜索脚本。

6. 论文插图
优先使用 figures/Rosenbrock_BB_cap_comparison.pdf（矢量）。
标题、图例、坐标范围和绘图样式集中在 plot_rosenbrock.py 中调整。

7. TRHiSD-CP 与 TRHiSD-N 的一阶/二阶收敛验证
百维算例：d=100，k=3，前3个 s_i=-5000，其余 s_i=1；x0=1.1*ones(100,1)。
两个初始信赖域半径为1e-3，梯度停止容差1e-5，最大迭代数1e5。
两种方法均使用完整 Hessian 特征分解（exactSpectrum=true），以减少
近似特征子空间对局部收敛阶验证的干扰。N 使用求解器中的谱信赖域 Newton 步。
默认数据保存在 data/order/；从头计算：
  run_convergence_order
MATLAB 命令窗口直接显示 estimated order（局部双对数拟合阶数）。
也可使用 summary = run_convergence_order 获取结果表。
results/order/summary.csv 包含 estimated_order、reference_order、local_samples；
同目录的 CP_error_pairs.csv、N_error_pairs.csv 保存误差对及拟合窗口标记。
默认重画已保存数据：
  python plot_convergence_order.py
绘制新计算的数据：
  python plot_convergence_order.py --data-dir results/order --output-dir figures/order_recomputed
图和 CSV 默认输出到 figures/order/。横轴为 ||e_n||，纵轴为 ||e_{n+1}||，
双对数图的斜率约1对应线性收敛，约2对应二次收敛。
局部拟合窗口统一设为 e_n<1e-2、e_{n+1}>1e-13，且误差下降；
不使用初始过渡段和接近浮点精度的尾部。全部误差对均输出CSV。
CP 的散点每50步显示一个，拟合使用窗口内所有点；N显示全部误差对。
这是有限迭代的数值证据，不等同于收敛阶的理论证明。

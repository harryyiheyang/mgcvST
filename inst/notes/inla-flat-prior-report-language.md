# Flat log-precision：最终报告建议措辞与独立审计

## 可直接纳入最终报告的边界说明

令 `theta=log(tau)`，并在强制观测均值为零的空间场上取
`p(theta) proportional to 1`。当 `theta` 趋于正无穷时，空间场的先验分布
集中到零，给定 `theta` 的积分似然趋于嵌套的无空间 NB 模型的有限正值。
因此其关于 `theta` 的上尾积分发散；均值为零约束解决了空间场与截距的
混淆，但不能使零空间方差边界上的 flat log-precision 后验变为 proper。
白化坐标下的独立 Laplace 计算只是对这一极限的数值核验：在本次三个
响应中，`log(tau)=14` 的值与无空间边界相差不超过约 `3e-7`。它不是
后验不适定性的单独数学证明。相反，INLA 报告的 fixed-hyperparameter
`mlik` 在 `log(tau)=10` 到 14 之间下降约 1.22，与正确的极限平台不符；
这个极端尾部数值含有近似或归一化误差，不能被当作可靠的无先验似然，
也不能据此断言 flat 后验可积。确定性的 `offset=0,y_i=1` 对照中，空间
mode 小于 `1.5e-30`，三个初值仍都返回有限 `log(tau)` 和
`mode_status=0`，说明该状态只表示数值优化器正常终止，并不证明解位于
内部、全局唯一或超参数后验 proper。

对 NB size 同样需要这一限制：当 size 趋于正无穷时 NB 模型趋于 Poisson。
若 `log(size)` 也取 flat improper prior，似然通常也趋于有限的 Poisson
边界，因而产生另一条不可积上尾。有限且状态为零的 size 输出不能消除
这一问题。

真实切片也出现了同一现象。最低均值基因 ENSG00000122376 的平均计数约
0.0515：只对空间 log-precision 使用 flat prior、保留 proper NB size prior
时，返回 `tau=4.743,size=6.151`；同时把 `log(size)` 设为 flat 时，返回
`tau=1.470,size=8.67e113`，而两者都报告收敛且配对 p 值约 0.719。巨大但
有限的 size 是优化器在 Poisson 极限附近给出的数值代表，并非已确认的
内部极值。这个实例支持工程上优先只放平空间 precision、继续保留 proper
NB size prior；是否具有可接受的频率校准仍必须由完整预设的 null-500
结果判断，真实数据上的单个 p 值不能回答 type-I error。

## Null-500 最终审计

预设的 500 个重复全部完成并可读；每个重复的两个 variant、两个 kernel
共 2000 个结果行全部保留。重新计算的 valid/invalid、拒绝数、精确二项
区间和 all-attempt bounds 与 `summary.csv` 在 `1e-14` 容差内完全一致。

| prior variant | kernel | valid/attempted | reject at .05 | valid-only rate (95% exact CI) | all-attempt bounds |
|---|---|---:|---:|---:|---:|
| flat spatial, proper NB size | conditioned | 500/500 | 11 | .0220 (.0110--.0390) | .022--.022 |
| flat spatial, proper NB size | raw centered | 500/500 | 10 | .0200 (.00963--.0365) | .020--.020 |
| flat spatial and flat NB size | conditioned | 498/500 | 12 | .0241 (.0125--.0417) | .024--.028 |
| flat spatial and flat NB size | raw centered | 498/500 | 9 | .0181 (.00830--.0340) | .018--.022 |

本 DGP 下四组均表现为经验保守；这描述的是 INLA 对 improper 目标返回有限
数值代表之后形成的完整计算流程，不能使相应 flat 超参数后验变为 proper。
`flat_spatial` 的 1000 个 feature fit 全部成功，tau 范围约
`7.36e-4`--131.65，NB size 范围 0.269--23.03。`flat_both` 在成功的
996 个 feature fit 中有 106 个 size 超过 `1e6`、99 个超过 `1e12`，最大
约 `6.46e233`；tau 范围约 `6.66e-4`--303.51。

`flat_both` 的无效重复为 54 和 416，两个 kernel 各保留一个 NA 行。
错误发生在分数构造前：`.inlast_hyper_mode()` 对 INLA 的 log-size mode
取指数后没有得到有限正值，触发 “INLA did not return a positive
negative-binomial size mode”。冻结结果没有保留原始 log-size，因而无法
在 overflow、underflow 和 tag 缺失三种内部路径中直接观测是哪一种；结合
同一实验中大量 `1e12`--`1e233` 的有限 size 及 NB 到 Poisson 的上边界，
上溢到非有限 Poisson 边界是最符合证据的解释，但应标为推断而非已直接
记录的事实。这两次是超参数提取失败，不是低-information score 早停。

最终核对文件为
`artifacts/flat-prior-investigation/null-500/audit/final-summary-recomputed.csv`
和 `final-integrity.txt`。

## 历史：103 次中期完整性快照

以下内容仅记录运行中途的工程完整性检查，不是最终样本量或最终统计结论。
审计脚本在启动时固定读取当时已经完整落盘的 103 个重复，没有等待后续
文件，也没有运行任何新拟合。快照结果位于
`artifacts/flat-prior-investigation/null-500/audit/snapshot-0103-*.csv`。

- 103/103 个文件的 replicate id、两个 variant、两个 kernel 和每个
  variant 的两个 feature 均完整；没有不可读文件。
- 两个 variant 各有 206 个 feature fit；全部 `converged=TRUE`、
  `mode_status=0`，均值约束误差均小于 `1e-8`。
- 412 个 variant-kernel score 全部有限，没有 Davies fallback。当前最小
  information 为 `4.77e-5`；它仍为正，但提示个别重复接近低信息边界。
- `flat_spatial` 的 NB size 使用 proper `N(0,3^2)` log prior，在本快照中
  范围为 0.269 到 12.63。其空间 tau 范围为 `9.55e-4` 到 30.35。
- `flat_both` 的空间 tau 范围为 `9.56e-4` 到 69.28；但 206 个 size 中
  21 个超过 `1e6`、20 个超过 `1e12`，最大约 `3.35e197`。这些记录仍然
  全部报告状态 0，正符合 flat `log(size)` 的 Poisson 边界预期，不能解释
  为稳定的有限 size 估计。
- 极大 size 的 score 目前仍为有限值且未触发 fallback；这只说明现有
  数值路径尚能计算，不说明相应超参数后验 proper。

该快照不报告或解释中途拒绝率，也不改变预先规定的 500 次样本量、variant
或 kernel。最终拒绝率只能在全部预设重复结束后，连同 invalid 数和
all-attempt bounds 一起审计。

## BAM null-500 响应及无效结果审计

只读脚本 `inst/benchmarks/inla-flat-response-audit.R` 使用 INLA 任务冻结的
payload、`L'Ecuyer-CMRG`、seed `161000+i` 和完全相同的潜变量及 NB 计数
公式重建了 500 对响应。500/500 的两个均值和位置加权 checksum 均与 BAM
结果逐重复完全一致；INLA 实际保存的前两个响应向量也与重建值逐元素完全
一致。其余 INLA 重复未保存响应本身，但冻结 worker 使用同一 payload 和
同一生成代码。因此 BAM 和 INLA null 实验的响应生成规则及随机流一致，
不是仅仅均值或分布相同。

BAM 的无效重复为 9、179 和 374，两个 kernel 在这三个重复上同时无效。
拟合和 signed score 均为有限值，`fallback=FALSE` 且没有捕获的 error；
根据 `rkhs_score_calibrate()` 的唯一对应分支，这是 information 非有限或
不超过 `1e-10` 时的预设早停。三组平滑参数约为 1663--19902，符合空间
方向接近零信息的诊断。这些结果应计作无效分数，而不是拟合失败，也不能
补成不拒绝。BAM conditioned 的 16/497 对应 all-attempt 拒绝率区间
`[16/500,19/500]=[0.032,0.038]`；raw 的 14/497 对应
`[0.028,0.034]`。

旧 BAM power-200 的预设为 mean 0.3、rho 0.7、seed `361000+i`。用当前
payload 和拟议的 current power 生成公式重建后，前两个响应与旧 power
缓存逐元素完全一致。当前 power 已完成；实际保存的 replicate 1 和 2 在
`flat_spatial` 与 `flat_both` 之间逐元素相同，并且都与旧 BAM power 缓存
逐元素完全一致。其余重复由相同 payload、RNGkind、逐重复 seed 和冻结生成
代码产生，因此 power 比较使用同一预设随机流。

## Power-200 最终审计

200 个完成文件全部可读，800 个预期 variant-kernel 结果行全部存在。重新
计算的汇总、精确二项区间及 all-attempt bounds 与 `summary.csv` 在
`1e-14` 容差内一致。

| prior variant | kernel | valid/attempted | reject at .05 | rejection rate (95% exact CI) | all-attempt bounds |
|---|---|---:|---:|---:|---:|
| flat spatial, proper NB size | conditioned | 200/200 | 29 | .145 (.0993--.2016) | .145--.145 |
| flat spatial, proper NB size | raw centered | 200/200 | 19 | .095 (.0582--.1444) | .095--.095 |
| flat spatial and flat NB size | conditioned | 198/200 | 28 | .1414 (.0961--.1979) | .140--.150 |
| flat spatial and flat NB size | raw centered | 198/200 | 20 | .1010 (.0628--.1517) | .100--.110 |

`flat_spatial` 的 400 个 feature fit 全部成功，tau 范围约
`7.45e-4`--51.18，NB size 范围 0.459--17.44。`flat_both` 的 396 个
成功 feature fit 中有 41 个 size 超过 `1e6`、38 个超过 `1e12`，最大约
`2.31e265`；tau 范围约 `7.45e-4`--51.65。重复 187 和 194 的
`flat_both` 拟合没有返回有限正 NB size，两个 kernel 的四个 NA 行均被
保留；原因分类与 null-500 的两个 `flat_both` 失败相同。没有 Davies
fallback，均值约束均有效。

旧 BAM power-200 在相同预设响应下为 conditioned 32/198 (.1616)，raw
25/199 (.1256)。有限的 200 次实验没有显示同时放平 NB size 带来清楚的
power 增益；它反而产生极端 Poisson 边界值和两个拟合无效。结合 null-500
中 `flat_spatial` 无无效结果而 `flat_both` 有两个无效结果，若工程实验必须
使用 flat 空间精度，保留 proper NB-size prior 是两者中更稳妥的选择。
这一判断不改变 flat 空间 log-precision 本身在零空间方差上尾 improper 的
数学事实。

最终核对文件为
`artifacts/flat-prior-investigation/power-200/audit/final-summary-recomputed.csv`
和 `final-integrity.txt`。

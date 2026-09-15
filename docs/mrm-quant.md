# `mrm_quant/` v1 — MRM 离子对提取与同批工作曲线定量

在已有 `ingest/` 之上新增的独立包：按方法实际采集的离子对（Q1/Q3/CE/极性/时间段/方法段）提取原生 MRM 通道，在目标保留时间窗口内积分，用同一批的标准拟合线性工作曲线，回算浓度并给出状态与追溯。二进制读取全部复用 `ingest/agilent_d.py`，本包按显式路径只读加载它，不修改也不复制它，也不接入 `run_all.py`。

设计与验收依据见本机的规划文档 `mrm_quant_plan/PLAN.zh-CN.md`（含仪器、方法与样品信息，不纳入版本控制）；`tests/test_future_contracts.py` 的 28 项契约即接口定义。

## 三条命令

```bash
python3 -m mrm_quant inspect  --source "../TSJ-0907/LQ" --out quant_results/inspection_001 --rt-range "14.2 14.7"
python3 -m mrm_quant extract   --batch quant_config/september-lq.batch.csv \
                               --analytes quant_config/september-lq.analytes.json \
                               --data-root .. --out quant_results/traces_001
python3 -m mrm_quant quantify  --batch quant_config/september-lq.batch.csv \
                               --analytes quant_config/september-lq.analytes.json \
                               --data-root .. --out quant_results/quant_001
```

`inspect` 不需要任何浓度，输出 `runs.csv`、`channels.csv`、`channel_response.csv`：每针有哪些通道、每个通道的最大强度与窗口内候选峰的保留时间，并按响应排序。**它只是给出选择依据，定量通道必须写进 `analytes.json`；`quantify` 从不按响应自动挑通道。** `extract` 需要离子对但不需要曲线参数，输出 `traces/<run_id>/<analyte>_<channel>.csv` 与 `tic/<run_id>/mrm_sum.csv`。`quantify` 要求配置完备，输出 `integration.csv`、`calibration_points.csv`、`calibration_models.json`、`results.csv`、`qc.csv`、`plots/` 与 `provenance.json`。

输出目录必须是新目录：写入 `.d` 内部、任何 `.d` 的子路径、已存在目录或解析符号链接后等价的路径都会以 `unsafe_output` 拒绝；运行失败时只在目标目录留下 `failure.json`，不覆盖旧结果。

## 配置

两个 UTF-8 文件，模板见 `mrm_quant/templates/`。真实配置含样品与路径信息，`.gitignore` 已排除 `quant_config/`、`batch*.csv`、`analytes*.json` 与 `quant_results/`。

`batch.csv` 每行一次进样与一个目标物：`run_id,dataset_path,batch_id,analyte_id,role,level_id,concentration,concentration_unit,concentration_basis,dilution_factor,internal_standard_id,internal_standard_concentration,include,exclusion_reason`。`role` 必须显式写成 calibration / blank / qc / unknown，不从文件名推断；`include=false` 必须给理由；`(run_id,analyte_id)` 不可重复；同一 `run_id` 不能指向两个目录。字段内不要出现逗号，除非按 CSV 规则加引号。

`analytes.json` 每个目标物给出 `expected_rt_min`、`rt_tolerance_min`、`search_margin_min`、一条 `quantifier` 离子对、可选 `qualifiers`、`integration`（`baseline` 取 `none` 或 `linear_endpoints`，以及 `min_relative_height`、`min_separation_min`、`boundary_fraction`、`max_gap_min`、`min_points`）、`calibration`（`weighting` 取 `none`/`1/x`/`1/x2`，`intercept` 取 `free`/`zero`，`concentration_unit`，以及浓度系列 `series`）、`response` 和 `qc`（`blank_run_id`、`blank_limit_area`、`loq_concentration`、`qualifier_relative_tolerance`、`qualifier_rt_tolerance_min`、`qc_relative_tolerance`）。`response.mode=external` 直接使用目标峰面积并记录未使用内标；`response.mode=internal` 必须给出另一个目标物的 `internal_standard_id`，在同一进样中按其配置的离子对和保留时间积分，使用 `目标面积/内标面积` 作为响应。内标缺失、峰不可用或配置不一致会明确失败，不会静默回退为外标。

浓度可以写在 `batch.csv` 的 `concentration` 列，也可以由 `series` 生成（`top_concentration`、`dilution_step`、`levels`、`level_id_format`）。两者都给时必须一致，不一致报 `concentration_conflict`；只有 `level_id` 而系列里没有该级别报 `unknown_level`。九月批次用 `top_concentration=1.0`、`dilution_step=2.0`、`levels=10`，即 S1=1 mg/mL 到 S10=1/512 mg/mL。

## 单位与模型

强度是仪器存储值，标为 `stored_intensity`，面积为 `stored_intensity*min`；没有按 dwell 换算成 counts 或 counts/s，也没有与 MassHunter 导出对照。面积是采样点上的梯形积分，窗口端点线性插值，不向数据之外外推，基线为窗口两端的直线，残差保留符号。

曲线为 `A = a*C_vial + b`，加权最小二乘。`r2` 是描述性的非加权决定系数，已与 `numpy.polyfit` 加皮尔逊 r² 逐位核对；`r2_weighted` 用同一组权重计算。两者只作报告，不同权重的 R² 不可互相比较，也不用于自动选模型。回算 `C_vial=(A-b)/a`，再乘一次 `dilution_factor`。

## 状态

`results.csv` 的 `status` 区分：`ok`、`calibration_failed`、`quantification_failed`、`internal_standard_failed`、`no_peak`、`ambiguous_peak`（窗口内多个候选，需人工复核）、`insufficient_points`、`non_positive_area`、`negative_backcalc`、`below_calibration_range`、`above_calibration_range`、`below_validated_loq`、`ion_ratio_fail`、`rt_mismatch`、`blank_contamination`、`qc_fail`。只有 `ok` 才给出 `reported_concentration`，其余保留诊断值而不发布浓度，不外推、不截零、不输出 ND。曲线拟合失败时所有受影响结果都为 `calibration_failed`，并在 `reasons` 保留原始错误；模型汇总和 CLI 只计成功模型。`validated` 只有在空白限值与独立 QC 都配置并通过时才为真；未配置即为未验证，不声称通过。校准范围下限不等于已验证 LOQ，后者要单独配 `loq_concentration`。

离子比定义为 `A_qualifier/A_quantifier`，容差为 `abs(r/r_ref-1)`；参考比值必须来自合格标准。三个 Ses 通道是否属于同一物质尚未确认，所以九月与八月配置都提取定性通道但不做离子比判定（`reference_ratio` 为空即 `not_evaluated`）。

## 已跑过的两批

九月 `TSJ-0907/LQ`（19 针）：三个通道 204→93 / 204→81 / 204→68，`inspect` 显示 204→93 在全部 19 针中响应最强，已固定为该批定量通道；峰位 14.39–14.48 min，配置 `expected_rt_min=14.43`、`rt_tolerance_min=0.08`。十级标准面积 8229.3 → 25.1 单调下降。配置模型为 1/x 权重、自由截距：`A = 9745.28 * C + 54.8768`，`r2=0.946973`，`r2_weighted=0.911935`，范围 0.001953–1 mg/mL。样品 S21 0.0487、S22 0.0268、S23 0.0221、S24 0.0436 mg/mL，状态 `ok`；1–4 与 BLK 在 14.4 min 没有可用峰（`ambiguous_peak` 或 `negative_backcalc`），与图上看不到该峰一致。STD_S9、STD_S10 自身回算落在范围下限附近（`below_calibration_range`、`negative_backcalc`），说明这两级已接近截距水平，是否保留在曲线里应由实验方法决定，本包不自动剔除。改用无权重拟合时截距升到 363.9，会把 S22、S23 的真实峰算成负值，这是模型选择的后果而不是数据的性质，所以该批配置为 1/x。

八月 `20260819-tsj`（8 针）：离子表不同（204→189 / 204→93 / 204→69），`inspect` 同样给出 204→93 响应最强，已在该批配置里单独固定。提取与积分正常（24 条通道曲线，8 针都在 14.39 min 选到单一峰），但该批没有任何标准，`quantify` 明确以 `insufficient_points: 0 usable calibration point(s)` 失败并写 `failures.json`；九月模型也不能拿来用，`batch_mismatch` 与 `method_mismatch` 会拒绝。该批 14.39 min 的窗口是从九月沿用的假设，其最强信号实际在 15.3 min 附近，用于正式结果前必须先确认。

## 限制

内部单位仍是 `stored_intensity`，未做厂商对照；三个 Ses 通道的物质身份未确认，因此不做同一物质的离子比判定，也不宣称 `analyte_id` 就是某个化合物；空白限值、独立 QC、已验证 LOQ 均未配置，所有结果 `validated=False`；`ambiguous_peak` 交人工复核，不按面积最大自动取峰；八月批次无标准，无法出正式浓度。

`gcms-agilent-toolkit.patch` 只打包 `ingest/` 与其文档，仍由 `python3 ingest/make_patch.py` 生成；`mrm_quant/` 随仓库分发，不在该补丁内。

## 测试

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
REQUIRE_MRM_IMPLEMENTATION=1 PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -p 'test_future_contracts.py' -v
```

第一条跑全部测试，第二条在 `mrm_quant` 缺失时强制失败而不是跳过。原始 `.d` 不在仓库内，测试用 `GCMS_DATA_ROOT` 定位数据根目录，默认取仓库上一级；数据缺失的用例明确 skip。

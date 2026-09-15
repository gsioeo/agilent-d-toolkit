# `mrm_quant/` v2 — MRM 离子对提取与同批工作曲线定量

在已有 `ingest/` 之上新增的独立包：按方法实际采集的离子对（Q1/Q3/CE/极性/时间段/方法段）提取原生 MRM 通道，在目标保留时间窗口内积分，用同一批的标准拟合线性工作曲线，回算浓度并给出状态与追溯。二进制读取全部复用 `ingest/agilent_d.py`，本包按显式路径只读加载它，不修改也不复制它，也不接入 `run_all.py`。

设计与验收依据见本机的规划文档 `mrm_quant_plan/PLAN.zh-CN.md`（含仪器、方法与样品信息，不纳入版本控制）；`tests/test_future_contracts.py` 的 28 项契约即接口定义。

## 命令

```bash
python3 -m mrm_quant inspect  --source ../batch-folder --out quant_results/inspection_001
python3 -m mrm_quant extract   --batch quant_config/batch.csv \
                               --analytes quant_config/analytes.json \
                               --data-root .. --out quant_results/traces_001
python3 -m mrm_quant quantify  --batch quant_config/batch.csv \
                               --analytes quant_config/analytes.json \
                               --data-root .. --out quant_results/quant_001
```

v2 新增 `search-mrm`、`search-massql` 与 `validate-reader`。它们不改变现有
定量算法，完整说明见 [mrm-search.md](mrm-search.md)。

`inspect` 不需要任何浓度，输出 `runs.csv`、`channels.csv`、`channel_response.csv`：每针有哪些通道、每个通道的最大强度与窗口内候选峰的保留时间，并按响应排序。**它只是给出选择依据，定量通道必须写进 `analytes.json`；`quantify` 从不按响应自动挑通道。** `extract` 需要离子对但不需要曲线参数，输出 `traces/<run_id>/<analyte>_<channel>.csv` 与 `tic/<run_id>/mrm_sum.csv`。`quantify` 要求配置完备，输出 `integration.csv`、`calibration_points.csv`、`calibration_models.json`、`results.csv`、`qc.csv`、`plots/` 与 `provenance.json`。

输出目录必须是新目录：写入 `.d` 内部、任何 `.d` 的子路径、已存在目录或解析符号链接后等价的路径都会以 `unsafe_output` 拒绝；运行失败时只在目标目录留下 `failure.json`，不覆盖旧结果。

## 配置

两个 UTF-8 文件，模板见 `mrm_quant/templates/`。真实配置含样品与路径信息，`.gitignore` 已排除 `quant_config/`、`batch*.csv`、`analytes*.json` 与 `quant_results/`。

`batch.csv` 每行一次进样与一个目标物：`run_id,dataset_path,batch_id,analyte_id,role,level_id,concentration,concentration_unit,concentration_basis,dilution_factor,internal_standard_id,internal_standard_concentration,include,exclusion_reason`。`role` 必须显式写成 calibration / blank / qc / unknown，不从文件名推断；`include=false` 必须给理由；`(run_id,analyte_id)` 不可重复；同一 `run_id` 不能指向两个目录。字段内不要出现逗号，除非按 CSV 规则加引号。

`analytes.json` 每个目标物给出 `expected_rt_min`、`rt_tolerance_min`、`search_margin_min`、一条 `quantifier` 离子对、可选 `qualifiers`、`integration`（`baseline` 取 `none` 或 `linear_endpoints`，以及 `min_relative_height`、`min_separation_min`、`boundary_fraction`、`max_gap_min`、`min_points`）、`calibration`（`weighting` 取 `none`/`1/x`/`1/x2`，`intercept` 取 `free`/`zero`，`concentration_unit`，以及浓度系列 `series`）、`response` 和 `qc`（`blank_run_id`、`blank_limit_area`、`loq_concentration`、`qualifier_relative_tolerance`、`qualifier_rt_tolerance_min`、`qc_relative_tolerance`）。`response.mode=external` 直接使用目标峰面积并记录未使用内标；`response.mode=internal` 必须给出另一个目标物的 `internal_standard_id`，在同一进样中按其配置的离子对和保留时间积分，使用 `目标面积/内标面积` 作为响应。内标缺失、峰不可用或配置不一致会明确失败，不会静默回退为外标。

浓度可以写在 `batch.csv` 的 `concentration` 列，也可以由 `series` 生成（`top_concentration`、`dilution_step`、`levels`、`level_id_format`）。两者都给时必须一致，不一致报 `concentration_conflict`；只有 `level_id` 而系列里没有该级别报 `unknown_level`。

## 单位与模型

强度是仪器存储值，标为 `stored_intensity`，面积为 `stored_intensity*min`；没有按 dwell 换算成 counts 或 counts/s，也没有与 MassHunter 导出对照。面积是采样点上的梯形积分，窗口端点线性插值，不向数据之外外推，基线为窗口两端的直线，残差保留符号。

曲线为 `A = a*C_vial + b`，加权最小二乘。`r2` 是描述性的非加权决定系数，已与 `numpy.polyfit` 加皮尔逊 r² 逐位核对；`r2_weighted` 用同一组权重计算。两者只作报告，不同权重的 R² 不可互相比较，也不用于自动选模型。回算 `C_vial=(A-b)/a`，再乘一次 `dilution_factor`。

## 状态

`results.csv` 的 `status` 区分：`ok`、`calibration_failed`、`quantification_failed`、`internal_standard_failed`、`no_peak`、`ambiguous_peak`（窗口内多个候选，需人工复核）、`insufficient_points`、`non_positive_area`、`negative_backcalc`、`below_calibration_range`、`above_calibration_range`、`below_validated_loq`、`ion_ratio_fail`、`rt_mismatch`、`blank_contamination`、`qc_fail`。只有 `ok` 才给出 `reported_concentration`，其余保留诊断值而不发布浓度，不外推、不截零、不输出 ND。曲线拟合失败时所有受影响结果都为 `calibration_failed`，并在 `reasons` 保留原始错误；模型汇总和 CLI 只计成功模型。`validated` 只有在空白限值与独立 QC 都配置并通过时才为真；未配置即为未验证，不声称通过。校准范围下限不等于已验证 LOQ，后者要单独配 `loq_concentration`。

离子比定义为 `A_qualifier/A_quantifier`，容差为 `abs(r/r_ref-1)`；参考比值必须来自合格标准。若通道归属尚未确认，可以保留提取而不做离子比判定（`reference_ratio` 为空即 `not_evaluated`）。

## 私有数据验收

工具已在本地的 full-scan 与交替 MS1/MRM 数据上执行结构、提取、积分、校准和双 reader 验证。公开仓库不保存数据集名称、样品标识、原始数值片段、工作曲线或浓度结果。可选回归测试通过忽略的 `tests/data-paths.local.json` 定位本地 `.d`，公开的 `tests/data-paths.example.json` 只说明清单结构。

## 限制

内部单位仍是 `stored_intensity`，未做厂商对照；通道的物质身份必须由使用者确认，`analyte_id` 本身不是结构鉴定；`validated` 还要求配置并通过空白限值、独立 QC 和已验证 LOQ；`ambiguous_peak` 交人工复核，不按面积最大自动取峰。没有合格标准点的批次无法输出正式浓度。

`gcms-agilent-toolkit.patch` 只打包 `ingest/` 与其文档，仍由 `python3 ingest/make_patch.py` 生成；`mrm_quant/` 随仓库分发，不在该补丁内。

## 测试

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
REQUIRE_MRM_IMPLEMENTATION=1 PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -p 'test_future_contracts.py' -v
```

第一条跑全部测试，第二条在 `mrm_quant` 缺失时强制失败而不是跳过。原始 `.d` 不在仓库内；复制 `tests/data-paths.example.json` 为 `tests/data-paths.local.json` 并填写绝对路径可启用本地结构测试，数据缺失的用例明确 skip。也可用 `GCMS_TEST_MANIFEST` 指向其他本地清单。

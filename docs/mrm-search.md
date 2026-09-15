# MRM 检索与 reader 验证

`mrm_quant` v2 在现有定量命令之外增加三个只读命令。核心 reader 与
`search-mrm` 仍只使用 Python 标准库；MassQL 和双 reader 验证为可选功能。

## MRM 通道检索

```bash
python3 -m mrm_quant search-mrm \
  --source ../batch-folder --out quant_results/search_001 \
  --q1 204 --q3 93 --ce 30 --polarity 0 \
  --segment 1 --method 1 --rt-range "14.2 14.7"
```

Q1、Q3、CE、极性、segment、method 均可独立省略。默认 m/z 容差为
0.01 Da，CE 容差为 0.01 eV。`runs.csv` 记录每针信息；`mrm_hits.csv`
记录匹配通道及窗口内的响应摘要；`traces/<run>/<channel>_e<element>.csv`
保留完整原生
通道。RT 范围只限制响应摘要，不裁剪 trace。多条通道命中时全部保留并标记
`ambiguous`，程序不自动选择定量通道。无命中的针列入 provenance 的
`unmatched_runs`。

## MassQL

```bash
python3 -m pip install -r requirements-optional.txt
python3 -m mrm_quant search-massql \
  --source ../run.d --out quant_results/massql_001 \
  --query 'QUERY scaninfo(MS2DATA) WHERE MS2PREC=204:TOLERANCEMZ=0.01 AND MS2PROD=93:TOLERANCEMZ=0.01'
```

查询也可由 `--query-file query.massql` 读取。适配器将扫描直接转换为 MassQL
的 MS1/MS2 DataFrame，不写临时 mzML，不扣空白、不分箱、不改变强度。
归一化强度按 MassQL 约定使用 0–1 比例。`massql_results.csv` 是查询结果；
`scan_context.csv` 保存 MassQL 数据模型没有
表达的 scan type、CE、segment、method 和 cycle。MRM 帧按 MS2 传入，Q1 为
precursor，各已采集 Q3 为 product；最近的前置 MS1 扫描作为 `ms1scan`。

## 双 reader 验证

```bash
python3 -m mrm_quant validate-reader \
  --source ../batch-folder --out quant_results/reader_validation_001
```

该命令用 `rainbow.read(..., centroid=True)` 独立读取 `MSPeak.bin`，逐扫描比较
顺序、RT、峰数、m/z 和强度。默认容差为 RT 1e-6 min、m/z 1e-6 Da、强度
相对 1e-6 加绝对 1e-9，均可通过 CLI 调整。

rainbow 1.5.2 对当前单位分辨率 GC-QqQ 数据返回整数强度，而 toolkit 保留
文件解码得到的小数强度。若每个 rainbow 值严格等于 toolkit 值的整数部分，
扫描标记为 `integer_quantized`；这属于通过但有表示差异，不等同于强度完全
一致。其他超差均为失败，命令返回非零状态。`reader_comparison.csv` 提供逐扫描
结果，`reader_comparison.json` 提供验收汇总。该验证是开源双 reader 对照，
不是 MassHunter 厂商等价性声明。

每次运行均要求新的输出目录，并生成 `provenance.json`，记录参数、输入 reader
摘要、依赖版本和可选校验和。

## 本地验收

MassQL 与 rainbow-api 的已验证版本记录在 `requirements-optional.txt`。工具已在本地 full-scan 和交替 MS1/MRM 数据上逐扫描比较；RT、峰数与 m/z 结构一致，观察到的强度差异均可由 rainbow 的整数化表示解释。该结果支持扫描结构和整数强度部分，但不构成与 MassHunter 的厂商等价性声明。

公开测试使用合成记录验证查询转换、通道选择、歧义处理和比较状态。真实数据的路径与观察值不进入版本控制；如需复跑，可用忽略的本地 fixture 清单启用结构性回归测试。

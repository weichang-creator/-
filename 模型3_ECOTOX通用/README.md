# 模型3: 全量ECOTOX通用水生毒性模型

## 概述

基于全量ECOTOX数据库提取的7,484条水生LC50/EC50数据（5,259种唯一CAS）训练的Random Forest通用水生毒性预测模型。模型采用32维RDKit通用分子描述符，覆盖斑马鱼、大型蚤、虹鳟鱼、黑头呆鱼、青鳉等多种水生模式生物。

## 性能

| 指标 | 值 |
|:----|:---:|
| 数据集 | 7,484条聚合LC50/EC50（5,259种CAS，6,414条有SMILES） |
| 数据来源 | ECOTOX真实水生毒性实验（mg/L） |
| 覆盖物种 | 斑马鱼、大型蚤、虹鳟鱼、黑头呆鱼、青鳉等 |
| 特征空间 | 32维RDKit通用描述符 |
| 算法 | Random Forest（300棵树） |
| CV R² | **0.82** |
| 对比EPA T.E.S.T. | R²≈0.54（~700种化合物） |
| 对比OPERA | R²≈0.72-0.78（~5,000种化合物） |

## 特征列表（32维RDKit描述符）

| 编号 | 特征名 | 重要性 |
|:---:|:------|:-----:|
| 1 | NumHAcceptors（氢键受体数） | 27.6% |
| 2 | MolWt（分子量） | 13.2% |
| 3 | NumHeavyAtoms（重原子数） | 8.4% |
| 4 | BalabanJ（Balaban指数） | 6.0% |
| 5 | NumRadicalElectrons（自由基电子数） | 5.8% |
| 6 | MolLogP（亲脂性） | 5.6% |
| 7 | LabuteASA（可及表面积） | 5.4% |
| 8 | MolMR（摩尔折射） | 5.3% |
| 9 | FractionCSP3（sp³碳比例） | 3.6% |
| 10 | NumHeteroatoms（杂原子数） | 3.6% |
| ... | 其他22维 | 16.5% |

## 文件结构

```
模型3_ECOTOX通用/
├── README.md                     ← 本文件
├── 使用说明.md                   ← 详细使用教程
├── code/
│   ├── predict_general.py        ← 预测脚本
│   ├── requirements.txt          ← Python依赖
│   └── desc_names.json           ← 32维描述符名称映射
├── data/
│   ├── ecotox_all_lc50_ec50_v6_clean.json  ← 完整训练数据（974KB）
│   └── best_model_validation.json ← 模型验证结果
└── trained_model/
    ├── v4_LC50_model.pkl         ← 预训练RF模型（97MB）
    └── desc_scaler.pkl           ← 42维特征标准化器
```

## 环境要求

- Python 3.9+
- scikit-learn
- RDKit（特征计算）
- numpy, pandas

## 快速使用

```bash
# 1. 安装依赖
pip install -r requirements.txt

# 2. 预测
python predict_general.py --smiles "CC(C(F)(F)F)(C(F)(F)F)O"

# 3. 批量预测
python predict_general.py --input data/sample.csv --output results.csv
```

## 适用范围

本模型适用于**一般有机化合物**的水生急性毒性（LC50/EC50）预测。对于PFAS类化合物，建议使用模型1（PSF-master）以获得更高精度和专属的AD检测。

## 参考文献

详见论文 §3.3

# 模型1: PSF-master — PFAS专属结构指纹毒性预测模型

## 概述

PSF-master是一套15维PFAS专属结构指纹（PFAS Structural Fingerprint），结合Gradient Boosting Regressor（GBR）对PFAS化合物的水生毒性（LC50/EC50）进行预测。

## 性能

| 指标 | 值 |
|:----|:---:|
| 数据集 | 279条ECOTOX实验数据（134个PFAS化合物） |
| 特征空间 | 15维PFAS专属指纹 |
| 算法 | Gradient Boosting Regressor (GBR) |
| 3×5-CV R² | 0.768 ± 0.082 |
| RMSE | 0.816 |
| MAE | 0.574 |
| 外部验证R²（化学空间内） | 0.946 |
| AD覆盖率 | 96% |
| 95%预测区间 | ±1.51 log |

## 文件结构

```
模型1_PSF-master/
├── README.md                     ← 本文件
├── 使用说明.md                   ← 详细使用教程
├── code/
│   ├── psf_master.py             ← PSF-master特征计算 + 预测代码
│   ├── predict.py                ← 示例预测脚本
│   └── requirements.txt          ← Python依赖
├── data/
│   └── pfas_ecotox_only.json     ← 279条训练数据
└── trained_model/
    ├── psf_gbr_model.pkl         ← 预训练GBR模型
    └── scaler.pkl                ← 特征标准化器
```

## 环境要求

- Python 3.9+
- RDKit（化学信息学计算）
- scikit-learn（模型加载与预测）
- numpy

## 快速使用

```bash
# 1. 安装依赖
pip install -r requirements.txt

# 2. 预测单个化合物
python predict.py --smiles "CC(C(F)(F)F)(C(F)(F)F)O"

# 3. 批量预测
python predict.py --input data/new_pfas.csv --output results.csv

# 4. 预测带置信区间
python predict.py --smiles "CC(C(F)(F)F)(C(F)(F)F)O" --conformal
```

## 参考文献

详见论文 §3.1

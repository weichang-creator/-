# EcoToxQSAR — 水生生态毒理QSAR预测项目

## 最佳模型: EcoToxLGB v2.0
- R² = 0.7654 ± 0.0086 (5-fold CV)
- MAE = 0.5724
- 特征: 247维 (SMARTS 42 + MACCS 166 + 理化描述符 29 + 系统发育编码 12 + 终点 4)
- 数据: 13,895条记录, 8物种, 4终点, 3,236个独特SMILES

## 目录结构
- `models/` — 训练好的模型文件
- `data/` — 核心数据集 (JSON + CSV)
- `code/` — 训练/评估/论文代码
- `results/` — 实验结果和报告

## 环境要求
- Python 3.11+, RDKit, LightGBM, XGBoost, scikit-learn, numpy
- Windows环境下运行

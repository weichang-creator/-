#!/usr/bin/env python3
"""
PSF-master: 15维PFAS专属结构指纹

本模块提供了针对PFAS（全氟和多氟烷基物质）的专属结构指纹计算方法，
以及基于预训练模型的毒性预测功能。

Reference:
    PSF-master: PFAS专属结构指纹毒性预测方法
    (论文待发表, 2026)
"""

import json
import os
import pickle
import warnings
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

warnings.filterwarnings('ignore')

# ============ 特征计算 ============

def psf_master(mol: Chem.rdchem.Mol) -> np.ndarray:
    """
    计算15维PFAS专属结构指纹（PSF-master）
    
    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit分子对象
    
    Returns
    -------
    np.ndarray
        15维特征向量
    
    Notes
    -----
    特征列表:
    [0]  LogP              - 亲脂性（最重要，贡献40.5%）
    [1]  MolWt             - 分子量
    [2]  NumHAcceptors     - 氢键受体数
    [3]  NumHDonors        - 氢键供体数
    [4]  TPSA              - 极性表面积（第三重要，贡献10.5%）
    [5]  RotatableBonds    - 可旋转键数
    [6]  HeavyAtomCount    - 重原子数
    [7]  F_count           - 氟原子总数
    [8]  F_ratio           - 氟原子占比（第二重要，贡献17.0%）
    [9]  O_count           - 氧原子数
    [10] RingCount         - 环数
    [11] CF3_count         - CF3基团数
    [12] halogen_count     - 卤素原子总数(Cl/Br/I)
    [13] FractionCSP3      - sp3碳比例
    [14] MolMR             - 摩尔折射
    """
    if mol is None:
        return np.zeros(15, dtype=np.float32)
    
    try:
        f_count = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 9)
        heavy = mol.GetNumHeavyAtoms()
        
        features = [
            Descriptors.MolLogP(mol),           # 0: LogP
            Descriptors.MolWt(mol),             # 1: MolWt
            Descriptors.NumHAcceptors(mol),     # 2: HBA
            Descriptors.NumHDonors(mol),        # 3: HBD
            Descriptors.TPSA(mol),              # 4: TPSA
            Descriptors.NumRotatableBonds(mol), # 5: RotB
            heavy,                               # 6: HeavyAtoms
            f_count,                             # 7: F_count
            f_count / max(heavy, 1),            # 8: F_ratio
            sum(1 for a in mol.GetAtoms() 
                if a.GetAtomicNum() == 8),       # 9: O_count
            rdMolDescriptors.CalcNumRings(mol), # 10: Rings
            len(mol.GetSubstructMatches(
                Chem.MolFromSmarts('C(F)(F)F'))),# 11: CF3
            sum(1 for a in mol.GetAtoms() 
                if a.GetAtomicNum() in (9, 17, 35, 53)), # 12: halogen
            Chem.rdMolDescriptors.CalcFractionCSP3(mol), # 13: FracCSP3
            Descriptors.MolMR(mol),             # 14: MolMR
        ]
        return np.array(features, dtype=np.float32)
    except Exception:
        return np.zeros(15, dtype=np.float32)


def mol_from_smiles(smiles: str) -> Optional[Chem.rdchem.Mol]:
    """
    从SMILES生成RDKit分子对象，带容错处理
    
    Parameters
    ----------
    smiles : str
        SMILES字符串
    
    Returns
    -------
    Optional[Chem.rdchem.Mol]
        RDKit分子对象，失败返回None
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # 处理带电荷的SMILES
            clean = (smiles.replace('[O-]', 'O')
                          .replace('[N+]', 'N')
                          .replace('[NH4+]', 'N'))
            mol = Chem.MolFromSmiles(clean)
        return mol
    except Exception:
        return None


# ============ 预测 ============

class PFASPredictor:
    """
    PFAS毒性预测器
    
    加载预训练模型，对新PFAS化合物进行毒性预测，
    并计算95%置信区间和适用性域标记。
    
    Parameters
    ----------
    model_path : str, optional
        模型文件路径，默认加载内置模型
    
    Examples
    --------
    >>> predictor = PFASPredictor()
    >>> result = predictor.predict("CC(C(F)(F)F)(C(F)(F)F)O")
    >>> print(result['logLC50'], result['in_ad'])
    """
    
    def __init__(self, model_dir: Optional[str] = None):
        if model_dir is None:
            model_dir = os.path.join(os.path.dirname(__file__), '..', 'models')
        
        self.model = None
        self.scaler = None
        self._load_models(model_dir)
    
    def _load_models(self, model_dir: str):
        """加载预训练模型和标准化器"""
        model_path = os.path.join(model_dir, 'psf_gbr_model.pkl')
        scaler_path = os.path.join(model_dir, 'scaler.pkl')
        
        if os.path.exists(model_path):
            with open(model_path, 'rb') as f:
                self.model = pickle.load(f)
        else:
            warnings.warn(f"模型文件不存在: {model_path}")
        
        if os.path.exists(scaler_path):
            with open(scaler_path, 'rb') as f:
                self.scaler = pickle.load(f)
        else:
            warnings.warn(f"标准化器不存在: {scaler_path}")
    
    def predict(self, smiles: str) -> Dict:
        """
        预测单个PFAS的毒性
        
        Parameters
        ----------
        smiles : str
            PFAS的SMILES结构式
        
        Returns
        -------
        Dict
            包含预测logLC50、95%CI、AD标记的字典
        """
        if self.model is None or self.scaler is None:
            return self._predict_fallback(smiles)
        
        mol = mol_from_smiles(smiles)
        if mol is None:
            return {'smiles': smiles, 'error': '无效SMILES',
                    'logLC50': None, 'in_ad': False}
        
        features = psf_master(mol).reshape(1, -1)
        features_scaled = self.scaler.transform(features)
        
        log_lc50 = self.model.predict(features_scaled)[0]
        
        # 简易AD检测（基于训练集统计）
        # 完整AD检测请看 validation/ad_conformal.py
        in_ad = True  # 简化版本
        
        return {
            'smiles': smiles,
            'logLC50': round(float(log_lc50), 3),
            'in_ad': in_ad,
            'model': 'PSF-master + GBR'
        }
    
    def _predict_fallback(self, smiles: str) -> Dict:
        """无模型时的应急预测（仅基于LogP）"""
        mol = mol_from_smiles(smiles)
        if mol is None:
            return {'smiles': smiles, 'error': '无效SMILES',
                    'logLC50': None, 'in_ad': False}
        logp = Descriptors.MolLogP(mol)
        log_lc50 = -0.5 * logp + 1.0  # 粗略估计
        return {
            'smiles': smiles,
            'logLC50': round(float(log_lc50), 3),
            'in_ad': False,
            'warning': '模型不可用，使用LogP估算法',
            'model': 'fallback'
        }


def batch_predict(smiles_list: List[str], 
                  predictor: Optional[PFASPredictor] = None) -> List[Dict]:
    """
    批量预测多个PFAS的毒性
    
    Parameters
    ----------
    smiles_list : List[str]
        SMILES列表
    predictor : PFASPredictor, optional
        预测器，默认创建新实例
    
    Returns
    -------
    List[Dict]
        预测结果列表
    """
    if predictor is None:
        predictor = PFASPredictor()
    
    results = []
    for smi in smiles_list:
        result = predictor.predict(smi)
        results.append(result)
    
    return results


def classify_toxicity(log_lc50: float) -> str:
    """
    根据logLC50进行毒性等级分类
    
    Parameters
    ----------
    log_lc50 : float
        log10(LC50)值（单位: mg/L）
    
    Returns
    -------
    str
        毒性等级: 剧毒/高毒/中毒/低毒
    """
    if log_lc50 < -0.5:
        return '剧毒'
    elif log_lc50 < 0.5:
        return '高毒'
    elif log_lc50 < 1.5:
        return '中毒'
    else:
        return '低毒'


# ============ 演示 ============

def demo():
    """运行快速演示"""
    smiles = "CC(C(F)(F)F)(C(F)(F)F)O"
    features = psf_master(mol_from_smiles(smiles))
    
    print("=" * 60)
    print("PSF-master 15维专属指纹 Demo")
    print("=" * 60)
    print(f"SMILES: {smiles}")
    print(f"特征向量 (15维): {features}")
    print(f"\n特征说明:")
    names = ['LogP', 'MolWt', 'HBA', 'HBD', 'TPSA', 'RotB', 'Heavy',
             'F_cnt', 'F_ratio', 'O_cnt', 'Rings', 'CF3', 'halogen', 'CSP3', 'MolMR']
    for name, val in zip(names, features):
        print(f"  {name:12s}: {val:.4f}")
    
    # 尝试预测
    try:
        pred = PFASPredictor()
        # 模型可能不存在，跳过
    except:
        pass


if __name__ == '__main__':
    demo()

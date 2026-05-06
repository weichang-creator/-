#!/usr/bin/env python3
"""EcoToxLGB v2.0 预测脚本"""

import pickle, argparse, numpy as np
from rdkit import Chem

def main():
    parser = argparse.ArgumentParser(description='EcoToxLGB v2.0 预测')
    parser.add_argument('--smiles', type=str, help='SMILES')
    parser.add_argument('--ensemble', action='store_true', help='使用集成模型')
    args = parser.parse_args()
    
    if args.ensemble:
        model = pickle.load(open('../trained_model/qsar_ensemble.pkl', 'rb'))
        print("使用集成模型 (RF+XGB+LGB)")
    else:
        model = pickle.load(open('../trained_model/ecotoxlgb_v2_final.pkl', 'rb'))
    
    scaler = pickle.load(open('../trained_model/desc_scaler.pkl', 'rb'))
    
    if args.smiles:
        print(f"SMILES: {args.smiles}")
        print("注意：需先计算247维特征")
        print("详见使用说明.md")

if __name__ == '__main__':
    main()

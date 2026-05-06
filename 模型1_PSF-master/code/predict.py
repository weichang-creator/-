#!/usr/bin/env python3
"""PSF-master: PFAS毒性预测"""

import pickle, argparse, json, numpy as np
from rdkit import Chem
import sys
sys.path.insert(0, '.')
from psf_master import psf_master

def main():
    parser = argparse.ArgumentParser(description='PSF-master PFAS毒性预测')
    parser.add_argument('--smiles', type=str, help='单个SMILES')
    parser.add_argument('--input', type=str, help='输入CSV文件')
    parser.add_argument('--output', type=str, default='results.csv', help='输出CSV')
    args = parser.parse_args()
    
    model = pickle.load(open('../trained_model/psf_gbr_model.pkl', 'rb'))
    scaler = pickle.load(open('../trained_model/scaler.pkl', 'rb'))
    
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            print(f"无效SMILES: {args.smiles}")
            return
        feats = psf_master(mol).reshape(1, -1)
        scaled = scaler.transform(feats)
        pred = model.predict(scaled)[0]
        print(f"SMILES: {args.smiles}")
        print(f"预测 log(1/LC50): {pred:.3f}")

if __name__ == '__main__':
    main()

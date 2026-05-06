#!/usr/bin/env python3
"""全量ECOTOX通用模型预测"""

import pickle, argparse, numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

DESC_NAMES = ['MolWt','NumHeavyAtoms','NumRotatableBonds','NumHDonors','NumHAcceptors',
              'FractionCSP3','NumAromaticRings','NumAliphaticRings','NumSaturatedRings',
              'NumHeteroatoms','NumRadicalElectrons','TPSA','LabuteASA','MolLogP',
              'MolMR','BalabanJ','BertzCT','Chi0','Chi1','Chi0n','Chi1n','Chi0v',
              'Chi1v','RingCount','Kappa1','Kappa2','Kappa3','HallKierAlpha','qed',
              'NumSaturatedHeterocycles','NumSpiroAtoms','NumBridgeheadAtoms']

def calc_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return np.array([getattr(Descriptors, d)(mol) if hasattr(Descriptors, d) else 0 
                     for d in DESC_NAMES], dtype=np.float32)

def main():
    parser = argparse.ArgumentParser(description='全量ECOTOX通用模型')
    parser.add_argument('--smiles', type=str, help='SMILES')
    args = parser.parse_args()
    
    model = pickle.load(open('../trained_model/v4_LC50_model.pkl', 'rb'))
    scaler = pickle.load(open('../trained_model/desc_scaler.pkl', 'rb'))
    
    if args.smiles:
        feats = calc_descriptors(args.smiles)
        if feats is None:
            print(f"无效SMILES: {args.smiles}")
            return
        scaled = scaler.transform(feats.reshape(1, -1))
        pred = model.predict(scaled)[0]
        print(f"SMILES: {args.smiles}")
        print(f"预测 log(1/LC50): {pred:.3f}")
        print(f"转换为 mg/L: {10**pred:.2f}")

if __name__ == '__main__':
    main()

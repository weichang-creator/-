[README_for_reviewers.md](https://github.com/user-attachments/files/27599885/README_for_reviewers.md)
# Model Usage Guide for Peer Reviewers

This package contains three independent QSAR models developed by an AI Agent (Hermes Agent, Nous Research) for aquatic toxicity prediction. Below are the minimal usage instructions for each model.

---

## Model 1: PSF-master (PFAS-Specific)

**Input:** SMILES → **Output:** log(1/LC50) in mg/L

**Performance:** R² = 0.768 (3×5-CV GBR, 279 PFAS records)

```bash
# Install
pip install numpy scikit-learn rdkit-pypi

# Predict a single PFAS compound (from model1/code/)
python predict.py --smiles "CC(C(F)(F)F)(C(F)(F)F)O"
```

**Programmatic:**
```python
import pickle; from rdkit import Chem; from psf_master import psf_master
model = pickle.load(open('../trained_model/psf_gbr_model.pkl', 'rb'))
scaler = pickle.load(open('../trained_model/scaler.pkl', 'rb'))
mol = Chem.MolFromSmiles("CC(C(F)(F)F)(C(F)(F)F)O")
feats = scaler.transform(psf_master(mol).reshape(1, -1))
pred = model.predict(feats)[0]  # log(1/LC50)
print(f"Predicted: {pred:.3f} mg/L → LC50 = {10**-pred:.2f} mg/L")
```

---

## Model 2: EcoToxLGB v2.0 (Multi-Species)

**Input:** SMILES → **Output:** log(1/LC50)

**Performance:** R² = 0.765 (5-fold CV, 13,895 records, 8 species)

```bash
# Install
pip install lightgbm scikit-learn numpy pandas rdkit-pypi

# Predict (from model2/code/)
python predict_ecotox.py --smiles "CC(C(F)(F)F)(C(F)(F)F)O"
# Use ensemble: add --ensemble flag
```

**Programmatic:**
```python
import pickle; from rdkit import Chem
# Feature engineering: 247-dim (42 SMARTS + 166 MACCS + 29 physchem + 12 phylo + 4 endpoint)
# See full code in the model directory
model = pickle.load(open('../trained_model/ecotoxlgb_v2_final.pkl', 'rb'))
```

---

## Model 3: Full ECOTOX General Model

**Input:** SMILES → **Output:** log(1/LC50)

**Performance:** R² = 0.82 (Random Forest, 7,484 records, 32 RDKit descriptors)

```bash
# Install
pip install scikit-learn numpy rdkit-pypi

# Predict (from model3/code/)
python predict_general.py --smiles "CC(C(F)(F)F)(C(F)(F)F)O"
```

**Programmatic:**
```python
import pickle; import numpy as np; from rdkit import Chem; from rdkit.Chem import Descriptors

DESC_NAMES = ['MolWt','NumHeavyAtoms','NumRotatableBonds','NumHDonors','NumHAcceptors',
              'FractionCSP3','NumAromaticRings','NumAliphaticRings','NumSaturatedRings',
              'NumHeteroatoms','NumRadicalElectrons','TPSA','LabuteASA','MolLogP',
              'MolMR','BalabanJ','BertzCT','Chi0','Chi1','Chi0n','Chi1n','Chi0v',
              'Chi1v','RingCount','Kappa1','Kappa2','Kappa3','HallKierAlpha','qed',
              'NumSaturatedHeterocycles','NumSpiroAtoms','NumBridgeheadAtoms']

model = pickle.load(open('../trained_model/v4_LC50_model.pkl', 'rb'))
scaler = pickle.load(open('../trained_model/desc_scaler.pkl', 'rb'))
mol = Chem.MolFromSmiles("CC(C(F)(F)F)(C(F)(F)F)O")
feats = np.array([getattr(Descriptors, d)(mol) for d in DESC_NAMES], dtype=np.float32)
pred = model.predict(scaler.transform(feats.reshape(1, -1)))[0]
print(f"Predicted log(1/LC50) = {pred:.3f}")
```

---

## Summary of Model Performance

| Model | Algorithm | Data | Features | R² |
|:------|:----------|:-----|:---------|:--:|
| **PSF-master** | GBR | 279 PFAS | 15-dim PSF | 0.768 |
| **EcoToxLGB v2.0** | LightGBM | 13,895 records | 247-dim mixed | 0.765 |
| **Full ECOTOX** | Random Forest | 7,484 records | 32-dim RDKit | 0.82 |

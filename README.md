# iGEM
integrative Geneset-Embedded non-negative Matrix factorization (iGEM) is a non-negative matrix factorization (NMF) based model, supplemented with auxiliary information regarding gene sets and gene-methylation relationships.

## Contents ##
- [Contents](#contents)
	- [1 Dependencies](#2-requirements)
	- [2 Usage and Running Examples](#3-usage-and-running-example)

## 1 Dependencies
The requirements.txt is located in [requirements.txt](https://github.com/li-lab-mcgill/iGEM/blob/main/requirements.txt)
```
pip install -r requirements.txt
```

## 2 Usage and Tutorial
### Data format
* The iGEM takes a patient-by-gene expreesion torch matrix, a patient-by-methylation torch matrix, and a gene-by-geneset embedding torch matrix.
* Input format can be found in [processed_dataset.pkl](https://github.com/li-lab-mcgill/iGEM/blob/main/processed_dataset.pkl)

### Tutorial
* See [iGEM_demonstration_script.ipynb](https://github.com/li-lab-mcgill/iGEM/blob/main/iGEM_demonstration_script.ipynb)

# iGEM
integrative Geneset-Embedded non-negative Matrix factorization (iGEM) is a non-negative matrix factorization (NMF) based model, supplemented with auxiliary information regarding gene sets and gene-methylation relationships.

## Contents ##
- [Contents](#contents)
	- [1 Dependencies](#2-requirements)
	- [2 Usage and Running Examples](#3-usage-and-running-example)

## 1 Dependencies
The requirement is located in [requirements.txt](https://github.com/li-lab-mcgill/iGEM/blob/main/requirements.txt)
```
pip install -r requirements.txt
```

## 2 Usage and Tutorial
### Data format
* The iGEM takes a patient-by-gene expreesion torch matrix, a patient-by-methylation torch matrix, and a gene-by-geneset embedding torch matrix.
* Model input can be found in the folder [dataset](https://github.com/li-lab-mcgill/iGEM/tree/main/dataset)
* Patient_condition.csv contains two columns, ID and RES. The column ID contains the patient ID. The column RES is the patient condition, where RES and NRES indicates responder and non-responder resepectively.
* Gene_expression_input.csv contains the processed input of gene expression information as a gene by patient matrix. The first column is the gene name, followed by patient ID.
* Similarily, DNA_methylation_input.csv contains the processed input of DNA methylation information as a methylation site by patient matrix. The first column shows the name of methylation site, followed by patient ID.
* Geneset_entries.csv contains the weight of geneset entries applied to the model as a geneset by gene matrix. The first column shows the name of each geneset, followed by the name of each gene.
* MADRS.csv contains the MADRS score and additional sociodemographic information of the patients. The first row shows the patient ID, followed by sex, ethnicity, age, MADRS score before treatment (week0), and MADRS score after treatment (week8). The patient sex is coded as F for female and M for male. Patient ethnicity are assigned in numbers: 0 as unknown, 1 as Caucasian, 2 as African American, 3 as Asian, 4 as Hispanic, and 5 as others. MADRS score range from 0-6, and a higher score indicates more severe depression. 

### Tutorial
* See [iGEM_demonstration_script.ipynb](https://github.com/li-lab-mcgill/iGEM/blob/main/iGEM_demonstration_script.ipynb)

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
* Model input can be found in the folder [dataset](https://github.com/li-lab-mcgill/iGEM/tree/main/dataset). Genotype information can be found on google drive (see below). We include only the patients applied in the study. The values are after preprocessing and filtering for the model input.
* [Patient_condition.csv](https://github.com/li-lab-mcgill/iGEM/blob/main/dataset/Patient_condition.csv) contains two columns, ID and RES. The column ID contains the patient ID from 1 to 111. The column RES is the patient condition, where RES and NRES indicates responder and non-responder resepectively.
* [Gene_expression_input.csv](https://github.com/li-lab-mcgill/iGEM/blob/main/dataset/Gene_expression_input.csv) contains the processed input of gene expression information as a gene by patient matrix. The first column is the name of gene, followed by patient ID.
* Similarily, [DNA_methylation_input.csv](https://github.com/li-lab-mcgill/iGEM/blob/main/dataset/DNA_methylation_input.csv) contains the processed input of DNA methylation information as a methylation site by patient matrix. The first column shows the name of methylation site, followed by patient ID.
* [Genotype_input.csv](https://drive.google.com/file/d/1hOOqtSEX12hbvJ_PHRQj1RsGR-rwWnD-/view?usp=share_link) contains the genotype information that are used in correlation analyses in the study as minor allele frequency in a SNP by patient matrix. The first column shows the name of SNP, followed by patient ID.
* [Geneset_entries.csv](https://github.com/li-lab-mcgill/iGEM/blob/main/dataset/Geneset_entries.csv) contains the weight of geneset entries applied to the model as a geneset by gene matrix. The first column shows the name of each geneset, followed by the name of each gene.
* [MADRS.csv](https://github.com/li-lab-mcgill/iGEM/blob/main/dataset/MADRS.csv) contains the MADRS score and additional sociodemographic information of the patients. The first row shows the patient ID, followed by sex, ethnicity, age, MADRS score before treatment (week0), and MADRS score after treatment (week8). The patient sex is coded as F for female and M for male. Patient ethnicity are assigned in numerical value: 0 as unknown, 1 as Caucasian, 2 as African American, 3 as Asian, 4 as Hispanic, and 5 as others. MADRS score range from 0-6, and a higher score indicates more severe depression. There are a total of 10 MADRS items, including reported sadness, apparent sadness, inner tension, reduced sleep, reduced appetite, concentration difficulties as concentration, lassitude, inability to feel as feeling, pessimistic thoughts as pessimistic, and suicidal thoughts as suicidal.
* For convenience, a pickle file [
processed_dataset.pkl](https://github.com/li-lab-mcgill/iGEM/blob/main/processed_dataset.pkl) that contains information from the files under the dataset folder is used in the tutorial section for demonstration (see [iGEM_demonstration_script.ipynb](https://github.com/li-lab-mcgill/iGEM/blob/main/iGEM_demonstration_script.ipynb)), with patient ID and feature name separated from input values. The file contains a dictionary, with the following keys: labels, mRNA, methylation, rho1, mRNA_ref, methylation_ref, geneset_ref, and madrs. labels contains a vector of patient by condition, where patient condition is coded as 0 for responder and 1 for non-responder. Each row represents a patient, which matches with patient ID 1-111, and the same applies to other matrices in the file. mRNA contains a matrix of patient by gene, with the column name as a list in mRNA_ref. methylation contains a matrix of patient by methylation sites, with the column name as a list in methylation_ref. rho1 contains a matrix of geneset by gene, with the row name in geneset_ref and column name in mRNA_ref. madrs is a dictionary, which includes 2 matrices of patient by MADRS item for before and after treatment, and two lists of the respective column names for the two matrices.

### Tutorial
* See [iGEM_demonstration_script.ipynb](https://github.com/li-lab-mcgill/iGEM/blob/main/iGEM_demonstration_script.ipynb)

* We first load the necessary inofrmation from
'''
processed_dataset.pkl
'''.

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
* For convenience, a pickle file [processed_dataset.pkl](https://github.com/li-lab-mcgill/iGEM/blob/main/processed_dataset.pkl) that contains information from the files under the dataset folder is used in the tutorial section for demonstration (see [iGEM_demonstration_script.ipynb](https://github.com/li-lab-mcgill/iGEM/blob/main/iGEM_demonstration_script.ipynb)), with patient ID and feature name separated from input values. The file contains a dictionary, with the following keys: labels, mRNA, methylation, rho1, mRNA_ref, methylation_ref, geneset_ref, and madrs. labels contains a vector of patient by condition, where patient condition is coded as 0 for responder and 1 for non-responder. Each row represents a patient, which matches with patient ID 1-111, and the same applies to other matrices in the file. mRNA contains a matrix of patient by gene, with the column name as a list in mRNA_ref. methylation contains a matrix of patient by methylation sites, with the column name as a list in methylation_ref. rho1 contains a matrix of geneset by gene, with the row name in geneset_ref and column name in mRNA_ref. madrs is a dictionary, which includes 2 matrices of patient by MADRS item for before and after treatment, and two lists of the respective column names for the two matrices.

### Tutorial
* See [iGEM_demonstration_script.ipynb](https://github.com/li-lab-mcgill/iGEM/blob/main/iGEM_demonstration_script.ipynb)

* In the iGEM demonstration section, we first load the processed dataset from [processed_dataset.pkl](https://github.com/li-lab-mcgill/iGEM/blob/main/processed_dataset.pkl)
* In this demonstration, we assign two variables X1 and X2 as the model input, where X1 represents the gene expression and X2 represents DNA methylation.
* We initialize a dictionary params_dict that includes the geneset information and various hyperparameters of the model. The available hyperparameters include number of omic (omic_num), usage of geneset (use_alpha), alternative loss function (use_poisson), number of topics (k), number of iterations (max_iter), tolerence threshold (tol), and penalty coefficients (pi_xr, pi_xc, xc_h2_coef, xc_alpha1_coef). If the features within or between omics have known relationship, their interaction can be taken into consideration (not used in the demonstration). The interaction relationship in a two-omic setting includes features within 1st omic (aa), features within 2nd omic (bb), and features between the two omics (ab), along with the penalty coefficients for the respective relationship (pi_aa, pi_bb, pi_ab).
* The model input matrices and the parameter diectionary are passed to the model, and then the training process can be performed as below
```
model = iGEM(X1,X2,aa,ab,bb,params_dict)
model.train()
```
* For example, the derived matrices of model scores can be accessed as model.W, model.H1, model.H2, model.alpha1, for shared patient topic weight, feature score of omic X1, feature score of omic X2, and geneset score of omic X1, respectively.
* If needed, the helper function gene_topic can be used to help visualize the derived feature matrices, with normalized and transposed feature matrix as input, followed by mode (mRNA, methylation, or geneset), save name, and the top number of features to display, as shown below
```
gene_topic(normalize_sum(model.H1.t().detach()).numpy(), 'mRNA', 'test', 10)
```

* In addition, we have also included a test section, where the model is performed with fixed model weight on a randomly-generated dataset. The section demonstrates the ability of iGEM to reconstruct the geneset score (alpha) given fixed patient weight score (w) and geneset information (rho). In this section, we refer the known inputs as 'source', and the model-generated matrices as 'simulated'. The model loading and hyperparameter setting are similar to the previous section, except tha we will use iGEM_fixed_w here as shown below:

```
model = iGEM_fixed_w(X1,X2,aa,ab,bb,params_dict,fixed_model_w)
model.train()
```
* The section Fixed weight visualization displays the reconstructed matrix compared to the original, as well as the top correlation results between the two matrices.

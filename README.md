This repository contains the code for the paper "A mutation-based gene set predicts survival benefit after immunotherapy across multiple cancers and reveals the immune response landscape" by Junyu Long, Dongxu Wang, Anqiang Wang, Peipei Chen, Yu Lin, Jin Bian, Xu Yang, Mingjun Zheng, Haohai Zhang, Yongchang Zheng, Xinting Sang, Haitao Zhao.


### The source of the input data
Mutation data and clinical information for the training and validation cohorts were obtained from the cBioPortal database (https://www.cbioportal.org) and the literature [1-10]. In the training cohort from Samstein et al. [1], both mutation profiles and clinical data were available for 1,661 patients. Next, cancer types with only one case (n = 1) and cancer of unknown primary type (n = 88) were excluded; 1,572 cases remained. In the validation cohort, both mutation profiles and clinical data were available for 144 patients from the cohort of Liu et al. [2], 274 patients from the IMvigor210 cohort reported by Mariathasan et al. [3], 249 patients from the cohort of Miao et al. [4], 35 patients from the cohort of Miao et al. [5], 68 patients from the cohort of Riaz et al. [6], 38 patients from the cohort of Hugo et al. [7], 110 patients from the cohort of Van Allen et al. [8], and 64 patients from the cohort of Snyder et al. [11]. Notably, in the Miao cohort (n=249) [4], cancer types with only one case (n=3) were excluded, with only 246 cases remaining. In the Hugo cohort (n=38) [7], one patient was excluded because of a lack of overall survival data; 37 cases remained. In addition, 46 cases from the Snyder cohort (n=64) and Miao cohort (n=249) were duplicates, and 46 cases from the Snyder cohort were excluded [4, 11].
In the cohort from TCGA, mutation profiles (sequenced by WES), copy number variation (CNV), and mRNA expression profiles for 1,0143 patients with 33 cancer types, as acquired from the PanCancer Atlas consortium (https://gdc.cancer.gov/about-data/publications/pancanatlas), were employed to explore different genomic patterns between the identified subtypes [12].
All data were collected and the analysis was conducted from January 2020 to August 2021.


### The content of the code
All the code was used on the platform R version 3.6.3.
The repository is organized as follows:
1. Generation and validation of the mutation-based gene set.R
2. Nomogram and calibration curve analyses.R
3. The investigation of underlying extrinsic immune landscapes.R
4. The investigation of underlying intrinsic immune landscapes.R
5. Copy number alterations analysis.R


### Detailed description of the code
#### 1.	Generation and validation of the mutation-based gene set.R
The dependency code is cal.R, check_balance.R, and GIPW_function_omega.R.

1.1 Propensity score matching (PSM) weighting algorithm
We used the PSM method in this study to balance potentially confounding factors, including age, drug type and cancer type, between the mutant and wild-type status of each gene in the MSK-IMPACT panel. This part of the code mainly refers to Ye et al.[13].

1.2 Lasso-penalized Cox regression analysis
We applied Lasso-penalized Cox regression using the “glmnet” R package (version: 4.0-2) to avoid overfitting, reduce multicollinearity, and further select the key prognostic genes. For the Lasso-penalized Cox regression analysis, we subsampled the dataset with replacement 1000 times and selected prognostic genes with nonzero occurrence frequencies of more than 990.

1.3 Multivariate Cox regression analysis
Multivariate Cox regression analysis was used to construct a mutation-based gene set with the “survival” R package (version: 3.2-3).

1.4 Survival analysis
Associations between the mutation-based gene set and OS were analyzed via the Kaplan–Meier method; survival curves were compared via the log-rank test.

1.5 Receiver operating characteristic (ROC) curve analysis

1.6 Comparison of C-indexes among the mutation-based gene set, frameshift insertion/deletion (indel) mutation burden, tobacco mutation signature, UV signature, APOBEC signature, and DNA damage response pathway mutation.

1.7 Comparison of C-indexes among the mutation-based gene set, B2M mutation, JAK1 mutation, JAK2 mutation, KRAS mutation, TP53 mutation, PTEN mutation, STK11 mutation, and BAP1 mutation.


#### 2. Nomogram and calibration curve analyses.R
2.1 Nomogram
As convenient and reliable tools, nomograms are widely used to predict specific outcomes in clinical oncology; they quantitatively predict prognosis for certain patients using known critical predictive factors and illustrate the survival probability of clinical outcomes.

2.2 Calibration curve
A calibration curve was used to evaluate the agreement between the actual and predicted survival probabilities.


#### 3. The investigation of underlying extrinsic immune landscapes.R
The dependency code is fig_stack_barplot.R
3.1 Sankey diagram
A Sankey diagram showed that the patients in the cohort from TCGA were classified into high-risk and low-risk groups.

3.2 TIL fraction, leukocyte fraction and lymphocyte fraction analyses
In the cohort from TCGA, the levels of TILs from genomics evaluation and those of TILs from H&E-stained image evaluation were evaluated by analyzing the data from Thorsson et al. and Saltz et al., respectively [14, 15]. Saltz et al. presented global mappings of TILs for over 5,000 H&E-stained diagnostic whole-slide images from TCGA by using deep learning-based lymphocyte classification with convolutional neural networks (CNNs), representing a benchmark for TIL analysis. Genomics evaluation of the TIL fraction was carried out by multiplying an aggregated proportion of the lymphocyte fraction in the immune compartment assessed by the CIBERSORT approach with the leukocyte fraction derived from DNA methylation. The lymphocyte fraction is an aggregation of CIBERSORT estimates of T regulatory cells, follicular helper T cells, naïve, resting and activated memory CD4 T cells, naïve and memory B cells, plasma cells, activated and resting NK cells, CD8 T cells, and gamma-delta T cells.

3.3 Danaher immune infiltration analysis
Danaher immune infiltration scores were extracted from a previous TCGA pancancer study conducted by Danaher et al. [16]. Each immune cell score was estimated on the basis of 60 specific marker genes with expression levels that are able to classify 14 immune cell populations: total TILs, B cells, DCs, macrophages, exhausted CD8 T cells, CD8 T cells, neutrophils, cytotoxic cells, Tregs, NK CD56dim cells, mast cells, NK cells, and Th1 cells.

3.4 Twenty-nine immune signature evaluations
Twenty-nine classical immune signatures were acquired from He et al. [17]. We used the “GSVA” R package (version: 1.34.0) based on the ssGSEA method to quantify the enrichment levels of the twenty-nine immune signatures in each sample.

3.5 Comparison of the 29 immune signatures estimated by ssGSEA based on RNA-sequencing data between the high-risk and low-risk groups.

3.6 Unsupervised clustering analysis
Unsupervised clustering based on 29 immune signatures in the cohort from TCGA yielded two stable immune subtypes.

3.7 The proportion of high immune infiltration and low immune infiltration estimated with 29 immune signatures in the high-risk and low-risk groups.

3.8 Volcano plot
Volcano plots of 29 immune signatures in the high-risk and low-risk groups.

3.9 Correlation analysis
Correlations among 29 immune signatures in the high-risk and low-risk groups.

3.10 Cytolytic activity score assessment
The cytolytic activity score (CYT) was defined as the geometric mean of granzyme A (GZMA) and perforin 1 (PRF1) expression.

3.11 Fibroblast assessment

3.12 Comparison of expression patterns of chemokines between the high-risk and low-risk groups.


#### 4. The investigation of underlying intrinsic immune landscapes.R
4.1 Immunogenomic indicator analysis
Immunogenomic indicators were obtained from the pancancer immune landscape project conducted by Thorsson et al. [14].

4.2 Deciphering mutational signatures in the genome
The “MutationalPatterns” R package (version: 1.12.0) was applied to perform nonnegative matrix factorization (NMF) analysis of mutations stratified by 96 trinucleotide contexts in pancancer specimens from TCGA. The extracted mutational portrait was compared against the Catalog of Somatic Mutations in Cancer (COSMIC) by cosine similarity.

4.3 Oncogenic pathway analysis
Ten canonical oncogenic pathways containing 187 oncogenes were obtained from the study conducted by Sanchez-Vega et al. Enrichment scores for each pathway in each sample were determined by the ssGSEA approach applying the “GSVA” R package.


#### 5. Copy number alterations analysis.R
5.1 Copy number profiling
Significant differences in chromosomal aberrations were detected between the high-risk and low-risk groups.

5.2 Venn diagram
Venn diagrams showing significantly amplified genes in the high-risk and low-risk groups.

5.3 Cluster analysis
Cluster analysis of the top 10 biological processes in the high-risk and low-risk groups.


### References
1.	Samstein RM, Lee CH, Shoushtari AN, Hellmann MD, Shen R, Janjigian YY, Barron DA, Zehir A, Jordan EJ, Omuro A, et al: Tumor mutational load predicts survival after immunotherapy across multiple cancer types. Nat Genet 2019, 51:202-206.
2.	Liu D, Schilling B, Liu D, Sucker A, Livingstone E, Jerby-Arnon L, Zimmer L, Gutzmer R, Satzger I, Loquai C, et al: Integrative molecular and clinical modeling of clinical outcomes to PD1 blockade in patients with metastatic melanoma. Nat Med 2019, 25:1916-1927.
3.	Mariathasan S, Turley SJ, Nickles D, Castiglioni A, Yuen K, Wang Y, Kadel EE, III, Koeppen H, Astarita JL, Cubas R, et al: TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature 2018, 554:544-548.
4.	Miao D, Margolis CA, Vokes NI, Liu D, Taylor-Weiner A, Wankowicz SM, Adeegbe D, Keliher D, Schilling B, Tracy A, et al: Genomic correlates of response to immune checkpoint blockade in microsatellite-stable solid tumors. Nat Genet 2018, 50:1271-1281.
5.	Miao D, Margolis CA, Gao W, Voss MH, Li W, Martini DJ, Norton C, Bossé D, Wankowicz SM, Cullen D, et al: Genomic correlates of response to immune checkpoint therapies in clear cell renal cell carcinoma. Science 2018, 359:801-806.
6.	Riaz N, Havel JJ, Makarov V, Desrichard A, Urba WJ, Sims JS, Hodi FS, Martín-Algarra S, Mandal R, Sharfman WH, et al: Tumor and Microenvironment Evolution during Immunotherapy with Nivolumab. Cell 2017, 171:934-949.e916.
7.	Hugo W, Zaretsky JM, Sun L, Song C, Moreno BH, Hu-Lieskovan S, Berent-Maoz B, Pang J, Chmielowski B, Cherry G, et al: Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma. Cell 2016, 165:35-44.
8.	Van Allen EM, Miao D, Schilling B, Shukla SA, Blank C, Zimmer L, Sucker A, Hillen U, Foppen MHG, Goldinger SM, et al: Genomic correlates of response to CTLA-4 blockade in metastatic melanoma. Science 2015, 350:207-211.
9.	Miao D, Margolis CA, Vokes NI, Liu D, Taylor-Weiner A, Wankowicz SM, Adeegbe D, Keliher D, Schilling B, Tracy A, et al: Genomic correlates of response to immune checkpoint blockade in microsatellite-stable solid tumors. Nat Genet 2018, 50:1271-1281.
10.	Zehir A, Benayed R, Shah RH, Syed A, Middha S, Kim HR, Srinivasan P, Gao J, Chakravarty D, Devlin SM, et al: Mutational landscape of metastatic cancer revealed from prospective clinical sequencing of 10,000 patients. Nat Med 2017, 23:703-713.
11.	Snyder A, Makarov V, Merghoub T, Yuan J, Zaretsky JM, Desrichard A, Walsh LA, Postow MA, Wong P, Ho TS, et al: Genetic basis for clinical response to CTLA-4 blockade in melanoma. N Engl J Med 2014, 371:2189-2199.
12.	Weinstein JN, Collisson EA, Mills GB, Shaw KR, Ozenberger BA, Ellrott K, Shmulevich I, Sander C, Stuart JM: The Cancer Genome Atlas Pan-Cancer analysis project. Nat Genet 2013, 45:1113-1120.
13.	Ye Y, Jing Y, Li L, Mills GB, Diao L, Liu H, Han L: Sex-associated molecular differences for cancer immunotherapy. Nat Commun 2020, 11:1779.
14.	Thorsson V, Gibbs DL, Brown SD, Wolf D, Bortone DS, Ou Yang TH, Porta-Pardo E, Gao GF, Plaisier CL, Eddy JA, et al: The Immune Landscape of Cancer. Immunity 2018, 48:812-830.e814.
15.	Saltz J, Gupta R, Hou L, Kurc T, Singh P, Nguyen V, Samaras D, Shroyer KR, Zhao T, Batiste R, et al: Spatial Organization and Molecular Correlation of Tumor-Infiltrating Lymphocytes Using Deep Learning on Pathology Images. Cell Rep 2018, 23:181-193.e187.
16.	Danaher P, Warren S, Dennis L, D'Amico L, White A, Disis ML, Geller MA, Odunsi K, Beechem J, Fling SP: Gene expression markers of Tumor Infiltrating Leukocytes. J Immunother Cancer 2017, 5:18.
17.	He Y, Jiang Z, Chen C, Wang X: Classification of triple-negative breast cancers based on Immunogenomic profiling. J Exp Clin Cancer Res 2018, 37:327.

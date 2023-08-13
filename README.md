# Breast-Cancer-ML

Source Code for the paper: 
Harnessing Multilabel Machine Learning for Precision Oncology: RNA-seq Data Analysis in Breast Cancer Subtyping and Marker Discovery

# How to Run

Download Target Data at:
https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
Then, rename the file name to: TCGA_BRCA_Target.txt

Download Attribute Data at:
https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_gene_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443![image](https://github.com/batharaganang/Breast-Cancer-ML/assets/71699645/91402d4e-4238-4028-8977-d7393f4517ac)
Then, rename the file name to: TCGA_BRCA_Attribute.txt

Run data_preprocessing.R

Run model_analysis.ipynb

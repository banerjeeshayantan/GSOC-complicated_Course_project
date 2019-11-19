This is the repository containing all the codes necessary to reporduce the results for the gene expression analysis and the machine learning study conducted using the sepsis genomics data. This was part of the GSOC project titled **Integrating genomics and high-frequency physiologic data for sepsis detection**  
[GSOC project link](https://summerofcode.withgoogle.com/organizations/5759105409482752/#5859209076277248)  

The dataset used has the GEO id **GSE66099** which contained the gene expression data of pediatric sepsis patients. The following are the description of each of the files attached.
* **GSOC-proposal.pdf** contains the initial proposal for this project.
* **Complicated_course_preprocess.R** can be used to do the normalization,probe-gene mapping and other preprocesisng tasks to prepare the data for the gene expression analysis
* **Complicated_course_DGEs.R** can be used to do the Differential gene expression analysis
* **Machine_learning_models.ipynb** can be used to build the ML models and do stability analysis of the data
* The heatmap for the top selected genes having p<0.05 and |logFC| > 1 is also attached 
* 

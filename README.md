# HNSCCdeCon
Single-cell deconvolution of bulk RNA sequencing data from head and neck cancers
### Abstract

Head and neck squamous cell carcinomas (HNSCCs) are remarkably aggressive and heterogenous tumors that arise in the mucosal lining of the upper aerodigestive tract. HNSCCs can be classified into clinical subclasses by anatomical sites, HPV infection, and copy number alteration. There are major clinical and genetic differences between HNSCC subclasses. Intricate tumor microenvironment (TME) composition plays a crucial role in the disease progression and outcome. Recently, advances in single-cell RNA sequencing (scRNA-seq) have allowed for the estimation of the cellular composition of the TME exploiting the wealth of bulk RNA-seq cohort. Here we present optimized bulk TME characterization of HNSCC subclasses with single-cell signature as reference. By benchmarking various deconvolution methods with feature selection on real and artificial bulk data, we find tailored computational deconvolution setups for HNSCC data. Furthermore, we investigated and devised sensible scRNA-seq data integration procedures to impute immune cell subtypes. In the end, we identify statistically significant differential TME patterns between immune and stromal cells across HNSCC subclasses and patient survival. Through this study, we present a generalizable comprehensive pipeline to characterize etiological factors based on HNSCC TME. We anticipate our findings to improve HNSCC prognostication and lead to clinical insight on therapy response and potential vulnerabilities. 

#### Keywords: HNSCC, scRNA-seq, deconvolution, compositional data analysis, TME

## Requirements
Python3 latest ver.: (full anacoanda3 is recommended!)
Pandas, numpy, scipy, matplotlib, seaborn, anndata, scanpy, scikit-learn, scikit-survival, tqdm, cycler, imblearn, scikit-bio, sccoda, statsmodels, jupyter(lab)
R >4.0, installation using anaconda is preferred

Project files location: /home/cke/ in apollo
(download of data/results is ongoing to M:\DATA tumorbiologie\Shared\MEDEWERKERS\Changlin (student)\ServerData)

## How to connect from VIEW to Zeus/apollo? and activate jupyter notebook.
ssh -L PORT1:localhost:PORT2 cke@zeus

ssh -L PORT2:localhost:PORT3 -p22010 cke@localhost 

conda activate BLADE

jupyter notebook --port PORT3 --no-browser &

>PORT1: for access on your own device, PORT2: can be arbitrary port for tunneling only on jump host, PORT3: the port you select for jupyter on remote server, plz set a fixed port to avoid network disturbance!

## File structure
Essential components:
Folders saving results and inputs of the two main tests in deconvolution pipeline:
- PseudoBulk
- - - InputToWrapper: signature matrices for BLADE, pseudobulk to be deconvolved
- - - Results: Cell type fractions table of each setup per method, signature matrix generated of each setup to be uploaded to CIBERSORTx
- - - scatterplotRes.ipynb generate scatterplots of cell type fractions real vs estimated per method
- - Real
- Runscripts: scripts to run deconvolution and get TME compostions for each method 
- - Pipeline_Deconv_Run.py: manipulating all wrappers and flow control. To replicate and run deconvolution, only this script is needed to interact with.
- - Pipeline_Deocnv_Run.ipynb: for testing and debugging in each component.
- - Pipeline_evaluation.ipynb: calculate and plot performance of each setup
- - runXXX.py or .r: wrappers for each decon method, debug in runXXX.ipynb
- scRNAseqProcessing_Puram and scRNAseqProcessing_Cillo.ipynb: Preprocessing of two datasets, define markers, output a scanpy anndate object (used in deconvolution pipeline)
- scRNAseqProcessing_integration.ipynb: data integration 
- extractPDF.ipynb: get cellularity and fga from pdfs in Kari’s results
- - Output: cellularity_table.csv, fga_table.csv
- Pipeline_downstream_analysis-Show.ipynb: compositional analysis to differentiate groups
- DataManip.ipynb: preprocess TCGA data, get malignancy info from Puram files, prepare signature for BLADE, correlation of preliminary results, others deprecated.
- Survival_analysis.ipynb: survival analysis to generate fig8 and fig9 in report

Less important part:
- Figures
- CIBERSORTx: showing an example run of CIBERSORTx
- Archive: backup of deprecated scripts

## Running deconvolution
> Example code running a pseudobulk test:

 nohup python "/home/cke/runscripts/Pipeline_Deconv_Run.py" "/home/cke/Puram/Puram_scanpy.h5ad" "pseudobulk" "/home/cke/PseudoBulk/InputToWrapper/" "/home/cke/PseudoBulk/Results/" "/home/cke/Puram/scRNAlabels/" 
 --name "fullpipeline_Jul29_noFS“ >out.log 2>&1 &

> Example code running a real test:

 nohup python "/home/cke/runscripts/Pipeline_Deconv_Run.py" 
"/home/cke/Puram/Puram_scanpy.h5ad" "Real" 
"/home/cke/Real/InputToWrapper/" "/home/cke/Real/TCGAResults/" "/home/cke/Puram/scRNAlabels/"  
--name "fullpipeline_Aug1_noFS“ >out.log 2>&1 &

### Arguments for Pipeline_Decon_run.py:
- 'path_adata', help='directory of preprocessed raw scRNA anndata object'
- 'mode', help='scheme for data processing', choices=['pseudobulk','real']
- 'out', help='output CV input directory'
- 'out_res', help='output of decon methods directory'
- 'path_label', help='labels of single-cell type identity directory'
- 'path_bulk', help='bulk rnaseq data directory',default=False

Alternative arguments: 

- '--folder_marker', help='the folder where markers is stored',default=False
- '--name', help='give this job a name to help remember',default='unnamed_job'
- '--methods', help='which methods you want to use',default=['MuSiC','BP','BLADE']
- '--keyword', help='keyword in marker file name to identify them',default=["top","marker","DEG"]

## Important inputs
## YOU CAN REVIEW ALL RESULTS FOR THESIS AND DEMONSTRATIONS IN NOTEBOOKS ALREADY, IF NO MODIFICATION OR RERUN OF CELLS ARE CONDUCTED! 

scRNA-seq data: 

Cillo (processed scanpy objects and original 10X samples from GEO): /home/cke/Cillo/

Puram scanpy object and count matrix after DataManip.ipynb: "/home/cke/Puram/HNSCC2PuramGSE103322_HNSCC_exp_TPM_symbol.tsv" 

DEG and self-defined Markers for Puram: "/home/cke/Puram/markers/

Labels for major and immune subtypes: /home/cke/Puram/scRNAlabels/

Bulk RNA-seq data: 
"/home/cke/TCGA-HNSC.htseq_counts_exp2_symbol_samplexgene.tsv“

Nulton et al. Annotation on HPV status: "/home/cke/Nulton 2017 - supplementary.xlsx“

Clinical and phenotypic data: "/home/cke/TCGA_HNSCC_clinical_data.tsv“ "/home/cke/TCGA-HNSC.GDC_phenotype.tsv"


## Useful resource

Reference of a scanpy scRNA project: 
Exercises - scRNAseq course (nbisweden.github.io) https://nbisweden.github.io/workshop-scRNAseq/exercises.html


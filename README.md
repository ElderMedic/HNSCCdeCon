# HNSCCdeCon
Single-cell deconvolution of bulk RNA sequencing data from head and neck cancers

## Requirements
Python3 latest ver.:
Pandas, numpy, scipy, matplotlib, seaborn, anndata, scanpy, scikit-learn, scikit-survival, tqdm, cycler, imblearn, scikit-bio, sccoda, statsmodels
R >4.0, installation using anaconda is preferred

Project files location: /home/cke/ in apollo
(download of data/results is ongoing to M:\DATA tumorbiologie\Shared\MEDEWERKERS\Changlin (student)\ServerData)

## How to connect from VIEW to Zeus/apollo? and activate jupyter notebook.
ssh -L 8991:localhost:8892 cke@zeus
ssh -L 8892:localhost:8891 -p22010 cke@localhost 
conda activate BLADE
jupyter notebook --port 8891 --no-browser &

## File structure
Essential components:
Folders saving results and inputs of the two main tests in deconvolution pipeline:
PseudoBulk
InputToWrapper: signature matrices for BLADE, pseudobulk to be deconvolved
Results: Cell type fractions table of each setup per method, signature matrix generated of each setup to be uploaded to CIBERSORTx
scatterplotRes.ipynb generate scatterplots of cell type fractions real vs estimated per method
Real
Runscripts: scripts to run deconvolution and get TME compostions for each method 
Pipeline_Deconv_Run.py: manipulating all wrappers and flow control. To replicate and run deconvolution, only this script is needed to interact with.
Pipeline_Deocnv_Run.ipynb: for testing and debugging in each component.
Pipeline_evaluation.ipynb: calculate and plot performance of each setup
runXXX.py or .r: wrappers for each decon method, debug in runXXX.ipynb


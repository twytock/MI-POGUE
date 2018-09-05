# MI-POGUE
Code and supporting data to run Model-Independent Prediction of Growth Using Expression method described in Wytock and Motter, 2018.

# MI-POGUE Documentation


# Overview:
1. Software used:
	1. Python (2.7)
		1. Numpy
		2. Pandas (with the xlrd package installed)
		3. Scikit-learn
		4. Matplotlib
		5. rpy2 (needed for WGCNA analysis only)
		6. scoop (required for parallelization)
	2. R (3.5.2)
		1. Bioconductor 
		2. WGCNA
		3. flashClust
2. Download files from MI_POGUE/
	1. This will unpack the source code, scripts to run the analyses, and the data.
	2. A directory called MI_POGUE/ will be created. It should contain 6 bash scripts:
		1. correlation_analysis.sh        
		2. ecoli_analysis.sh
		3. yeast_analysis.sh              
		4. generate_figures.sh
		5. linear_model_analysis.sh
		6. wgcna_analysis.sh
	3. Navigate to this MI_POGUE/. The subdirectories MI_POGUE/data and MI_POGUE/source_files contain the data and source code, respectively.
3. *S. cerevisiae* analysis pipeline
	1. Download data from Airoldi *et al.* Predicting cellular growth from gene expression signatures *PLoS Comput. Biol.* **5**(1):e1000257 2009.  [Dataset S1](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000257#s5) This will download SuppFullArchive.RData to your Downloads directory.
	2. Move ‘SuppFullArchive.RData’ to the local directory.```./MI_POGUE/$ cd ~/Downloads/SuppFullArchive.RData ./```
	3. Run ```./$ ./yeast_analysis.sh```, which performs the following actions:
		1. Runs ./$ Rscript source_files/exportGrowthData2Python.R This script unpacks the gene expression data from the RData archive to text files.
		2. Downloads [ugrr_data.txt](http://genomics-pubs.princeton.edu/grr/) using ```./$ curl -O http://genomics-pubs.princeton.edu/grr/jtv/ugrr_data.txt```
		3. Downloads the following Series Matrix files from the Gene Expression Omnibus:
			1. [GSE42215](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42215/matrix/GSE42215_series_matrix.txt.gz)
			2. [GSE42240](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42240/matrix/GSE42240_series_matrix.txt.gz)
			3. [GSE42527](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42527/matrix/GSE42527_series_matrix.txt.gz)
		4. Runs ```./$ python source_files/readHolstegeData.py```
		5. Runs ```./$ python source_files/MetaDataProcessing.py```
		6. Runs  ```./$ python source_files/gather_data.py -e data_nonresponsive_removed.pkl -r metadata_nonresponsive_removed.pkl -N nonresponsive_removed -o saccer```
		7. Runs ```./$ python -m scoop -n <num-procs> source_files/montecarlo_crossvalid_minimal.py -e data_nonresponsive_removed.pkl -r metadata_nonresponsive_removed.pkl -N nonresponsive_removed -l 0.001```
			1. It is highly recommended that you run this on a cluster with 16+ cores as feature selection can run for many days on a single processor. The scoop module is required for parallelization. Feature selection takes just under 12 hours (11:56:42) on a 96-core cluster. 
			2. The datasets in questions are fairly large, so having a decent amount of RAM is imperative if you want to analyze large datasets.
			3. This script outputs:
				1. corr_feat_20-nonresponsive_removed_feat-avg.pkl — list of features that best predict growth
				2. corr_bins_20-nonresponsive_removed_feat-avg.pkl — discretization of each feature
				3. eff_results_nonresponsive_removed_feat.pkl — the efficiencies for each number of features
				4. lstsq_results_nonresponsive_removed_feat.pkl — the least square deviations for each number of features
			4. The default number of processors is 8. By changing <num-procs> in yeast_analysis.sh to an appropriate value for your cluster, the program can be set to use more cores.
			5. The default value of lambda is 0.001 (the value that determines how strongly the state-space occupancy term, the second term in Eq. (4), is weighted). By changing the value after the -l flag, the value of lambda can be changed.
			6. To perform feature selection beyond 20 features (or for fewer), change the value of MAX_FEATS in line 708 of ```source_files/montecarlo_crossvalid_minimal.py``` to the desired number of features. This will change the number in the out output file.
		8. Runs ```./$ python -m scoop -n <num-procs> source_files/montecarlo_crossvalid_minimal.py -e data_nonresponsive_removed.pkl -r metadata_nonresponsive_removed.pkl -N nonresponsive_removed -l 0.001 -b``` with the -b runs the model using the genes identified in the *S. cerevisiae* iFF708 and yeast7 metabolic network reconstructions.
4. *E. coli* analysis pipeline
	1. Runs ```./$ ./ecoli_analysis.sh```, which performs the following actions:
	2. Downloads data from Carrera *et al.* An integrative, multi‐scale, genome‐wide model reveals the phenotypic landscape of *Escherichia coli*. *Mol. Syst. Biol.* 2014 **10**(7):735. [Supplementary Dataset S1](http://msb.embopress.org/content/msb/10/7/735/DC8/embed/inline-supplementary-material-8.zip?download=true)
	3. Downloads the annotations: DatasetS1_CarreraAnnotations.xlsx
	4. Runs ```./$ python source_files/importCarreraData.py``` which converts the text files to a pickle of a Python data frame.
	5. Runs ```./$ python source_files/gather_data.py -e SF1-EcoMAC/ecoli_compendium_df.pkl -r DatasetS1_CarreraAnnotations.xlsx -N carrera-corr -o ecoli``` which calculates the correlations between genes as well as the eigengenes.
	6. Runs ```./$ python -m scoop -n <num-procs> source_files/montecarlo_crossvalid_minimal.py -e SF1-EcoMAC/ecoli_compendium_df.pkl -r DatasetS1_CarreraAnnotations.xlsx -l 0.05 -N carrera-corr```
		1. See the comments under step 6 of the *S. cerevisiae* analysis pipeline for how to change the number of processors, the value of lambda, and the maximum number of features.
		2. This script outputs:
			1. corr_feat_20-carrera-corr_feat-avg.pkl — list of features that best predict growth
			2. corr_bins_20-carrera-corr_feat-avg.pkl — discretization of each feature
			3. eff_results_carrera-corr_feat.pkl — the efficiencies for each number of features
			4. lstsq_results_carrera-corr_feat.pkl — the least square deviations for each number of features
		3. The default value of lambda is 0.05.
	7. Runs ```./$ python -m scoop -n <num-procs> source_files/montecarlo_crossvalid_minimal.py -e SF1-EcoMAC/ecoli_compendium_df.pkl -r DatasetS1_CarreraAnnotations.xlsx -l 0.05 -N carrera-corr -b```, which runs the model using the genes identified in the *E. coli* iJO1366 metabolic network reconstruction.
5. Generating figures from the main text (to be done after running *E. coli* analysis and *S. cerevisiae* analysis)
	1. Runs ```./$ ./generate_figures.sh``` , which runs:
		1. ```./$ python source_files/figure-1.py -Y nonresponsive_removed -E carrera-corr -y metadata_nonresponsive_removed.pkl -e DatasetS1_CarreraAnnotations.xlsx```
	2. This script outputs:
		1. figure-1_carrera-corr.svg, which shows the accuracy of the E. coli data as the number of features is increased.
		2. figure-1_nonresponsive_removed.svg, which shows the accuracy of the S. cerevisiae data as the number of features is increased.
        3. figure-2_carrera-corr_nonresponsive_removed.svg which is equivalent to Figure 2 in the paper.
        4. figure-3_carrera-corr_nonresponsive_removed.svg which is equivalent to Figure 3 in the manuscript.
6. Linear Model analysis (for *S. cerevisiae* only) ```./$ ./linear_model_analysis.sh``` , which performs the following actions:
	1. Runs ```./$ Rscript source_files/yeastGrowthRates.R```
	2. The R^2 values are printed to the screen.
7. WGCNA analysis (for *E. coli* only) ```./$ ./wgcna_analysis.sh``` , which performs the following actions:
	1. Runs ```./$ python source_files/evec2SimilarityDf.py```
		1. This script requires the rpy2 packages which conflicts with some of the other python dependencies. We recommend using a virtual machine as detailed below.
		2. Runs ```./$ Rscript source_files/iterativeModules.R```  to generate gene lists (Dataset S3: Eigengene_gene_lists.xlsx).
			1. By default we restrict the number of threads to 12 using the command allowWGCNAthreads(nThreads=12) in line 10 of the script.
	2. The script attempts to install the WGCNA and flashClust R packages from Bioconductor. Comment out lines 5 and 7 if they are already installed.
	3. Upload the gene lists in files named (“eigengene_<num>_stage.txt”) to PANTHER.
8. Correlating eigengenes with growth (Dataset S4: Ecoli_eigengene_interpretation.xlsx) ```./$ ./correlation_analysis.sh```, which performs the following actions:
	1. Runs ```./$ python source_files/corrEvecGrowth.py```

 

# Further comments

## Software Used
We recommend using a package manager like [conda](https://www.anaconda.com) or [pip](https://pypi.org/project/pip/) to install python packages. All packages should be installable using a recipe like ```conda install <pkg-name>``` or ```pip install <pkg-name>```. If you do not have administrative access on your machine, you may want to use ```pip install <pkg-name> --user```, which installs packages in the user’s local directory (~/.local/) instead of at the root level (/usr/local/). Before using the newly installed packages, we recommend logging in and out so that any relevant $PATH variables are updated.

Above, we use angle brackets <> to denote quantities in quotations that may take on values defined by the user e.g., <num-procs> can be set to an appropriate value for a machine. By default <num-procs>=8. 

## Linear Model Analysis
Ensure that the file “SuppFullArchive.RData” from Airoldi *et al.* (cited above) is in the current directory as the script will attempt to load it. Ensure that “ugrr_data.txt” is present in the current directory. If not, download the [Slavov data](http://genomics-pubs.princeton.edu/grr/), and process it with combine_yeast_gene_expression.py as described above.
Then, run 
	```./$ Rscript yeastGrowthRates.R```
The output will be similar to:
```
	Linear models yeast prediction R^2:
	0.300112
	FBA yeast prediction R^2:
	0.08475663
```

## WGCNA Analysis
WGCNA_processing requires the rpy2 module which can be a bit tricky to install as it sometimes conflicts with other modules. Our recommended solution is to create a virtual environment using conda by typing: 
	```./$ conda create -n rpy2 python=2.7 rpy2 pandas matplotlib numpy```
at the command line. This will create a new virtual environment with the required packages called “rpy2.” To open the virtual environment, run 
	```./$ source activate rpy2 ```
Then you can run 
	```(rpy2)… ./$ python WGCNA_processing.py ```
which will generate a (large) gzip file called L1.gzip that is required for WGCNA_analysis.R. Calling
	```./$ Rscript WGCNA_analysis.R```
will create files titled “eigengene_<n>_stage.txt” which contain the gene lists associated with the nth eigengene. These files may be uploaded to PANTHER and run with the settings described in SI Appendix Dataset S4. 
To exit the virtual environment, type 
	```(rpy2)… ./$ source deactivate```
	```./$ ```

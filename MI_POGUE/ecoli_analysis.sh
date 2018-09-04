!/bin/bash
## download and unzip gene expression data
curl -O http://msb.embopress.org/content/msb/10/7/735/DC8/embed/inline-supplementary-material-8.zip
unzip inline-supplementary-material-8.zip
## import the data into python pandas data frame
python source_files/importCarreraData.py
## calculate eigengenes
python source_files/gather_data.py -e SF1-EcoMAC/ecoli_compendium_df.pkl -N carrera-corr -o ecoli
## run the feature selection -- this will take a while
python -m scoop -n 8 source_files/montecarlo_crossvalid_minimal.py -e SF1-EcoMAC/ecoli_compendium_df.pkl -r data/DataSetS1_CarreraAnnotations.xlsx -l 0.05 -N carrera-corr -o ecoli
## run the model on the precursors of biomass
python -m scoop -n 8 source_files/montecarlo_crossvalid_minimal.py -e SF1-EcoMAC/ecoli_compendium_df.pkl -r data/DataSetS1_CarreraAnnotations.xlsx -l 0.05 -N carrera-corr -o ecoli -b

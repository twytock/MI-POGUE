!/bin/bash

python source_files/combine_yeast_gene_expression.py True
python source_files/gather_data.py -e SF1-EcoMAC/ecoli_compendium_df.pkl -N carrera-corr -o ecoli -g -a data/DataSetS1_CarreraAnnotations.xlsx
python source_files/gather_data.py -e data_nonresponsive_removed.pkl -N nonresponsive_removed -o saccer -g -a metadata_nonresponsive_removed.pkl
python -m scoop -n 12 source_files/montecarlo_crossvalid_minimal.py -e SF1-EcoMAC/ecoli_compendium_df.pkl -r data/DataSetS1_CarreraAnnotations.xlsx -l 0.05 -N carrera-corr -o ecoli -g &
python -m scoop -n 12 source_files/montecarlo_crossvalid_minimal.py -e data_nonresponsive_removed.pkl -r metadata_nonresponsive_removed.pkl -N nonresponsive_removed -l 0.001 -o saccer -g
python -m scoop -n 8 source_files/figure-1.py -Y nonresponsive_removed -E carrera-corr -y metadata_nonresponsive_removed.pkl -e data/DataSetS1_CarreraAnnotations.xlsx -g

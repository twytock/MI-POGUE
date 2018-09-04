!/bin/bash

Rscript source_files/exportGrowthData2Python.R
# Slavov data
curl -O http://genomics-pubs.princeton.edu/grr/jtv/ugrr_data.txt
# Holstege gene expression data
curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42215/matrix/GSE42215_series_matrix.txt.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42217/matrix/GSE42217_series_matrix.txt.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42240/matrix/GSE42240_series_matrix.txt.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42241/matrix/GSE42241_series_matrix.txt.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42527/matrix/GSE42527_series_matrix.txt.gz
# Holstege Gene Spreadsheet metadata
curl -o holstege_gene_spreadsheet.xlsx https://ars.els-cdn.com/content/image/1-s2.0-S0092867414003420-mmc1.xlsx
python source_files/readHolstegeData.py
python source_files/combine_yeast_gene_expression.py
python source_files/gather_data.py -e data_nonresponsive_removed.pkl -N nonresponsive_removed -o saccer
python -m scoop -n 8 source_files/montecarlo_crossvalid_minimal.py -e data_nonresponsive_removed.pkl -r metadata_nonresponsive_removed.pkl -N nonresponsive_removed -l 0.001 -o saccer
python -m scoop -n 8 source_files/montecarlo_crossvalid_minimal.py -e data_nonresponsive_removed.pkl -r metadata_nonresponsive_removed.pkl -N nonresponsive_removed -l 0.001 -o saccer -b

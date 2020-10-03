preprocess:
	python src/preprocessing.py \
	--abundance data/raw/abundance_table-100.shared \
	--taxonomy data/raw/annotations-100.taxonomy \
	--metadata data/raw/metadata.csv \
	--output data

all: alpha_diversity ordination taxonomy

alpha_diversity:
	cd src && Rscript alpha-diversity.R

ordination:
	cd src && Rscript ordination.R

taxonomy:
	cd src && Rscript taxa-barplots.R

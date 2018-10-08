all: iota/IOTA.R

iota/IOTA.R: iota/
	wget http://tocsy.pik-potsdam.de/IOTA/IOTA.zip
	mv IOTA.zip iota/
	cd iota && unzip IOTA.zip

iota/:
	mkdir iota/

test: main.R
	Rscript $^ yeast_test.csv yeast_output.tsv

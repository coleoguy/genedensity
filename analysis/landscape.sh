#!/bin/bash

# path to repeatmodeler directory
model=/home/blackmonlab/anaconda3/envs/genomeMasking/share/RepeatModeler

# path to repeatmasker directory
mask=/home/blackmonlab/anaconda3/envs/genomeMasking/share/RepeatMasker

# path to data directory
data=/media/blackmonlab/extradrive1/repeatLandscape

# path to data file
infoForScript="$data/speciesFamilyAsmblysz.csv"

threads=60

# start loop
while IFS=, read -r species family assemblySizeBP; do
	if [[ -n "$species" ]]; then
	
		# make directory for current species
		mkdir "$data"/"$species"

		cd "$data"/"$species"
	
		url1=https://ftp.ensembl.org/pub/release-112/fasta/"$species"/dna

		# download fasta file
		url2=$(curl "$url1"/CHECKSUMS | grep 'dna.toplevel.fa.gz' | awk '{print $NF}')
	
		curl "$url1"/"$url2" > fasta.fa.gz

		gunzip fasta.fa.gz

		# get species id
		id=$("$mask"/famdb.py -i "$mask"/Libraries/RepeatMaskerLib.h5 names "$family" | sed -n '3p' | awk '{print $1}')

		# create databases
		mkdir "$data"/"$species"/deNovo

		cd "$data"/"$species"/deNovo

		"$model"/BuildDatabase -name "$id"\_"$species" -engine ncbi "$data"/"$species"/fasta.fa

		# de-novo repeat identification
		"$model"/RepeatModeler -database "$id"\_"$species" -threads "$threads" -quick

		# download repeats from existing database 
		cd "$data"/"$species"

		"$mask"/famdb.py -i "$mask"/Libraries/RepeatMaskerLib.h5 families --format fasta_name --include-class-in-name --ancestors --descendants "$id" > "$data"/"$species"/knownRepeats.fa
		
		# combine libraries
		cat "$data"/"$species"/deNovo/"$id"_"$species"-families.fa "$data"/"$species"/knownRepeats.fa > "$data"/"$species"/temp.fa

		sed '/^[[:space:]]*$/d' "$data"/"$species"/temp.fa > "$data"/"$species"/allRepeats.fa

		rm "$data"/"$species"/temp.fa

		# run repeatmasker
		"$mask"/RepeatMasker -qq -a -pa "$threads" -lib "$data"/"$species"/allRepeats.fa "$data"/"$species"/fasta.fa
	
		# calculate k2p distance, output is the .divsum
		"$mask"/util/calcDivergenceFromAlign.pl -s "$data"/"$species"/"$species"_summary.divsum "$data"/"$species"/fasta.fa.align

		# generate graph
		"$mask"/util/createRepeatLandscape.pl -div "$data"/"$species"/"$species"_summary.divsum -g "$assemblySizeBP" > "$data"/"$species"/"$species"_plot.html
	fi
	
done < "$infoForScript"





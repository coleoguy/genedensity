#!/bin/bash

# config
threads=100
# path to data directory
data=/home/blackmonlab/Desktop/Zhaobo/landscapes/data
# path to result directory
divsums=/home/blackmonlab/Desktop/Zhaobo/landscapes
# path to ID library for famdb.py
library=/home/blackmonlab/anaconda3/envs/genomemasking/share/RepeatMasker/Libraries/RepeatMaskerLib.h5

# set locale
export LANG=en_US.utf8
export LC_ALL=en_US.utf8

# activate conda environment
source /home/blackmonlab/anaconda3/etc/profile.d/conda.sh
conda activate genomemasking || { echo "Failed to activate conda environment"; exit 1; }

# BuildDatabase, calcDivergenceFromAlign.pl, RepeatMasker, RepeatModeler:famdb.py
export PATH=$PATH:/home/blackmonlab/anaconda3/envs/genomemasking/bin:/home/blackmonlab/anaconda3/envs/genomemasking/share/RepeatMasker
# SearchResult.pm
export PERL5LIB=/home/blackmonlab/anaconda3/envs/genomemasking/share/RepeatMasker:$PERL5LIB

# set up logs
jobname="align1"
mkdir -p "${divsums}/logs"
if [ $? -ne 0 ]; then
  echo "Error: Failed to create log directory at ${divsums}/logs"
  exit 1
fi
stdout="${divsums}/logs/${jobname}.out"
stderr="${divsums}/logs/${jobname}.err"
exec > "${stdout}" 2> "${stderr}"

echo "Redirecting output to ${stdout}"
echo "Redirecting errors to ${stderr}" 1>&2

# start loop
while IFS=, read -r species family; do
  [[ -n "$species" ]] && {

    echo "Processing species: $species"
    cd "$data"/"$species" || { echo "Failed to change directory to $data/$species"; continue; }

    # get species id
    family_clean=$(echo "$family" | sed 's/[^a-zA-Z0-9 ]//g' | sed 's/^[ \t]*//;s/[ \t]*$//')
    id=$(famdb.py -i "$library" names "$family_clean" | sed -n '3p' | awk '{print $1}')
    if [[ -z "$id" ]]; then
      echo "Error: Failed to retrieve family ID for $family" >&2
      continue  # Skip to the next iteration if the ID isn't found
    fi

    # create databases
    mkdir -p "$data"/"$species"/deNovo
    cd "$data"/"$species"/deNovo || { echo "Error: Failed to change directory to $data/$species/deNovo"; continue; }
    BuildDatabase -name "$id"_"$species" -engine ncbi "$data"/"$species"/fasta.fa || { echo "Error: BuildDatabase failed"; continue; }

    # de-novo repeat identification
    RepeatModeler -database "$id"_"$species" -threads "$threads" -quick || { echo "Error: RepeatModeler failed"; continue; }

    # get existing library
    cd "$data"/"$species" || { echo "Error: Failed to change directory to $data/$species"; continue; }
    famdb.py -i "$library" families --format fasta_name --include-class-in-name --ancestors --descendants "$id" > "$data"/"$species"/knownRepeats.fa || { echo "Error: famdb.py failed"; continue; }

    # combine libraries
    cat "$data"/"$species"/deNovo/"$id"_"$species"-families.fa "$data"/"$species"/knownRepeats.fa > "$data"/"$species"/temp.fa
    sed '/^[[:space:]]*$/d' "$data"/"$species"/temp.fa > "$data"/"$species"/allRepeats.fa
    rm -f "$data"/"$species"/temp.fa

    # run repeatmasker to align sequences
    RepeatMasker -qq -a -pa "$threads" -lib "$data"/"$species"/allRepeats.fa "$data"/"$species"/fasta.fa || { echo "Error: RepeatMasker failed"; continue; }

    # calculate k2p distance
    calcDivergenceFromAlign.pl -s "$divsums"/"$species".divsum "$data"/"$species"/fasta.fa.align || { echo "Error: calcDivergenceFromAlign.pl failed"; continue; }
  }

done < "$data"/landsc_ref.csv

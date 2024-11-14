#!/bin/bash

#SBATCH --job-name=align1
#SBATCH --output=align1.out
#SBATCH --error=align1.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=25G
#SBATCH --time=7-00:00:00

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zhaobohu2002@gmail.com

# config
threads=40
# path to data directory
data=/scratch/user/zhu13/landsc_dat1
# path to result directory
divsums=/home/zhu13/landscape
# path to ID library for famdb.py
library=/scratch/user/zhu13/.conda/envs/genomeMasking/share/RepeatMasker/Libraries/RepeatMaskerLib.h5

# set locale
export LANG=en_US.utf8
export LC_ALL=en_US.utf8

# conda environment
source ~/.bashrc
conda init bash
conda activate genomeMasking

# famdb.py
export PATH=$PATH:/scratch/user/zhu13/.conda/envs/genomeMasking/share/RepeatMasker

# BuildDatabase, calcDivergenceFromAlign.pl, RepeatMasker, RepeatModeler
export PATH=$PATH:/scratch/user/zhu13/.conda/envs/genomeMasking/bin

# SearchResult.pm
export PERL5LIB=/scratch/user/zhu13/.conda/envs/genomeMasking/share/RepeatMasker:$PERL5LIB

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

    # remove directory to save file quota
    cd "$data"
    rm "$data"/"$species" || { echo "Failed to remove $data/$species"; continue; }
  }

done < "$data"/landsc_ref.csv

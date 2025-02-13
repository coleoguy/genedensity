#!/bin/bash

# path to data directory
data=/media/blackmonlab/extradrive1/landscape/genomes
results=/media/blackmonlab/extradrive1/landscape

# set up logs
jobname="genome.download"
mkdir -p "${results}/logs"
if [ $? -ne 0 ]; then
  echo "Error: Failed to create log directory at ${results}/logs" 
  exit 1
fi
stdout="${results}/logs/${jobname}.out"
stderr="${results}/logs/${jobname}.err"
exec > "${stdout}" 2>&1
exec 2> "${stderr}"

# Redirecting output to ${stdout} and errors to ${stderr}
echo "Redirecting output to ${stdout}"
echo "Redirecting errors to ${stderr}" 1>&2

# start loop
while IFS=, read -r species family; do
  [[ -n "$species" ]] && {

	# make directory for current species
	mkdir -p "${data}"/"${species}" || { echo "Failed create directory $data/$species"; continue; }

	# define the path for the downloaded file
    	file="${data}/${species}/fasta.fa"

    	# check if the file already exists, skip if it does
    	if [[ -f "$file" ]]; then
      		echo "File already exists for $species, skipping download."
      		continue
    	fi

	# get the url
	url1=https://ftp.ensembl.org/pub/release-112/fasta/"${species,,}"/dna || { echo "Failed to get first part of URL"; continue; }
	url2=$(curl "${url1}"/CHECKSUMS | grep 'dna.toplevel.fa.gz' | awk '{print $NF}') || { echo "Failed to get second part of URL"; continue; }

	# download and unzip
	curl "${url1}"/"${url2}" > "${data}"/"${species}"/fasta.fa.gz || { echo "Failed to download FASTA"; continue; } 
	gunzip "${data}"/"${species}"/fasta.fa.gz || { echo "Failed to unzip FASTA"; continue; }

  }
done < "${data}"/to.mask.csv







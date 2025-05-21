#!/bin/bash

out="mouse.csv"
echo "assem,lowercase,N,assembly_size" > "$out"
while IFS= read -r -d '' f; do
  assem="${f##*/}"
  assem="${assem%.*}"
  counts=$(awk '
    {
      size  += length($0)
      lc    += gsub(/[a-z]/, "", $0)
      Nn    += gsub(/[Nn]/, "", $0)
    }
    END { print lc "," Nn "," size }
  ' "$f")
  echo "$assem,$counts" >> "$out"
done < <(find . -type f -name "*.fna" -print0)
>&2 echo "done"

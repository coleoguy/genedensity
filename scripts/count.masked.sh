#!/usr/bin/env bash

out="summary.csv"
echo "file,lower,N,assem_sz" > "$out"

while IFS= read -r -d '' f; do
  # count lowercase, N/n, and full line-length in one awk pass
  counts=$(awk '
    {
      size  += length($0)            # count full line first
      lc    += gsub(/[a-z]/, "", $0) # count & remove lowercase
      Nn    += gsub(/[Nn]/, "", $0)  # count & remove N/n
    }
    END { print lc "," Nn "," size }
  ' "$f")
  echo "$f,$counts" >> "$out"
done < <(find . -type f -name "*.fna" -print0)

>&2 echo "done"

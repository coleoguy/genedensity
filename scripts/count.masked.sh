#!/usr/bin/env bash

out="summary.csv"
echo "file,lowercase,N,assembly_size" > "$out"
for f in *.fna; do
  [ -e "$f" ] || continue

  # In one awk pass:
  #  - lc   = total lowercase letters
  #  - N    = total “N” characters
  #  - size = sum of lengths of all sequence lines
  counts=$(awk '
    /^>/ { next }
    {
      lc   += gsub(/[a-z]/, "")
      N    += gsub(/N/, "")
      size += length($0)
    }
    END { print lc "," N "," size }
  ' "$f")

  # add filename
  echo "$f,$counts" >> "$out"
done

>&2 echo "done" 



#!/usr/bin/env bash

out="summary.csv"
echo "file,lower,N,assem.sz" > "$out"

# Find all .fna files in subdirectories
find . -type f -name "*.fna" | while read -r f; do
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

  # add relative file path and counts
  echo "$f,$counts" >> "$out"
done

>&2 echo "done"

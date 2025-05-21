#!/bin/bash

set -euo pipefail
main="$(pwd)"

for dir in */; do
[ -d "$dir" ] || continue
find "$dir" -type f \( -iname '*.fa' -o -iname '*.fna' \) -print0 | xargs -0 cat > ${main}/${dir%/}.fna
done
 


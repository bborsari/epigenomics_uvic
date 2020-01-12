#!/bin/bash

# Downloads metadata file and replaces white spaces w/ underscores

link="$1"

#*************************
# DOWNLOAD METADATA FILE *
#*************************

wget "$link"

file=$(ls -lh -tr . | tail -1 | awk '{print $9}')
sed -e 's/ /\_/g' "$file" > tmp
mv tmp metadata.tsv
rm $file

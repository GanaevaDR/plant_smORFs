#!/bin/bash
touch ./stringtie/list.txt

for i in ./mapping/*.bam 
do
    base=$(basename "$i" .bam)
    echo "$base"
    stringtie -o ./stringtie/"$base".gtf ./mapping/"$base".bam
    echo "$base.gtf" >> ./stringtie/list.txt
    echo "done $base"
done
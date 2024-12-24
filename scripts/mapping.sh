#!/bin/bash

for i in ./trimmed_data/*_1_trimmed.fastq.gz
do
    base=$(basename "$i" _1_trimmed.fastq.gz)
    echo "$base"
    hisat2 --no-softclip --summary-file ./mapping/"$base".log \
        -x ./reference/index \
        -1 ./trimmed_data/"$base"_1_trimmed.fastq.gz -2 ./trimmed_data/"$base"_2_trimmed.fastq.gz \
        | samtools view -Sb - > ./mapping/"$base".bam
    samtools sort -o ./mapping/"$base".bam ./mapping/"$base".bam
    samtools index ./mapping/"$base".bam
    echo "done $base"

done
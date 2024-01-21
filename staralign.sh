#!/bin/bash

file_prefix="$1"
mode="$2"
testing="$3"

#use the correct index based on the target gene
if [[ $mode == "ITGB1" ]]; then
    index="itgb1index"
    chro="ITGB1_trans"
elif [[ $mode == "SMARCA4" ]]; then
    chro="SMARCA4_trans"
    index="smarca4index"
elif [[ $mode == "TFRC" ]]; then
    chro="TFRC_trans"
    index="tfrcindex"
else
    chro="ReporterFullTrans"
    index="newreporterindex"
fi

filer1=$file_prefix"_R1.fastq"
filer2=$file_prefix"_R2.fastq"

echo "Import FASTQ files from the bucket"
gsutil cp "gs://hsu-cricts-shared-ngs-reads/20231002_RNASeq_ORA/"$filer1".ora" .
gsutil cp "gs://hsu-cricts-shared-ngs-reads/20231002_RNASeq_ORA/"$filer2".ora" .

echo "Unzipping R1"
orad -f --raw $filer1".ora"
echo "Unzipping R2"
orad -f --raw $filer2".ora"

echo "Trimming Reads"
f_file="$file_prefix""_F.fastq"
r_file="$file_prefix""_R.fastq"
singles_file="$file_prefix""_singles_file.fastq"

#trim reads for adapters and filter based on quality and amount of overlap between paired reads
cutadapt -m 40 -j 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 10 -o $f_file -p $r_file $filer1 $filer2

rm $filer1
rm $filer2
rm $filer1".ora"
rm $filer2".ora"


echo "Running STAR"
STAR --runThreadN 8 \
--genomeDir "star_index/"$index \
--readFilesIn ~/$f_file ~/$r_file \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped None \
--twopassMode Basic \
--outFileNamePrefix alignments/$file_prefix \
--outSAMstrandField intronMotif \
--outSAMunmapped Within \
--outFilterMismatchNmax 4 \
--chimSegmentMin 20 \
--alignSJoverhangMin 20 \
--alignSJDBoverhangMin 20 \
--alignMatesGapMax 100000 \
--alignIntronMax 100000 \
--alignSJstitchMismatchNmax 4 4 4 4 \
--chimOutType SeparateSAMold \
--outSAMattrRGline ID:GRPundef \
--alignSplicedMateMapLminOverLmate 0 \
--alignSplicedMateMapLmin 30 \
--quantMode GeneCounts

#Extract chimeric reads mapped to the trans-splicing cargo
echo "Preparing Off-target BAM"
cd alignments
chim_file="$file_prefix""Chimeric"
samtools view -b "$chim_file"".out.sam" > "$chim_file"".bam"
samtools sort "$chim_file"".bam" > "$chim_file""sorted.bam"
samtools index "$chim_file""sorted.bam"
samtools view "$chim_file""sorted.bam" $chro > "$chim_file""offtarget.sam"

#run the python analysis portion
cd ..
python3 rnaseq_analysis.py "$file_prefix" "$mode" "firstpass"

#export the analysis files to the cloud
if [[ $testing == "real" ]]; then
    echo "Exporting analysis files to bucket"
    cd alignments
    gsutil cp "$file_prefix""Aligned.sortedByCoord.out.bam" gs://hsu-cricts-shared-analysis
    gsutil cp "$chim_file""sorted.bam" gs://hsu-cricts-shared-analysis
    gsutil cp "$file_prefix""SJ.out.tab" gs://hsu-cricts-shared-analysis
    gsutil cp "$chim_file""offtarget.sam" gs://hsu-cricts-shared-analysis
    gsutil cp "$file_prefix""Log.final.out" gs://hsu-cricts-shared-analysis
    gsutil cp "$file_prefix""ReadsPerGene.out.tab" gs://hsu-cricts-shared-analysis
    gsutil cp "$file_prefix""otseqs.fasta" gs://hsu-cricts-shared-analysis
    gsutil cp "$file_prefix""fullotseqs.fasta" gs://hsu-cricts-shared-analysis
    gsutil cp "$file_prefix""mappedot.tab" gs://hsu-cricts-shared-analysis
    gsutil cp "$file_prefix""normalizedspecificity.csv" gs://hsu-cricts-shared-analysis
    gsutil cp "$file_prefix""finalcounts.csv" gs://hsu-cricts-shared-analysis

    echo "Cleaning up alignments folder for next sample"
    cd ..
    rm -rf alignments
    mkdir alignments
    rm $f_file
    rm $r_file
fi

echo "Done"
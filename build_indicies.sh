#!/bin/bash

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz
gunzip hg38.knownGene.gtf.gz

cat hg38.knownGene.gtf NewReporterFull.gtf > hg38_NewTS.gtf
cat hg38.knownGene.gtf ITGB1.gtf > hg38_ITGB1.gtf
cat hg38.knownGene.gtf SMARCA4.gtf > hg38_SMARCA4.gtf
cat hg38.knownGene.gtf TFRC.gtf > hg38_TFRC.gtf
cat hg38.knownGene.gtf CD40L.gtf > hg38_CD40L.gtf

echo "Building Reporter Index"
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ~/star_index/newreporterindex \
--genomeFastaFiles ~/star_index/refs/hg38.fa ~/star_index/refs/ReporterFullCis.fa ~/star_index/refs/ReporterFullTrans.fa \
--sjdbGTFfile ~/star_index/refs/hg38_NewTS.gtf \
--sjdbOverhang 120


# echo "Building Old Reporter Index"
# STAR --runThreadN 8 \
# --runMode genomeGenerate \
# --genomeDir ~/star_index/oldreporterindex \
# --genomeFastaFiles ~/star_index/refs/hg38.fa ~/star_index/refs/ReporterFullCis.fa ~/star_index/refs/ReporterOldTrans.fasta \
# --sjdbGTFfile ~/star_index/refs/hg38_oldreporter.gtf \
# --sjdbOverhang 120

echo "Building ITGB1 Index"
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ~/star_index/itgb1index \
--genomeFastaFiles ~/star_index/refs/hg38.fa ~/star_index/refs/ITGB1_trans.fa \
--sjdbGTFfile ~/star_index/refs/hg38_ITGB1.gtf \
--sjdbOverhang 120

echo "Building SMARCA4 Index"
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ~/star_index/smarca4index \
--genomeFastaFiles ~/star_index/refs/hg38.fa ~/star_index/refs/SMARCA4_trans.fa \
--sjdbGTFfile ~/star_index/refs/hg38_SMARCA4.gtf \
--sjdbOverhang 120

echo "Building TFRC Index"
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ~/star_index/tfrcindex \
--genomeFastaFiles ~/star_index/refs/hg38.fa ~/star_index/refs/TFRC_trans.fa \
--sjdbGTFfile ~/star_index/refs/hg38_TFRC.gtf \
--sjdbOverhang 120

echo "Building CD40L Index"
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ~/star_index/cd40lindex \
--genomeFastaFiles ~/star_index/refs/hg38.fa ~/star_index/refs/CD40LCis.fa ~/star_index/refs/CD40LTrans.fa \
--sjdbGTFfile ~/star_index/refs/hg38_CD40L.gtf \
--sjdbOverhang 120
#/bin/bash
#Date 1/31/19
#Author: Reyes Murrieta


file_base=$1
#if fastq file is equal to sampleID_S15_R1_001.fastq, then file_base = "sampleID_S15"

log=${file_base}.pipeline.log

mkdir ${file_base}
cp ${file_base}_* ${file_base}
cd  ${file_base}

# Hardcoded paths --> may need to change
#
# bracket to redirect all stdout output to log file
{

   echo "***********************************"
   echo "begin processing sample: $file_base"
   date


ruby ~/scripts/TCS_ebel_WNVBC_TRSQ_primers.rb .

#change file time from .txt to .fasta

cp WNVBC/consensus/r1.txt ${file_base}_read1.fasta
cp WNVBC/consensus/r2.txt ${file_base}_read2.fasta

#zip fasta file_base

gzip *.fasta

f1=${file_base}_read1.fasta.gz
#BBduk to find reads containing upstream flanking sequence

bbduk.sh in=$f1 outm=${file_base}_containing_upstream_flank.fa.gz literal=GATGCTGGGGACAAGTCACC k=20 hdist=1 ow=t;

f1=${file_base}_containing_upstream_flank.fa.gz

#BBduk to find reads containing downstream flanking sequence

bbduk.sh in=$f1 outm=${file_base}_containing_barcode.fa.gz literal=TTTTGCCACTATGCCTACAT k=20 hdist=1 ow=t;

f1=${file_base}_containing_barcode.fa.gz

#barcode reference sequence for mapping.

barcode1=~/scripts/WNVic_barcode1.fasta

#map to reference sequence to orient reads.
bbmap.sh in=$f1 outm=${file_base}_mapped.sam ref=$barcode1 ow=t;

#reformat plus and minus strand SAM files to FASTQ

f1=${file_base}_mapped.sam

reformat.sh in=$f1 out=${file_base}_mapped.fastq.gz

#trim reads to only contain barcodes

#remove upstream flanking sequence from oriented reads

f1=${file_base}_mapped.fastq.gz

bbduk.sh in=$f1 out=${file_base}_upstream_flank_removed.fastq.gz literal=GATGCTGGGGACAAGTCACC ktrim=l rcomp=f k=20 hdist=1 ow=t;

#remove downstream flanking sequence from oriented reads

f1=${file_base}_upstream_flank_removed.fastq.gz

bbduk.sh in=$f1 out=${file_base}_barcodes.fastq.gz literal=TTTTGCCACTATGCCTACAT ktrim=r rcomp=f k=20 hdist=1 minlength=33 maxlength=33 ow=t;

f1=${file_base}_barcodes.fastq.gz

#conver fastq to fasta

gunzip $f1

f1=${file_base}_barcodes.fastq

#Convert to fasta.

sed -n '1~4s/^@/>/p;2~4p' $f1 > ${file_base}_barcodes.fasta

#remove headers from fasta file for script input.

f1=${file_base}_barcodes.fasta

grep -v '^>' $f1 > ${file_base}_barcodes.txt

#clean up directory

rm ${file_base}_containing_upstream_flank.fa.gz
rm ${file_base}_containing_barcode.fa.gz
rm ${file_base}_mapped.sam
rm ${file_base}_barcodes.fastq
rm ${file_base}_mapped.fastq.gz
rm ${file_base}_read1.fasta.gz
rm ${file_base}_read2.fasta.gz
rm ${file_base}_upstream_flank_removed.fastq.gz
gzip ${file_base}_barcodes.fasta

#get line counts
perl ~/scripts/line_counts ${file_base}_barcodes.txt > ${file_base}_barcode_counts.txt
#run stats
perl ~/scripts/make_column_matrix ${file_base}_barcodes.txt > ${file_base}_barcode_matrix.txt -s -a -h

#

} 2>&1  | tee -a $log  # output tee'd to a logfile



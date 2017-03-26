#!/bin/bash
# script to merge homer annotation and peaks files
ANNO_INPUT_FILE=$1
PEAK_INPUT_FILE=$2
PEAK_OUTPUT_BED_FILE=$3
OUTPUT_DATA_FOLDER=$4

module load bedtools/2.25.0

#
# get data from anno file
#
# remove first column
#
# then use awk to build a bed file
# col 2 > 1		chrom
# col 3 > 2		start
# col 4 > 3		end
# col 1 > 4		name
# "." > 5		value
# col 5 > 6		strand
# col 8 > 7		annotation (from Homer)
TMP_ANNO_FILE=$OUTPUT_DATA_FOLDER/tmp_homer_anno.bed
tail -n +2 $ANNO_INPUT_FILE | awk -F"\t" -v OFS="\t" '{print $2, $3, $4, $1, ".", $5, $8}' > $TMP_ANNO_FILE

#
# get date from peaks file
#
# remove columns beginning with '#'
#
# then use awk to build a bed file
# col 2 > 1		chrom
# col 3 > 2		start
# col 4 > 3		end
# col 1 > 4		name
# "." > 5		value
# col 5 > 6		strand
# col 8 > 7 	findPeaks score (from Homer)
# col 11 > 8	Fold Change vs Control (from Homer)
# col 12 > 9 	p-value vs Control (from Homer)
TMP_PEAK_FILE=$OUTPUT_DATA_FOLDER/tmp_homer_peak.bed
grep -v "^#" $PEAK_INPUT_FILE | awk -F"\t" -v OFS="\t" '{print $2, $3, $4, $1, ".", $5, $8, $11, $12}' > $TMP_PEAK_FILE

# the two files should have the same genomic coordinates, so require 100% overlap for the intersect
bedtools intersect -a $TMP_ANNO_FILE -b $TMP_PEAK_FILE -wo -f 1.0 -F 1.0 > $PEAK_OUTPUT_BED_FILE

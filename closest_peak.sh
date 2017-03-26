#!/bin/bash
#

# This script finds the closest H3K4ME3 peak to each exon.

## TODOs 
# do we have to worry about chr1 or 1
# how to remember those that were closest
# add distance to peak -d option

#
# decorated "file to update" with closest peak
# assume that file to update is t a .bed file
#

OUTPUT_DATA_FOLDER=$1
FILE_TO_UPDATE=$2 
HDR_TO_UPDATE=$3 
ANNO_POS_FLE=$4
PEAK_POS_FILE=$5 
DATA_SET_NAME=$6

module load bedtools/2.21.0

# convert homer .pos file to a .bed file
# 1) remove header from pos file
# 2) use awk to get chrom (pos col 2), start (pos col 3), end (pos col 4) 
#    name (pos col 1), peak score (col 6), then in the bed file 6 to 12 are "." for empty
#	 then bed file col 13 is Annotation (pos col 8), bed file col 14 is Detail Annotation (pos col 9)
# 3) 

tmp1_bed="$OUTPUT_DATA_FOLDER"/closest_peak_tmp1.bed

tail -n +2 $PEAK_POS_FILE | awk -F"\t" -v OFS="\t" '{print $2, $3, $4, $1, $6, $8, $9}' > $tmp1_bed

tmp2_bed="$OUTPUT_DATA_FOLDER"/closest_peak_tmp2.bed

#
# use bedtools closest to find the closest H3K4ME3 peak to each exon following the rules below
# 1) if the exon is on the "-" strand the closest peak is the peak on the 3` end
# 2) if the exon is on the "+" strand the closest peak is the peak on the 5` end
# options used...
# -a is the file with exons
# -b is the file with the H3K4ME3 peaks
# -D a include distance in output, upstream features are reported as negative, 
# and use strand of feature in file a to determine what is upstream and downstream
# -id ignore downstream
bedtools closest -a $FILE_TO_UPDATE -b $tmp1_bed -t first -D a -id > $tmp2_bed

mv $tmp2_bed $FILE_TO_UPDATE
rm $tmp2_bed
rm $tmp1_bed

# construct header
echo -e "peak_chrom\tpeak_start\tpeak_end\tpeak_id\tpeak_score\tpeak_annotation\tpeak_annotation_detail" > $tmp1_bed
paste $HDR_TO_UPDATE $tmp1_bed > $tmp2_bed && mv $tmp2_bed $HDR_TO_UPDATE

rm $tmp2_bed
rm $tmp1_bed

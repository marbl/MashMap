#!/bin/bash

FASTA_FAI=$1
PAF=$2
COVERAGE=$3

cat $FASTA_FAI | awk -v OFS='\t' '{print($1,"0",$2)}' > $PAF.sequences.bed
cat \
    <(cat $PAF | awk -v OFS='\t' '{print $1, $3, $4, "", "", $5}') \
    <(cat $PAF | awk -v OFS='\t' '{print $6, $8, $9, "", "", "+"}') \
    | bedtools sort | bedtools merge > $PAF.query+target.bed

echo "#seq.name" coverage | tr ' ' '\t' > $PAF.coverage.txt
bedtools intersect -a $PAF.sequences.bed -b $PAF.query+target.bed -wo > $PAF.overlap.bed
awk 'BEGIN{FS=OFS="\t"}{
    if(NR==FNR){
        len[$1]=$3-$2; coverage[$1]=0; 
    } else {
        coverage[$1]+=$NF
    }
    } END{
    for(seq in len){
        printf("%s\t%f\n", seq, coverage[seq] / len[seq])
    }
}' $PAF.sequences.bed $PAF.overlap.bed >> $PAF.coverage.txt

cat \
    <(head -n 1 $PAF.coverage.txt) \
    <(sed '1d' $PAF.coverage.txt | sort -k 2,2nr -k 1,1) | column -t

awk -v threshold=$COVERAGE 'NR > 1 && $2 < threshold {
    print "Low coverage for sequence " $1 " with coverage " $2;
    flag = 1
} END {
    if (flag) exit 1
}' $PAF.coverage.txt

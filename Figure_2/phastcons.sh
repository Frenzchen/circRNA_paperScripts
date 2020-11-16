#!/usr/bin bash


# get phastcons for circrna splice sites

bigwig="/archive/cig/kaessmann/fgruhl/genomes/mouse/mm10/mm10_GRCm38/mm10.60way.phastCons60wayPlacental.bw"

while IFS=$'\t' read -r -a line; do
	~/bin/bigWigToBedGraph -chrom=${line[0]} -start=${line[1]} -end=${line[2]} $bigwig  tmp.bed;
	
	cat tmp.bed >> whole.bed;
done < "AccDon_exons.bed"
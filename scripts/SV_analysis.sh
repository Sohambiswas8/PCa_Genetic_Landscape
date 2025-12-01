#!/bin/bash
# a bash script for identifying gnomAD annotated structural variants

zcat gnomad.v4.1.sv.sites.bed.gz | cut -f1,2,3,5 | awk 'BEGIN{OFS="\t"} $1!~/^#/ {print $1,$2,$3,$4}' | sort -k1,1V -k2,2n > gnomad4.sv.bed

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' WES.vt_indels.vcf.gz \
  | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $2, $2, $3, $4, $5}' \
  | sort -k1,1V -k2,2n \
  | bedtools intersect -a - -b gnomad4.sv.bed -wa -wb \
  | awk 'BEGIN{OFS="\t"; print "CHROM","POS","ID","REF","ALT","END","SVTYPE"} {print $1,$4,$5,$6,$7,$9,$10}' \
  > matched_vcf_with_gnomad_svinfo.tsv
mv matched_vcf_with_gnomad_svinfo.tsv WES.vt_indels.point.bed 
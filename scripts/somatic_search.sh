#!/bin/bash
# a search strategy to identify prostate cancer (PCa) somatic SNVs from variant calls
# download the Somatic variant list from ClinVar (https://www.ncbi.nlm.nih.gov/clinvar) by the keyword "Prostate Cancer"
# the ClinVar somatic list ----> clinvar_somatic_variants.txt

cut -f8,9,12,16 clinvar_somatic_variants.txt \
     > clinvar_somatic_variants_dbsnp.txt     # Extract only Chromosome (GRCh38Chromsome), position (GRCh38Location), rsID (dbSNP ID) and clinvar classification

# the output might contain some range for short variants
# the next section demarkate the start and end positions for small variants

########################################################################### (FORMAT) #############################################################################

awk 'BEGIN{FS=OFS="\t"}
NR==1 {
  print "Chr","start","End","ID"
  next
}
{
  # read fields (robust to missing $3)
  chr = $1
  loc = $2
  id  = ($3=="" ? "." : $3)

  # trim whitespace
  gsub(/^[ \t]+|[ \t]+$/, "", chr)
  gsub(/^[ \t]+|[ \t]+$/, "", loc)
  gsub(/^[ \t]+|[ \t]+$/, "", id)

  # extract numbers from loc; handles "12345", "12345-12346", "12345 - 12346", etc.
  if (loc == "" ) {
    start = ""; end = ""
  } else {
    n = split(loc, a, /[^0-9]+/)
    # the split above may create empty elements; find first two numeric tokens
    start=""; end=""
    for(i=1;i<=n;i++){
      if(a[i] ~ /^[0-9]+$/) {
        if(start=="") start=a[i]
        else if(end=="") end=a[i]
      }
    }
    if(start=="" ) start = loc      # fallback: keep original text
    if(end=="" ) end = start
  }

  print chr, start, end, id
}' clinvar_somatic_variants_dbsnp.txt > clinvar_somatic_chr_start_end_id.txt

# FINAL: search for matched somatic records from the formatted clinvar list to the target VCF file
# search key: chr:pos:id
# if the rsID is missing in the target VCF then automatically reverts back to chr:pos and record all possible REF and ALT alleles in VCF format
#######################################################################################################################################################################

awk -F'\t' 'BEGIN { OFS = FS = "\t" }
NR==FNR {
    # Skip header in TSV if present (assumes header label "Chr"; if not add "chr")
    if (FNR==1 && $1 ~ /^Chr$/) next
    chr = $1; sub(/^chr/, "", chr)      # normalize TSV chr (strip leading "chr")
    pos = $2
    rsid = ($4 == "" ? "." : $4)       # use "." for missing RSID
    clinmap[chr "_" pos] = rsid
    next
}
# Now processing VCF file
/^##/ { print; next }                  # print meta-lines unchanged
/^#CHROM/ { print; next }              # print header line once
{
    chr = $1
    tmp = chr; sub(/^chr/, "", tmp)    # normalize VCF chr for matching
    key = tmp "_" $2

    # Only print variants that are present in the TSV (match by chr:pos)
    if (key in clinmap) {
        if (clinmap[key] != ".") {
            $3 = clinmap[key]         # populate ID only if TSV has rsID
        }
        print
    }
    # else skip (do not print unmatched variant lines)
}' clinvar_somatic_chr_start_end_id.tsv WES_LPS.vt.vcf > WES_LPS_somatic_vt.vcf
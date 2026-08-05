#!/usr/bin/env bash
set -euo pipefail

# ====================================================================================
# Configuration – adjust the paths to the environment
# ====================================================================================
READS_DIR="reads"                             # directory containing *_R{1,2}.fastq.gz
QC_DIR="fastqc_raw"                           # QC results will be placed here
TRIMMED_DIR="reads_trim"                      # trimmed reads will be placed here
TRIMMED_QC_DIR="fastqc_trim"                  # QC of trimmed reads will be placed here
TRIMMED_REPORT_DIR="trim_report"              # fastp reports will be placed here
SAM_DIR="sam"                                 # SAM files will be placed here
BAM_DIR="bam"                                 # BAM files will be placed here
SORT_BAM_DIR="sort_bam"                       # sorted BAM files will be placed here
MPILEUP_DIR="mpileup"                         # mpileup files will be placed here
VARSCAN_DIR="varscan"
BCFTOOLS_VCF_DIR="bcftools_vcf"
FREEBAYES_VCF_DIR="freebayes_vcf"             # VarScan results will be placed here
VT_VCF_DIR="vt_vcf"                           # VarScan results will be placed here
CONSENSUS_VCF_DIR="consensus_vcf"             # VarScan results will be

REFERENCE="/home/Data/hg38/hg38.fa"           # reference genome (with .fai and .dict)
BOWTIE2_INDEX="/home/Data/hg38/hg38"          # BOWTIE2 index path
THREADS=20                                    # number of CPU threads to use

# =====================================================================================
# Timer – record start time and print elapsed on exit
# =====================================================================================
start_time=$(date +%s)

finish() {
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    hours=$((elapsed / 3600))
    minutes=$(( (elapsed % 3600) / 60 ))
    seconds=$((elapsed % 60))
    printf "\nTotal pipeline execution time: %02d:%02d:%02d\n" "$hours" "$minutes" "$seconds"
}
trap finish EXIT

# ============================================================================
# Create required subdirectories
# ============================================================================
mkdir -p fastqc_raw fastqc_trim reads_trim sam bam sort_bam mpileup varscan bcf bcftools_vcf freebayes_vcf vt_vcf consensus_vcf clinvar_results

# ============================================================================
# Loop over all R1 files and process each sample
# ============================================================================
for r1 in "${READS_DIR}"/*_1.fastq.gz; do
    # Extract sample name (e.g., S16 from S16_R1.fastq.gz)
    sample=$(basename "$r1" _1.fastq.gz)
    r2="${READS_DIR}/${sample}_2.fastq.gz"

    if [[ ! -f "$r2" ]]; then
        echo "ERROR: Missing _2 file for sample $sample" >&2 
        exit 1
    fi
    
    echo "=== Processing sample: $sample ==="

    # --------------------------------------------------------------------------------
    # 1. Quality check of raw reads
    # --------------------------------------------------------------------------------
    fastqc -o "${QC_DIR}" -t "${THREADS}" "$r1" "$r2"
    echo "Quality check completed for sample: $sample"

    # --------------------------------------------------------------------------------
    # 2. Trim adapters and low-quality bases with fastp
    # --------------------------------------------------------------------------------
    fastp --in1 "$r1" --in2 "$r2" \
          --out1 "${TRIMMED_DIR}/${sample}_1.trimmed.fastq.gz" \
          --out2 "${TRIMMED_DIR}/${sample}_2.trimmed.fastq.gz" \
          --detect_adapter_for_pe \
          --thread 16 \
          --verbose \
          -R "Sample ${sample} fastp report" \
          -h "${TRIMMED_REPORT_DIR}/${sample}_fastp.html" \
          -j "${TRIMMED_REPORT_DIR}/${sample}_fastp.json"
    echo "Trimming completed for sample: $sample"

    # --------------------------------------------------------------------------------
    # 3. Quality check of trimmed reads
    # --------------------------------------------------------------------------------
    fastqc -o "${TRIMMED_QC_DIR}" -t "${THREADS}" \
           "${TRIMMED_DIR}/${sample}_1.trimmed.fastq.gz" \
           "${TRIMMED_DIR}/${sample}_2.trimmed.fastq.gz"
    echo "Quality check of trimmed reads completed for sample: $sample"

    # --------------------------------------------------------------------------------
    # 4. Align with bwa-mem2, convert to BAM, and sort (all in a pipe)
    # --------------------------------------------------------------------------------
    bowtie2 -x ${BOWTIE2_INDEX} -1 ${READS_DIR}/${sample}_1.fastq.gz -2 ${READS_DIR}/${sample}_2.fastq.gz --rg-id "00${sample}" --rg "SM:${sample}" --rg "PL:ILLUMINA" \
    --threads ${THREADS} -t | \
    samtools view -bS -@ "${THREADS}" -o "${BAM_DIR}/${sample}.bam"
    samtools sort -@ "${THREADS}" -o "${SORT_BAM_DIR}/${sample}.sorted.bam" "${BAM_DIR}/${sample}.bam"
    echo "Alignment and sorting completed for sample: $sample"

    # ------------------------------------------------------------------------------
    # 5. Index the sorted BAM
    # ------------------------------------------------------------------------------
    samtools index "${SORT_BAM_DIR}/${sample}.sorted.bam"
    echo "Indexing completed for sample: $sample"

    # =============================
    # 4) VarScan (SNP & INDEL)
    # =============================
    echo "-------------------------------"
    echo "[${sample}] VarScan SNP call..."
    echo "-------------------------------"

    echo "Generating mpileup for sample: $sample and calling SNPs with VarScan..."
    samtools mpileup -f "${REFERENCE}" "${SORT_BAM_DIR}/${sample}.sorted.bam" | \
        varscan mpileup2snp \
            --min-coverage 8 --min-reads2 2 --min-var-freq 0.01 --min-freq-for-hom 0.75 \
            --strand-filter 1 --min-avg-qual 15 --p-value 0.05 --output-vcf 1 > varscan/${sample}.varscan.snp.vcf || echo "VarScan SNP failed for $sample"
  
    echo "---------------------------------"
    echo "[${sample}] VarScan INDEL call..."
    echo "---------------------------------"
    echo "Generating mpileup for sample: $sample and calling INDELs with VarScan..."
    samtools mpileup -f "${REFERENCE}" "${SORT_BAM_DIR}/${sample}.sorted.bam" | \
        varscan mpileup2indel \
            --min-coverage 8 --min-reads2 2 --min-var-freq 0.01 --min-freq-for-hom 0.75 \
            --strand-filter 1 --min-avg-qual 15 --p-value 0.05 --output-vcf 1 > varscan/${sample}.varscan.indel.vcf || echo "VarScan INDEL failed for $sample"
    
    # --------------------------------------------------------------------------------
    # 7. BCFtools variant calling
    # --------------------------------------------------------------------------------
    echo "-------------------------------"
    echo "[${sample}] bcftools calling..."
    echo "-------------------------------"
    bcftools mpileup -Ou -f $REFERENCE sort_bam/${sample}.sorted.bam | \
        bcftools call -mv -Ob -o bcf/${sample}.bcftools.bcf || echo "bcftools mpileup/call failed"
    bcftools view bcf/${sample}.bcftools.bcf > bcftools_vcf/${sample}.bcftools.vcf || echo "bcftools view failed"

    # --------------------------------------------------------------------------------
    # 8. FreeBayes variant calling
    # --------------------------------------------------------------------------------
    echo "------------------------"
    echo "[${sample}] freebayes..."
    echo "------------------------"
    /home/software/freebayes-1.3.10-linux-amd64-static -f $REFERENCE sort_bam/${sample}.sorted.bam > freebayes_vcf/${sample}.freebayes.vcf || echo "freebayes failed"

    # --------------------------------------------------------------------------------
    # 9. VT variant calling
    # --------------------------------------------------------------------------------
    echo "--------------------------"
    echo "[${sample}] VT discover..."
    echo "--------------------------"
    vt discover -b sort_bam/${sample}.sorted.bam -s ${sample} -r $REFERENCE -o vt_vcf/${sample}.vt.vcf || echo "vt discover failed"

    echo "=== Finished sample: $sample ==="


    # ---------------------------------------------------------------------------------
    # 10. Sort Variants
    # ---------------------------------------------------------------------------------
    echo "--------------------------------------------"
    echo "[${sample}] Sorting VCFs before consensus..."
    echo "--------------------------------------------"
    echo "[${sample}] varscan sort..."
    awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' varscan/${sample}.varscan.snp.vcf > varscan/${sample}.varscan.sorted.vcf || echo "varscan sort failed"
    echo "[${sample}] bcftools sort..."
    awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' bcftools_vcf/${sample}.bcftools.vcf > bcftools_vcf/${sample}.bcftools.sorted.vcf || echo "bcftool sort failed"
    echo "[${sample}] freebayes sort..."
    awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' freebayes_vcf/${sample}.freebayes.vcf > freebayes_vcf/${sample}.freebayes.sorted.vcf || echo "bcftool sort failed"
    echo "[${sample}] vt sort..."
    awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' vt_vcf/${sample}.vt.vcf > vt_vcf/${sample}.vt.sorted.vcf || echo "vt sort failed"

    # ----------------------------------------------------------------------------------
    # 11. Consensus Calling
    # ----------------------------------------------------------------------------------
    echo "-----------------"
    echo "Consensus Calling"
    echo "-----------------"
    echo "-------------------------------------"
    echo "Count of varscan SNVs:'"
    awk '!/^#/' varscan/${sample}.varscan.snp.vcf | wc -l
    echo "Count of varscan indels:'"
    awk '!/^#/' varscan/${sample}.varscan.indel.vcf | wc -l
    echo 'Count of bcftools variants:'
    awk '!/^#/' bcftools_vcf/${sample}.bcftools.vcf | wc -l
    echo "Count of freebayes variants"
    awk '!/^#/' freebayes_vcf/${sample}.freebayes.vcf | wc -l
    echo 'Count of vt variants:'
    awk '!/^#/' vt_vcf/${sample}.vt.vcf | wc -l
    echo "-------------------------------------"

    echo "-------------------------------------------"
    echo "Freebayes and bcftools Consensus"
    echo "-------------------------------------------"
    awk 'NR == FNR {a[$1 $2 $4 $5];next} $1 $2 $4 $5 in a' freebayes_vcf/${sample}.freebayes.vcf bcftools_vcf/${sample}.bcftools.vcf > consensus_vcf/${sample}.freebayes_bcftools.vcf
    echo "Count of freebayes and bcftools consensus variants:'"
    awk '!/^#/' consensus_vcf/${sample}.freebayes_bcftools.vcf | wc -l
    echo "-------------------------------------------"
    echo "Freebayes, bcftools and vt Consensus"
    echo "-------------------------------------------"
    awk 'NR == FNR {a[$1 $2 $4 $5];next} $1 $2 $4 $5 in a' consensus_vcf/${sample}.freebayes_bcftools.vcf vt_vcf/${sample}.vt.vcf > consensus_vcf/${sample}.freebayes_bcftools_vt.vcf
    echo "Count of freebayes, bcftools and vt consensus variants:'"
    awk '!/^#/' consensus_vcf/${sample}.freebayes_bcftools_vt.vcf | wc -l
    echo "-------------------------------------------"
    echo "Freebayes, bcftools, vt and varscan snp Consensus"
    echo "-------------------------------------------"
    awk 'NR == FNR {a[$1 $2 $4 $5];next} $1 $2 $4 $5 in a' consensus_vcf/${sample}.freebayes_bcftools_vt.vcf varscan/${sample}.varscan.snp.vcf > consensus_vcf/${sample}.freebayes_bcftools_vt_varsnp.vcf
    echo "Count of Freebayes, bcftools, vt and varscan snp  consensus variants:'"
    awk '!/^#/' consensus_vcf/${sample}.freebayes_bcftools_vt_varsnp.vcf | wc -l
    echo "---------------------------------------------------------"
    echo "Freebayes, bcftools, vt and varscan indel Consensus"
    echo "-------------------------------------------"
    awk 'NR == FNR {a[$1 $2 $4 $5];next} $1 $2 $4 $5 in a' consensus_vcf/${sample}.freebayes_bcftools_vt.vcf varscan/${sample}.varscan.indel.vcf > consensus_vcf/${sample}.freebayes_bcftools_vt_varind.vcf
    echo "Count of Freebayes, bcftools, vt and varscan indel  consensus variants:'"
    awk '!/^#/' consensus_vcf/${sample}.freebayes_bcftools_vt_varind.vcf | wc -l
    echo "---------------------------------------------------------"

    # ----------------------------------------------------------------------------------
    # 12. ClinVar validation
    # ----------------------------------------------------------------------------------

    echo "-----------------------------------------------------------------"
    echo "Finding valid Clinvar variant id's corresponding to the consensus"
    echo "-----------------------------------------------------------------"

    awk 'NR == FNR {a[$1 $2 $4 $5];next} (("chr")$1 $2 $4 $5 in a)' consensus_vcf/${sample}.freebayes_bcftools_vt_varsnp.vcf /home/Data/clinvar.vcf > clinvar_results/${sample}.clinvar_validated_variants.vcf
    echo "------------------------------------------------------------------------------------"
    echo "Total number of validated Clinvar entries corresponding to your list of variants is:"
    wc -l clinvar_results/${sample}.clinvar_validated_variants.vcf
    echo "------------------------------------------------------------------------------------"
    echo "--------------------------------"
    echo "Downstream Pathgenicity Analysis"
    echo "--------------------------------"
    grep -w "CLNSIG" clinvar_results/${sample}.clinvar_validated_variants.vcf > clinvar_results/${sample}.clnid_available.vcf
    echo "----------------------------------------------------------------"
    echo "Total number of validated Clinvar entries with CLNSIG available:"
    wc -l clinvar_results/${sample}.clnid_available.vcf
    echo "----------------------------------------------------------------"
 
    grep -w "CLNSIG=Pathogenic" clinvar_results/${sample}.clnid_available.vcf > clinvar_results/${sample}.pathogenic.vcf
    echo "------------------------------------------"
    echo "Number of pathogenic variants in ${sample}"
    wc -l clinvar_results/${sample}.pathogenic.vcf
    echo "------------------------------------------"
    grep -w "CLNSIG=Likely_pathogenic" clinvar_results/${sample}.clnid_available.vcf > clinvar_results/${sample}.likely_pathogenic_variants.vcf
    echo "------------------------------------"
    echo "Number of likely pathogenic variants"
    wc -l clinvar_results/${sample}.likely_pathogenic_variants.vcf
    echo "------------------------------------"
    grep -w "CLNSIG=Pathogenic/Likely_pathogenic" clinvar_results/${sample}.clnid.vcf > clinvar_results/${sample}.pathogenic_and_likely_pathogenic.vcf
    echo "-----------------------------------------------"
    echo "Number of pathogenic_likely_pathogenic variants"
    wc -l clinvar_results/${sample}.pathogenic_and_likely_pathogenic.vcf
    echo "-----------------------------------------------"
    grep -w "CLNSIG=Uncertain_significance" clinvar_results/${sample}.clnid.vcf > clinvar_results/${sample}.uncertain_significance.vcf
    echo "--------------------------------------------"
    echo "Number of variants of uncertain significance"
    wc -l clinvar_results/${sample}.uncertain_significance.vcf
    echo "--------------------------------------------"
    grep -w "CLNSIG=Conflicting_classifications_of_pathogenicity" clinvar_results/${sample}.clnid.vcf > clinvar_results/${sample}.conflicting_variants.vcf
    echo "--------------------------------------------"
    echo "Number of variants of Conflicting classifications of pathogenicity"
    wc -l clinvar_results/${sample}.conflicting_variants.vcf
    echo "--------------------------------------------"
done

echo "All samples processed successfully."


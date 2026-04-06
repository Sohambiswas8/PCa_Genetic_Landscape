#!/usr/bin/env bash
set -u
shopt -s nullglob

# --- Editing maybe required to match the local system ---
REFERENCE="ref/hg38.fa"
BOWTIE2_INDEX="ref/hg38"          # bowtie2 index basename
# -------------------------------------

# tools required: fastqc, trim_galore, bowtie2, samtools, varscan, bcftools, vt, freebayes
# create directories
# mkdir -p sam bam sort_bam mpileup varscan bcftools_vcf bcf vt_vcf freebayes_vcf \
#         fastqc_raw reads_trim fastqc_trim

# sample list: LPS1..LPS9 plus LPSUK1
samples=( LPS{2..9} LPSUK1 )

start_all=$SECONDS

for sample in "${samples[@]}"; do
  echo
  echo "==== Processing sample: $sample ===="
  tstart=$SECONDS

  r1="reads/${sample}_R1.fastq.gz"
  r2="reads/${sample}_R2.fastq.gz"

  # check input fastqs exist
  if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "WARNING: input FASTQ missing for ${sample}. Expected files:"
    echo "  $r1"
    echo "  $r2"
    echo "Skipping $sample"
    continue
  fi

  # ----------------------------
  # 1) FastQC on raw reads
  # ----------------------------
  echo "[${sample}] Running FastQC (raw reads)..."
  fastqc ${r1} ${r2} -t 10 -o fastqc_raw || echo "FastQC (raw) failed for ${sample} (continuing)"

  # ----------------------------
  # 2) Fastp (paired-end trimming) + FastQC on trimmed reads
  # ----------------------------
  echo "[${sample}] Running Fastp (paired + post-trim FastQC)..."
  
  trimdir="reads_trim"
  # run Fastp paired mode; keeps output as <input_basename>_val_1.fq.gz and _val_2.fq.gz
  /home/soham/fastp --in1 ${r1} --in2 ${r2} --out1 ${trimdir}/${sample}_R1_trimmed.fastq.gz --out2 ${trimdir}/${sample}_R2_trimmed.fastq.gz --detect_adapter_for_pe \
    --thread 8 -j ${trimdir}/${sample}_fastp_report.json -h ${trimdir}/${sample}_fastp_report.html  || {
    echo "Fastp failed for ${sample} — skipping to next sample"
    continue
  }

  # Locate trimmed files: Fastp names are usually: <sample>_R1_trimmed.fastq.gz and _R2_trimmed.fastq.gz
  # We'll pick the first match that contains sample and _R1/_R2.
  trimmed_r1=$(ls -1 "${trimdir}/${sample}"*R1*trimmed* 2>/dev/null | head -n1 || true)
  trimmed_r2=$(ls -1 "${trimdir}/${sample}"*R2*trimmed* 2>/dev/null | head -n1 || true)

  # fallback: try patterns without R1/R2 if above didn't match
  if [[ -z "$trimmed_r1" ]]; then
    trimmed_r1=$(ls -1 "${trimdir}"/*"${sample}"*val*1* 2>/dev/null | head -n1 || true)
  fi
  if [[ -z "$trimmed_r2" ]]; then
    trimmed_r2=$(ls -1 "${trimdir}"/*"${sample}"*val*2* 2>/dev/null | head -n1 || true)
  fi

  if [[ ! -f "$trimmed_r1" || ! -f "$trimmed_r2" ]]; then
    echo "ERROR: trimmed files not found for $sample. Expected in $trimdir. Skipping."
    ls -l "${trimdir}"/*"${sample}"* || true
    continue
  fi

  echo "[${sample}] trimmed files: ${trimmed_r1} , ${trimmed_r2}"
  echo "[${sample}] Running FastQC on trimmed reads..."

  fastqc ${trimmed_r1} ${trimmed_r2} -t 10 -o fastqc_trim || echo "FastQC (trimmed) failed for ${sample} (continuing)"

  echo "[${sample}] fastqc for trimmed reads done"
  # ----------------------------
  # 3) Align with bowtie2, sorting and indexing and pileup with samtools
  # ----------------------------
  echo "[${sample}] Aligning with bowtie2..."
  bowtie2 -x ${BOWTIE2_INDEX} -1 ${trimmed_r1} -2 ${trimmed_r2} --rg-id "00${sample}" --rg "SM:${sample}" --rg "PL:ILLUMINA" \
    --threads 20 -t -S "sam/${sample}.sam" || {
    echo "bowtie2 failed for ${sample} — skipping"
    continue
  }
  echo "[${sample}] bowtie2 done"

  # convert SAM -> BAM
  echo "[${sample}] Converting SAM to BAM..."
  samtools view -bS sam/${sample}.sam -o bam/${sample}.bam || { echo "samtools view failed"; continue; }
  echo "[${sample}] bam created: bam/${sample}.bam"

  # sort & index
  echo "[${sample}] Sorting and indexing BAM..."
  samtools sort bam/${sample}.bam -o sort_bam/${sample}.sorted.bam || { echo "samtools sort failed"; continue; }
  samtools index sort_bam/${sample}.sorted.bam || { echo "samtools index failed"; continue; }
  echo "[${sample}] sorted & indexed"

  # mpileup (create pileup file used by VarScan)
  echo "[${sample}] Creating mpileup..."
  samtools mpileup -B -q 1 -f $REFERENCE sort_bam/${sample}.sorted.bam > mpileup/${sample}.pileup || { echo "mpileup failed"; continue; }
  echo "[${sample}] mpileup created"

  # ----------------------------
  # 4) VarScan (SNP & INDEL)
  # ----------------------------

  echo "[${sample}] VarScan SNP call..."

  # mpileup2snp reads mpileup text; make sure VarScan jar is in PATH or varscan is callable
  varscan mpileup2snp mpileup/${sample}.pileup \
    --min-coverage 8 --min-reads2 2 --min-var-freq 0.01 --min-freq-for-hom 0.75 \
    --strand-filter 1 --min-avg-qual 15 --p-value 0.01 --output-vcf 1 > varscan/${sample}.varscan.snp.vcf || echo "VarScan SNP failed for $sample"

  echo "[${sample}] VarScan INDEL call..."
  varscan mpileup2indel mpileup/${sample}.pileup \
    --min-coverage 8 --min-reads2 2 --min-var-freq 0.01 --min-freq-for-hom 0.75 \
    --strand-filter 1 --min-avg-qual 15 --p-value 0.01 --output-vcf 1 > varscan/${sample}.varscan.indel.vcf || echo "VarScan INDEL failed for $sample"

  # ----------------------------
  # 5) Bcftools
  # ----------------------------

  echo "[${sample}] bcftools calling..."
  bcftools mpileup -Ou -f $REFERENCE sort_bam/${sample}.sorted.bam | \
    bcftools call -mv -Ob -o bcf/${sample}.bcftools.bcf || echo "bcftools mpileup/call failed"
  bcftools view bcf/${sample}.bcftools.bcf > bcftools_vcf/${sample}.bcftools.vcf || echo "bcftools view failed"

  
  # ----------------------------
  # 6) Freebayes
  # ----------------------------

  echo "[${sample}] freebayes..."
  /home/soham/freebayes-1.3.10-linux-amd64-static -f $REFERENCE sort_bam/${sample}.sorted.bam > freebayes_vcf/${sample}.freebayes.vcf || echo "freebayes failed"


  # ----------------------------
  # 7) VT
  # ----------------------------

  # echo "[${sample}] VT discover..."
  /home/software/vt/vt discover -b sort_bam/${sample}.sorted.bam -s ${sample} -r $REFERENCE -o vt_vcf/${sample}.vt.vcf || echo "vt discover failed"


  # Done with sample
  tend=$SECONDS
  took=$((tend - tstart))
  echo "[${sample}] Completed in $((took/60))m $((took%60))s"

done

total_end=$SECONDS
elapsed=$(( total_end - start_all ))
echo "All samples done. Total elapsed time: $((elapsed/60))m $((elapsed%60))s"


#!/bin/bash
# Genotype check for WANG2024_HEAT_GILL: are the 6 C. gigas RNA-seq libraries
# 6 individuals, technical replicates of fewer RNA samples, or pools?
#
# Run from the repository root of a checkout that has the comparison's Salmon
# output under results/ and the FASTQs under data/raw/ (as run on Pori,
# 2026-10-13, with minimap2 2.31, samtools 1.24, bcftools 1.24 from the aree
# conda env). Intermediate BAMs and VCFs go to a scratch directory; the two
# summary tables in this folder were produced from its calls.vcf.gz.
#   bash data/studies/WANG2024_HEAT_GILL/genotype_check/run_genocheck.sh [scratch_dir]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WT=$(cd "$HERE/../../../.." && pwd)
D=${1:-$(mktemp -d)}
SAL=$WT/results/rnaseq/WANG2024_HEAT_GILL/gigas_heat_12h_vs_control/rnaseq/salmon
RAW=$WT/data/raw/WANG2024_HEAT_GILL
REF=$WT/data/reference/GCF_963853765.1
NPAIRS=3000000
mkdir -p "$D"; cp "$HERE"/*.py "$D"/
cd "$D"

# 1. Top 400 genes by mean TPM, one transcript (the most expressed) per gene.
python3 "$D/pick_transcripts.py" "$SAL" "$REF/tx2gene.tsv" 400 > targets.txt
samtools faidx "$REF/rna.fna" 2>/dev/null || true
samtools faidx "$REF/rna.fna" -r targets.txt > targets.fa
samtools faidx targets.fa
echo "targets: $(grep -c '>' targets.fa) transcripts, $(awk '{s+=$2} END{print s}' targets.fa.fai) bp"

# 2. Map the first NPAIRS read pairs of each library.
for r in SRR24489461 SRR24489463 SRR24489465 SRR24489471 SRR24489481 SRR24489483; do
  minimap2 -t 6 -ax sr --secondary=no -R "@RG\tID:$r\tSM:$r" targets.fa \
    <(gzip -cd "$RAW/${r}_1.fastq.gz" | head -n $((NPAIRS*4))) \
    <(gzip -cd "$RAW/${r}_2.fastq.gz" | head -n $((NPAIRS*4))) 2>/dev/null |
    samtools view -b -q 20 -F 0x904 - | samtools sort -o "$r.bam" -
  samtools index "$r.bam"
  echo "$r mapped: $(samtools view -c "$r.bam")"
done

# 3. Joint variant calling with per-sample allele depths.
bcftools mpileup -f targets.fa -a AD,DP -d 100000 -q 20 -Q 20 -Ou SRR*.bam 2>/dev/null |
  bcftools call -mv -Oz -o calls.vcf.gz 2>/dev/null
bcftools index calls.vcf.gz
bcftools view -v snps -m2 -M2 calls.vcf.gz |
  bcftools query -H -f '%CHROM\t%POS[\t%AD]\n' > snps_ad.tsv
echo "biallelic SNPs called: $(($(wc -l < snps_ad.tsv) - 1))"

# 4. Compare allele fractions.
python3 "$D/compare_af.py" snps_ad.tsv "$WT/data/studies/WANG2024_HEAT_GILL/samplesheet.tsv" "$HERE"

Code to execute (assuming working environment with aws cli, R with bambu install, minimap2 and SAMtools):
set -euo pipefail
# Clone repository
git clone https://github.com/Rajvinder-d/dsa4262_task5.git
cd dsa4262_task5
mkdir -p inputs/fastq inputs/reference
# Download the four FASTQs
drna1=SGNex_A549_directRNA_replicate1_run1
drna2=SGNex_Hct116_directRNA_replicate6_run1
cdna1=SGNex_Hct116_cDNA_replicate1_run6
cdna2=SGNex_K562_cDNA_replicate1_run3

for sample in "$drna1" "$drna2" "$cdna1" "$cdna2"; do
    aws s3 cp --no-sign-request \
        "s3://sg-nex-data/data/sequencing_data_ont/fastq/$sample/$sample.fastq.gz" \
        "inputs/fastq/$sample.fastq.gz"
done

# Download the genome FASTA and annotation GTF
ensembl=https://ftp.ensembl.org/pub/release-91
fasta=inputs/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gtf=inputs/reference/Homo_sapiens.GRCh38.91.gtf

curl -fL --retry 3 \
    "$ensembl/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
    -o "$fasta"

curl -fL --retry 3 \
    "$ensembl/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz" \
    -o "$gtf.gz"

gzip -dc "$gtf.gz" > "$gtf"
# Scenario 1: annotated, from scratch, including indexing
nextflow run ./main.nf \
    --fasta "$fasta" \
    --gtf "$gtf" \
    --direct_rna "inputs/fastq/$drna1.fastq.gz,inputs/fastq/$drna2.fastq.gz" \
    --cdna "inputs/fastq/$cdna1.fastq.gz,inputs/fastq/$cdna2.fastq.gz" \
    --outdir ./results \
    -ansi-log false

test -s results/scenarios/annotated/execution/report.html

# Scenario 2: unannotated, reusing indexing/alignment/BAM conversion
nextflow run ./main.nf \
    --outdir ./results \
    --without_annotations \
    -resume \
    -ansi-log false

test -s results/scenarios/unannotated/execution/report.html

# Compare native Nextflow execution reports
python3 compare_reports.py \
    results/scenarios/annotated/execution/report.html \
    results/scenarios/unannotated/execution/report.html \
    | tee results/runtime_comparison.txt

#!/usr/bin/env bash
#SBATCH --job-name=bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=hailiu@sund.ku.dk
#SBATCH --nodes=1                          # Run all processes on a single node
#SBATCH --ntasks=1                         # Run a single task
#SBATCH --cpus-per-task=10                 # Number of CPU cores per task
#SBATCH --mem=32gb                         # Job memory request
#SBATCH --time=0-20:00:00
#SBATCH --output=DRB_TTchem_seq2_log%j.log

set -euo pipefail

########################################
### Author: Haiyue Liu
### Date: Initial version 15-12-2023
### Date: Bug fixed and pipeline optimized on 20-03-2026
### This is the standard pipeline for DRB/TTchem-seq2 data processing.
### Libraries: single-ended reversely-stranded reads with UMI (11-mer) sequences. We demultiplexed the UMI sequences as R2 which are attached to R1 as headers in our pipeline.
### We also use yeast spike-ins for normalization of the nascent RNA libraries. For this, we created a human-yeast combined genome reference for mapping.
### Main steps for data processing in this pipeline are:
### 1. Attach UMIs to reads header
### 2. Adapter trimming & fastqc
### 3. STAR alignment
### 4. Extract uniquely mapping reads
### 5. Deduplication (optional)
### 6. Removal of exon-intron-exon mapping reads
### 7. Split forward and reverse strand reads
### 8. Count mapped reads (library size)
### 9. Quantify gene counts
### 10. Calculate yeast spike-in size factors
### 11. Convert bam to bigwig format (spike-in normalized)
#############################################

#################################################
### CONFIGURATION
#################################################
### Configuration: this section is project-specific. You need to modify those lines based on your own project setup.
### Tools executables: all the tools are called via module load in this pipeline. If you use a different system, you could change those module load lines.
cores=10                                                    ### required
work_dir="/path/to/project/"                                ### required; project directory
script_dir="${work_dir}script/"                             ### required; script directory. The R function "size_factor.R" should be found there
sample_sheet="${work_dir}script/samplesheet.tsv"            ### required; one example tsv file can be found in this github repository
sequencing_type="SE"                                        ### required; options: "SE", "PE"
UMI=true                                                    ### required; options: true, false
UMI_sequence="NNNNNNNNNNN"                                  ### required if UMI=true; the length of UMIs can vary depends on the project
strandedness="reverse-stranded"                             ### required; options: "reverse-stranded", "stranded", "unstranded"
fastq_suffix=".fastq.gz"                                    ### required; options: ".fastq.gz", ".fq.gz"
index_dir="/path/to/star_index/"                            ### required; your STAR index directory (e.g. GRCh38.gencode.v43_sacCer3.108_star_2.7.9a_index)
gtf="/path/to/annotation.gtf"                               ### required; the annotation gtf file (e.g. GRCh38.gencode.v43_sacCer3.108.gtf)
chrsize="${index_dir}/chrNameLength.txt"                    ### required; the chromosome size file found in the STAR index directory
#################################################
#################################################

############################
### CHECK CONFIG PARAMETERS
############################
[ ! -d "${work_dir}" ]        && echo "ERROR: Working directory DOES NOT exist."        && exit 1
[ ! -f "${sample_sheet}" ]    && echo "ERROR: Samplesheet file DOES NOT exist."         && exit 1
[ ! -f "${script_dir}size_factor.R" ] && echo "ERROR: size_factor.R DOES NOT exist."   && exit 1  # FIX: was -d (directory test), must be -f (file test)
[ ! -d "${index_dir}" ]       && echo "ERROR: Index directory DOES NOT exist."          && exit 1
[ ! -f "${gtf}" ]             && echo "ERROR: Annotation gtf file DOES NOT exist."      && exit 1
[ ! -f "${chrsize}" ]         && echo "ERROR: Chromosome size file DOES NOT exist."     && exit 1

############################
### SAMPLE ID & NAME
############################
# FIX: use mapfile to populate true bash arrays (safe with any filename, correct [@] behaviour)
mapfile -t sample_ids   < <(awk 'NR>1 {print $2}' "${sample_sheet}")
mapfile -t sample_names < <(awk 'NR>1 {print $3}' "${sample_sheet}")

#############################
### DIRECTORIES
#############################
### As far as the work_dir is given, these are automatically created.
fastq_dir="${work_dir}fastq/"
fastqc_dir="${work_dir}fastqc/"
trimmed_fastq_dir="${work_dir}trimmed_fastq/"
tmp_dir="${work_dir}tmp/"
bam_dir="${work_dir}bam/"
featureCounts_dir="${work_dir}featureCounts/"
bigwig_dir="${work_dir}bigwig/"
analysis_dir="${work_dir}analysis/"
mkdir -p "${fastqc_dir}"
mkdir -p "${trimmed_fastq_dir}"
mkdir -p "${tmp_dir}"
mkdir -p "${bam_dir}"
mkdir -p "${featureCounts_dir}"
mkdir -p "${bigwig_dir}"
mkdir -p "${analysis_dir}"
cd "${work_dir}"

########## Data processing ###############

#########################################
### 1. Attach UMIs in R2 to header of R1
#########################################
module purge
module load dangpu_libs python/3.7.13 umi_tools/1.1.4

for sample_id in "${sample_ids[@]}"
do
  sample=$(awk -v s="${sample_id}" '$2==s {print $3}' "${sample_sheet}")  # FIX: removed useless cat | awk

  if [[ "${sequencing_type}" == "SE" ]]; then
    printf "Extract UMIs from R2 and attach them to R1 read headers for %s\n" "${sample}"  # FIX: echo "\n" -> printf
    R1_in="${fastq_dir}${sample_id}_R1_001${fastq_suffix}"
    umi_read="${fastq_dir}${sample_id}_R2_001${fastq_suffix}"
    R1_out="${fastq_dir}${sample}_R1_umi_attached${fastq_suffix}"
    umi_tools extract --extract-method=string --bc-pattern="${UMI_sequence}" \
      --stdin "${umi_read}" --read2-in="${R1_in}" --read2-out="${R1_out}"
  fi

  if [[ "${sequencing_type}" == "PE" ]]; then
    printf "Extract UMIs from R2 and attach them to R1 and R3 read headers for %s\n" "${sample}"  # FIX: echo "\n" -> printf
    R1_in="${fastq_dir}${sample_id}_R1_001.fastq.gz"
    R2_in="${fastq_dir}${sample_id}_R3_001.fastq.gz"
    umi_read="${fastq_dir}${sample_id}_R2_001.fastq.gz"
    R1_out="${fastq_dir}${sample}_R1_umi_attached.fastq.gz"
    R2_out="${fastq_dir}${sample}_R2_umi_attached.fastq.gz"
    umi_tools extract --extract-method=string --bc-pattern="${UMI_sequence}" \
      --stdin="${umi_read}" --read2-in="${R1_in}" --stdout="${R1_out}" --read2-stdout
    umi_tools extract --extract-method=string --bc-pattern="${UMI_sequence}" \
      --stdin="${umi_read}" --read2-in="${R2_in}" --stdout="${R2_out}" --read2-stdout
  fi
done

############################
### 2. QC & trim adapter
############################
module purge
module load anaconda3/2021.11
module load openjdk/13.0.1 perl/5.26.3 fastqc/0.11.9
module load pigz trimgalore/0.6.6
module load dangpu_libs python/3.7.13 cutadapt/4.1

if [[ "${UMI}" == true ]]; then name_suffix="umi_attached"; else name_suffix=""; fi

for sample in "${sample_names[@]}"
do
  echo "${sample} : FastQC -- Adapter trimming -- FastQC"
  if [[ "${sequencing_type}" == "SE" ]]; then
    fastqc -t "${cores}" -o "${fastqc_dir}" "${fastq_dir}${sample}_R1_${name_suffix}${fastq_suffix}"
    trim_galore --cores "${cores}" --basename "${sample}" \
      --output_dir "${trimmed_fastq_dir}" \
      --fastqc --fastqc_args "-o ${fastqc_dir} -t ${cores}" \
      "${fastq_dir}${sample}_R1_${name_suffix}${fastq_suffix}"
  fi
  if [[ "${sequencing_type}" == "PE" ]]; then
    fastqc -t "${cores}" -o "${fastqc_dir}" \
      "${fastq_dir}${sample}_R1_${name_suffix}${fastq_suffix}" \
      "${fastq_dir}${sample}_R2_${name_suffix}${fastq_suffix}"
    trim_galore --cores "${cores}" --paired --basename "${sample}" \
      --output_dir "${trimmed_fastq_dir}" \
      --fastqc --fastqc_args "-o ${fastqc_dir} -t ${cores}" \
      "${fastq_dir}${sample}_R1_${name_suffix}${fastq_suffix}" \
      "${fastq_dir}${sample}_R2_${name_suffix}${fastq_suffix}"
  fi
done

###########################
### 3. STAR alignment
###########################
module purge
module load gcc star/2.7.9a
module load samtools

rm -rf "${tmp_dir}"
mkdir -p "${tmp_dir}"

for sample in "${sample_names[@]}"
do
  echo "Mapping for ${sample}"
  # FIX: clear per-sample tmp dir before each run to avoid stale files
  rm -rf "${tmp_dir}${sample}_tmp"

  if [[ "${sequencing_type}" == "SE" ]]; then
    STAR \
      --runThreadN "${cores}" \
      --genomeDir "${index_dir}" \
      --readFilesIn "${trimmed_fastq_dir}${sample}_trimmed.fq.gz" \
      --readFilesCommand zcat \
      --outSAMattributes NH HI AS NM MD XS \
      --outSAMstrandField intronMotif \
      --outSAMmultNmax 1 \
      --outSAMunmapped None \
      --outFileNamePrefix "${bam_dir}${sample}." \
      --outSAMtype BAM SortedByCoordinate \
      --outBAMsortingBinsN 50 \
      --outTmpDir "${tmp_dir}${sample}_tmp"
  fi

  if [[ "${sequencing_type}" == "PE" ]]; then
    STAR \
      --runThreadN "${cores}" \
      --genomeDir "${index_dir}" \
      --readFilesIn "${trimmed_fastq_dir}${sample}_val_1.fq.gz" \
                    "${trimmed_fastq_dir}${sample}_val_2.fq.gz" \
      --readFilesCommand zcat \
      --outSAMattributes NH HI AS NM MD XS \
      --outSAMstrandField intronMotif \
      --outSAMmultNmax 1 \
      --outSAMunmapped None \
      --outFileNamePrefix "${bam_dir}${sample}." \
      --outSAMtype BAM SortedByCoordinate \
      --outBAMsortingBinsN 50 \
      --outTmpDir "${tmp_dir}${sample}_tmp"
      # FIX: was _val_1_trimmed.fq.gz / _val_2_trimmed.fq.gz — TrimGalore paired output is _val_1.fq.gz / _val_2.fq.gz
  fi

  samtools index -@ "${cores}" "${bam_dir}${sample}.Aligned.sortedByCoord.out.bam"
  rm "${bam_dir}${sample}.Log.out"
  rm "${bam_dir}${sample}.Log.progress.out"
done

########################################
### 4. Extract uniquely mapped reads
########################################
module purge
module load samtools/1.15.1

for sample in "${sample_names[@]}"
do
  echo "Extract unimappers for ${sample}"
  in_bam="${bam_dir}${sample}.Aligned.sortedByCoord.out.bam"
  out_bam="${bam_dir}${sample}.unimappers.bam"
  if [[ "${sequencing_type}" == "SE" ]]; then
    # FIX: removed -F 0x2 (excludes "properly paired" — flag is never set in SE data, misleading and unnecessary)
    samtools view --threads "${cores}" -q 255 "${in_bam}" -o "${out_bam}"
  fi
  if [[ "${sequencing_type}" == "PE" ]]; then
    samtools view --threads "${cores}" -q 255 -f 0x2 "${in_bam}" -o "${out_bam}"
  fi
  samtools index -@ "${cores}" "${out_bam}"
done

#######################################
### 5. Deduplication of UMIs (optional)
#######################################
module purge
module load parallel
module load dangpu_libs python/3.7.13 umi_tools/1.1.4
module load samtools/1.15.1

if [[ "${UMI}" == true ]]; then
  echo "Deduplication (running samples in parallel)"
  # Run all samples in parallel for speed; index each deduped bam afterwards
  find "${bam_dir}" -name "*.unimappers.bam" \
    | parallel -j "${cores}" \
        "umi_tools dedup --umi-separator='_' -I {} -S {.}.deduped.bam \
           --log={.}.dedup.log --output-stats={} && \
         samtools index {.}.deduped.bam"
fi

###############################################
### 6. Removal of exon-intron-exon reads
###############################################
module purge
module load samtools/1.15.1

for sample in "${sample_names[@]}"
do
  echo "Split unspliced and spliced reads for ${sample}"
  if [[ "${UMI}" == true ]]; then
    bam_name="${bam_dir}${sample}.unimappers.deduped"
  else
    bam_name="${bam_dir}${sample}.unimappers"
  fi
  samtools view --threads "${cores}" -he '[XS]'  "${bam_name}.bam" -o "${bam_dir}${sample}.spliced.bam"
  samtools view --threads "${cores}" -he '![XS]' "${bam_name}.bam" -o "${bam_dir}${sample}.unspliced.bam"
  samtools index "${bam_dir}${sample}.spliced.bam"
  samtools index "${bam_dir}${sample}.unspliced.bam"
done

###############################################################
### 7. Split reads transcribed from forward and reverse strands
###############################################################
module purge
module load samtools/1.15.1

for sample in "${sample_names[@]}"
do
  bam_name="${bam_dir}${sample}.unspliced"
  echo "Split forward and reverse strand ${sample}"

  if [[ "${sequencing_type}" == "SE" && "${strandedness}" == "reverse-stranded" ]]; then
    ### Forward strand reads -- R1 mapped to reverse strand (flag 16)
    samtools view -b -f 16  --threads "${cores}" "${bam_name}.bam" -o "${bam_name}.fwd.bam"
    samtools index -@ "${cores}" "${bam_name}.fwd.bam"
    ### Reverse strand reads -- mapped to forward strand (exclude flag 16)
    samtools view -b -F 16  --threads "${cores}" "${bam_name}.bam" -o "${bam_name}.rev.bam"
    samtools index -@ "${cores}" "${bam_name}.rev.bam"
  fi

  if [[ "${sequencing_type}" == "PE" && "${strandedness}" == "reverse-stranded" ]]; then
    ### Forward strand reads
    samtools view -b -f 80  --threads "${cores}" "${bam_name}.bam" -o "${bam_name}.R1.fwd.bam"
    samtools view -b -f 128 -F 16 --threads "${cores}" "${bam_name}.bam" -o "${bam_name}.R2.fwd.bam"
    samtools merge -l 9 --threads "${cores}" -f "${bam_name}.fwd.bam" "${bam_name}.R1.fwd.bam" "${bam_name}.R2.fwd.bam"
    samtools index -@ "${cores}" "${bam_name}.fwd.bam"
    ### Reverse strand reads
    samtools view -b -f 64  -F 16 --threads "${cores}" "${bam_name}.bam" -o "${bam_name}.R1.rev.bam"
    samtools view -b -f 144 --threads "${cores}" "${bam_name}.bam" -o "${bam_name}.R2.rev.bam"
    samtools merge -l 9 --threads "${cores}" "${bam_name}.rev.bam" "${bam_name}.R1.rev.bam" "${bam_name}.R2.rev.bam"
    samtools index -@ "${cores}" "${bam_name}.rev.bam"
    ### Remove intermediate files
    rm "${bam_name}.R1.fwd.bam" "${bam_name}.R2.fwd.bam"
    rm "${bam_name}.R1.rev.bam" "${bam_name}.R2.rev.bam"
  fi
done

#########################################
### 8. Count mapped reads (library size)
#########################################
module purge
module load samtools/1.15.1

reads_number="${work_dir}analysis/reads_number.txt"

if [[ "${UMI}" == true ]]; then
  printf '%s\t' 'sample_name' 'unimapper' 'deduped' 'unspliced' 'human_unspliced'
  printf '\n'
else
  printf '%s\t' 'sample_name' 'unimapper' 'unspliced' 'human_unspliced'
  printf '\n'
fi > "${reads_number}"

for sample in "${sample_names[@]}"
do
  echo "Count mapped reads for ${sample}"
  n_unimapper=$(samtools view --threads "${cores}" -c "${bam_dir}${sample}.unimappers.bam")
  if [[ "${UMI}" == true ]]; then
    bam_name="${bam_dir}${sample}.unimappers.deduped"
    n_deduped=$(samtools view --threads "${cores}" -c "${bam_name}.bam")
  else
    bam_name="${bam_dir}${sample}.unimappers"
  fi
  n_unspliced=$(samtools view --threads "${cores}" -c "${bam_dir}${sample}.unspliced.bam")
  n_human_unspliced=$(samtools idxstats "${bam_dir}${sample}.unspliced.bam" | awk '/^chr/ {s+=$3} END{print s}')

  if [[ "${UMI}" == true ]]; then
    echo "${sample}" | awk -v OFS="\t" \
      -v unimapper="${n_unimapper}" -v deduped="${n_deduped}" \
      -v unspliced="${n_unspliced}" -v human_unspliced="${n_human_unspliced}" \
      '{print $0, unimapper, deduped, unspliced, human_unspliced}' >> "${reads_number}"
  else
    echo "${sample}" | awk -v OFS="\t" \
      -v unimapper="${n_unimapper}" \
      -v unspliced="${n_unspliced}" -v human_unspliced="${n_human_unspliced}" \
      '{print $0, unimapper, unspliced, human_unspliced}' >> "${reads_number}"
  fi
done

#########################################
### 9. Quantify gene counts
#########################################
module purge
module load subread/2.0.3

for sample in "${sample_names[@]}"
do
  if [[ "${UMI}" == true ]]; then
    bam_name="${bam_dir}${sample}.unimappers.deduped"
  else
    bam_name="${bam_dir}${sample}.unimappers"
  fi

  if [[ "${strandedness}" == "unstranded" ]]; then
    strand_type=0
  elif [[ "${strandedness}" == "stranded" ]]; then
    strand_type=1
  else
    strand_type=2
  fi

  echo "featureCount for ${sample}"
  if [[ "${sequencing_type}" == "SE" ]]; then
    featureCounts -T "${cores}" -s "${strand_type}" -t gene -g gene_id \
      -a "${gtf}" -o "${featureCounts_dir}${sample}.gene.featureCounts.txt" "${bam_name}.bam"
  fi
  if [[ "${sequencing_type}" == "PE" ]]; then
    featureCounts -T "${cores}" -p --countReadPairs -s "${strand_type}" -t gene -g gene_id \
      -a "${gtf}" -o "${featureCounts_dir}${sample}.gene.featureCounts.txt" "${bam_name}.bam"
  fi
done

########################################
### 10. Calculate spike-in size factors
########################################
module purge
module load R/4.2.1

echo "Calculate size factors using DESeq2"
Rscript "${script_dir}size_factor.R" "${work_dir}" "${sample_sheet}" "${gtf}"

##################################
### 11. Convert bam to bigwig (spike-in normalized)
##################################
module purge
module load samtools/1.15.1
module load bedtools/2.30.0
module load GenomeToolset

size_factors="${work_dir}analysis/size_factors_deseq2.txt"
[ ! -f "${size_factors}" ] && echo "ERROR: Size factor file DOES NOT exist." && exit 1

for sample in "${sample_names[@]}"
do
  bam_name="${bam_dir}${sample}.unspliced"
  scale_factor=$(awk -v s="${sample}" '$1==s {print $2^-1}' "${size_factors}")
  echo "Convert bam to bigwig for ${sample}, scale factor: ${scale_factor}"

  ### Forward strand
  bedtools genomecov -ibam "${bam_name}.fwd.bam" -bg -split -strand - -scale "${scale_factor}" \
    | sort --parallel="${cores}" -k1,1 -k2,2n \
    > "${bigwig_dir}${sample}.spikein.normalized.fwd.bedgraph"
  bedGraphToBigWig "${bigwig_dir}${sample}.spikein.normalized.fwd.bedgraph" \
    "${chrsize}" "${bigwig_dir}${sample}.spikein.normalized.fwd.bw"

  ### Reverse strand
  bedtools genomecov -ibam "${bam_name}.rev.bam" -bg -split -strand + -scale "${scale_factor}" \
    | sort --parallel="${cores}" -k1,1 -k2,2n \
    > "${bigwig_dir}${sample}.spikein.normalized.rev.bedgraph"
  bedGraphToBigWig "${bigwig_dir}${sample}.spikein.normalized.rev.bedgraph" \
    "${chrsize}" "${bigwig_dir}${sample}.spikein.normalized.rev.bw"

  ### Remove intermediate bedgraph files
  rm "${bigwig_dir}${sample}.spikein.normalized.fwd.bedgraph"
  rm "${bigwig_dir}${sample}.spikein.normalized.rev.bedgraph"
done

##############################
### END
##############################

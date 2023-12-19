#!/usr/bin/env bash
#SBATCH --job-name=bash
#SBATCH --mail-type=END,fal
#SBATCH --mail-user=hailiu@sund.ku.dk
#SBATCH --nodes=1                          # Run all processes on a single node	
#SBATCH --ntasks=1                         # Run a single task
#SBATCH --cpus-per-task=10                 # Number of CPU cores per task
#SBATCH --mem=32gb                         # Job memory request
#SBATCH --time=0-20:00:00
#SBATCH --output=drb_ttchem_seq_log%j.log

########################################
### Author: Haiyue Liu
### Date: 15-12-2023
### This is the strandard pipeline for DRB/TTchem-seq2 data processing. 
### Libraries: single-ended reversely-stranded reads with UMI (11-mer) sequences. We demultiplexed the UMI sequences as R2 which are attached to R1 as headers in our pileline. 
### We also use yeast spike-ins for normalization of the nascent RNA libraries. For this, we creasted a huamn-yeast combine genome reference for mapping. 
### This main step for data processing in this pipeline are:
### 1. Attache UMIs to reads header
### 2. Adaptor trimming & fastqc 
### 3. STAR alignment
### 4. Extract uniquely mapping reads
### 5. Deduplication (optional)
### 6. Removal of exon-intron-exon mapping reads
### 7. Split forward and reverse strand reads
### 8. Count mapped reads (library size)
### 9. Quantify gene counts
### 10. Calculte yeast spike-in size factors 
### 11. Convert bam to bigwig format 
#############################################

#################################################
#################################################
### configuration & references
#################################################
#################################################    
### Configuration: this section is project-specific. You need to modify those lines based on your own project setup. 
### Tools executables: all the tools are called via module load in this pipeline. If you use a diffrent system, you could change those module load lines.
cores=10                                                                              ### required; 
work_dir="/maps/projects/dan1/people/mjh723/projects/DRB_TTchem_seq/"                 ### required; project directory; your fastq folder containing all raw reads should have been loaded into this folder 
script_dir="${work_dir}script/"                                                       ### required; script direcoty. The R function "size_factor.R" should be found there
sample_sheet="${work_dir}script/samplesheet.tsv"                                      ### required; one example tsv file can be found in this github repository. The column names and orders matter for the pipeline 
sequencing_type="SE"                                                                  ### required; options: "SE", "PE"
UMI=true                                                                              ### required; options: true, false
UMI_sequence="NNNNNNNNNNN"                                                            ### required if UMI=true; the length of UMIs can vary depends on the project
strandedness="reverse-stranded"                                                       ### required; options: "reverse-stranded", "stranded", "unstranded"
fastq_suffix=".fastq.gz"                                                              ### required; options: ".fastq.gz", ".fq.gz"                                                                 
### index & annotation
reference_dir="/maps/projects/dan1/people/mjh723/reference/gencode/"
# genome_fasta="${reference_dir}/genome/human_yeast/GRCh38_sacCer3.fa"
index_dir="${reference_dir}index/GRCh38.gencode.v43_sacCer3.108_star_2.7.9a_index"
gtf="${reference_dir}annotation/human_yeast/GRCh38.gencode.v43_sacCer3.108.gtf"
chrsize="${index_dir}/chrNameLength.txt"
#################################################
#################################################

############################
### check config parameters
############################
[ ! -d ${work_dir} ] && echo "Working directory DOES NOT exists." && exit
[ ! -f ${sample_sheet} ] && echo "Samplesheet file DOES NOT exists." && exit
[ ! -f ${script_dir}size_factor.R ] && echo "size_factor.R file DOES NOT exists." && exit
sample_ids=$(cat $sample_sheet | awk 'NR>1 {print $2}')
sample_names=$(cat $sample_sheet | awk 'NR>1 {print $3}')

#############################
### create sub-directories
#############################
### As far as the work_dir is given, these are automatically creasted.  
fastq_dir="${work_dir}fastq/"
fastqc_dir="${work_dir}fastqc/"
trimmed_fastq_dir="${work_dir}trimmed_fastq/"
tmp_dir="${work_dir}tmp/"
bam_dir="${work_dir}bam/"
featureCounts_dir="${work_dir}featureCounts/"
bigwig_dir="${work_dir}bigwig/"
analysis_dir="${work_dir}analysis/"
mkdir -p ${fastqc_dir}
mkdir -p ${trimmed_fastq_dir}
mkdir -p ${tmp_dir}
mkdir -p ${bam_dir}
mkdir -p ${featureCounts_dir}
mkdir -p ${bigwig_dir}
mkdir -p ${analysis_dir}
cd ${work_dir}

########## Data processing ###############

#########################################
### 1. Attach UMIs in R2 to header of R1
#########################################
module purge
module load dangpu_libs python/3.7.13 umi_tools/1.1.4
for sample_id in ${sample_ids[@]}
do
  sample=$(cat ${sample_sheet} | awk -v s=$sample_id '{ if($2==s) {print $3} }')
  if [[ "${sequencing_type}" == "SE" ]]; then
    echo "Extract UMIs from R2 and attach them to R1 read headers for $sample \n"
    R1_in=${fastq_dir}${sample_id}_R1_001${fastq_suffix}
    umi_in=${fastq_dir}${sample_id}_R2_001${fastq_suffix}
    R1_out=${fastq_dir}${sample}_umi_attached${fastq_suffix}
    #########################################
    ### extract UMI sequenc from read2 and add it to read headers of read1
    ### umi_tools is not multiple threading
    #########################################
    umi_tools extract --extract-method=string --bc-pattern=${UMI_sequence} --stdin ${umi_read} --read2-in=${R1_in} --read2-out=${R1_out} 
  fi 
  if [[ "${sequencing_type}" == "PE" ]]; then
    echo "Paired-end reads \n"
    echo "Extract UMIs from R2 and attach them to R1 and R3 read headers for $sample \n"
    R1_in=${fastq_dir}${sample_id}_R1_001.fastq.gz
    R2_in=${fastq_dir}${sample_id}_R3_001.fastq.gz
    umi_read=${fastq_dir}${sample_id}_R2_001.fastq.gz
    R1_out=${fastq_dir}${sample}_R1_umi_attached.fastq.gz
    R2_out=${fastq_dir}${sample}_R2_umi_attached.fastq.gz
    #########################################
    umi_tools extract --extract-method=string --bc-pattern=${UMI_sequence} --stdin=${umi_read} --read2-in=${R1_in} --stdout=${R1_out} --read2-stdout
    umi_tools extract --extract-method=string --bc-pattern=${UMI_sequence} --stdin=${umi_read} --read2-in=${R2_in} --stdout=${R2_out} --read2-stdout
  fi
done
############################
### 2. QC & trim adaptor
############################
module purge  
module load anaconda3/2021.11
module load openjdk/13.0.1 perl/5.26.3 fastqc/0.11.9
module load pigz trimgalore/0.6.6
module load dangpu_libs python/3.7.13 cutadapt/4.1
### match the strings between fastq.gz(fq.gz) and sample_id, can be empty string
if $UMI; then name_suffix="umi_attached"; else name_suffix=""; fi
for sample in ${sample_names[@]}
do
  echo "$sample : FastQC -- Adaptor trmming -- FastQC "
	if [[ "${sequencing_type}" == "SE" ]]; then
	  ### fastQC
    fastqc -t ${cores} -o ${fastqc_dir} ${fastq_dir}${sample}_R1_${name_suffix}${fastq_suffix}
	  ### trim adaotors & fastQC on trimmed reads
	  trim_galore --cores ${cores} --basename ${sample} --output_dir ${trimmed_fastq_dir} --fastqc --fastqc_args "-o ${fastqc_dir} -t ${cores}" ${fastq_dir}${sample}_R1_${name_suffix}${fastq_suffix}
    # ### rename files with sample name
    # mv ${fastqc_dir}${sample}_R1_${name_suffix}_fastqc.html ${fastqc_dir}${sample}_fastqc.html
    # mv ${fastqc_dir}${sample}_R1_${name_suffix}_fastqc.zip ${fastqc_dir}${sample}_fastqc.zip
    # mv ${trimmed_fastq_dir}${sample}_R1_${name_suffix}${fastq_suffix}_trimming_report.txt ${trimmed_fastq_dir}${sample}_trimming_report.txt
	fi
	if [[ "${sequencing_type}" == "PE" ]]; then
		### fastQC
    fastqc -t ${cores} -o ${fastqc_dir} ${fastq_dir}${sample}_R1_${name_suffix}${fastq_suffix} ${fastq_dir}${sample}_R2_${name_suffix}${fastq_suffix}
	  ### trim adaotors & fastQC on trimmed reads
	  trim_galore --cores ${cores} --paired --basename ${sample} --output_dir ${trimmed_fastq_dir} --fastqc --fastqc_args "-o ${fastqc_dir} -t ${cores}" ${fastq_dir}${sample}_R1_${name_suffix}${fastq_suffix} ${fastq_dir}${sample}_R2_${name_suffix}${fastq_suffix}
    # ### rename files with sample name
    # mv ${fastqc_dir}${sample}_R1_${name_suffix}_fastqc.html ${fastqc_dir}${sample}_R1_fastqc.html
    # mv ${fastqc_dir}${sample}_R2_${name_suffix}_fastqc.html ${fastqc_dir}${sample}_R2_fastqc.html
    # mv ${fastqc_dir}${sample}_R1_${name_suffix}_fastqc.zip ${fastqc_dir}${sample}_R1_fastqc.zip
    # mv ${fastqc_dir}${sample}_R2_${name_suffix}_fastqc.zip ${fastqc_dir}${sample}_R2_fastqc.zip
    # mv ${trimmed_fastq_dir}${sample}_R1_${name_suffix}${fastq_suffix}_trimming_report.txt ${trimmed_fastq_dir}${sample}_R1_trimming_report.txt
	  # mv ${trimmed_fastq_dir}${sample}_R2_${name_suffix}${fastq_suffix}_trimming_report.txt ${trimmed_fastq_dir}${sample}_R2_trimming_report.txt
  fi
done
###########################
### 3. STAR alignment
###########################
module purge
module load gcc star/2.7.9a
module load samtools
rm -r ${tmp_dir}
mkdir -p ${tmp_dir}
for sample in ${sample_names[@]}
do
  echo "Mapping for $sample"
	if [[ "${sequencing_type}" == "SE" ]]; then
    STAR \
  	  --runThreadN ${cores} \
  	  --genomeDir ${index_dir} \
  	  --readFilesIn ${trimmed_fastq_dir}${sample}_trimmed.fq.gz \
    	--readFilesCommand zcat \
    	--outSAMattributes NH HI AS NM MD XS \
  	  --outSAMstrandField intronMotif \
  	  --outSAMmultNmax 1 \
  	  --outSAMunmapped None \
  	  --outFileNamePrefix ${bam_dir}${sample}. \
    	--outSAMtype BAM SortedByCoordinate \
  	  --outTmpDir ${tmp_dir}${sample}
  fi
  if [[ "${sequencing_type}" == "PE" ]]; then
    STAR \
  	  --runThreadN ${cores} \
  	  --genomeDir ${index_dir} \
  	  --readFilesIn ${trimmed_fastq_dir}${sample}_val_1_trimmed.fq.gz ${trimmed_fastq_dir}${sample}_val_2_trimmed.fq.gz \
    	--readFilesCommand zcat \
    	--outSAMattributes NH HI AS NM MD XS \
  	  --outSAMstrandField intronMotif \
  	  --outSAMmultNmax 1 \
  	  --outSAMunmapped None \
  	  --outFileNamePrefix ${bam_dir}${sample}. \
    	--outSAMtype BAM SortedByCoordinate \
  	  --outTmpDir ${tmp_dir}${sample}
  fi
	### index bam files
	samtools index -@ ${cores} ${bam_dir}${sample}.Aligned.sortedByCoord.out.bam
	### remove immediate files
	rm ${bam_dir}${sample}.Log.out
	rm ${bam_dir}${sample}.Log.progress.out
done
########################################
### 4. Extract uniquely mapping reads
########################################
module purge
module load samtools/1.15.1
for sample in ${sample_names[@]}
do
  echo "Exact unimappers for $sample"
  in_bam=${bam_dir}${sample}.Aligned.sortedByCoord.out.bam
  out_bam=${bam_dir}${sample}.unimappers.bam
  if [[ "${sequencing_type}" == "SE" ]]; then
    samtools view --threads ${cores} -q 255 -F 0x2 ${in_bam} -o ${out_bam}
  fi
  if [[ "${sequencing_type}" == "PE" ]]; then
    samtools view --threads ${cores} -q 255 -f 0x2 ${in_bam} -o ${out_bam}
  fi
  samtools index -@ ${cores} ${out_bam}
done
#######################################
### 5. Deduplication of UMIs (optional)
#######################################
module purge
module load parallel
module load dangpu_libs python/3.7.13 umi_tools/1.1.4
module load samtools/1.15.1
### DRB_20minR_rep1.deduped.bam takes 21.90361h
if $UMI; then
  echo "deduplication"
  for sample in ${sample_names}
  do
    umi_tools dedup --umi-separator="_" -I ${bam_dir}${sample}.unimappers.bam -S ${bam_dir}${sample}.unimappers.deduped.bam --log=${bam_dir}${sample}_dedup.log --output-stats=${bam_dir}${sample}
  done
  ## To improve run speed, we can also run all samples in parallel
  # find ${bam_dir} -name "*.unimappers.bam" | parallel -j ${cores} umi_tools dedup --umi-separator="_" -I {} -S {.}.deduped.bam --log={.}.dedup.log --output-stats={}
  # find ${bam_dir} -name "*.deduped.bam" | parallel -j ${cores} samtools index -@ 1 {}
fi
###############################################
### 6. Removal of exon-intron-exon reads
###############################################
## In our data, we observed the reads that mapped across exon-intron-exons junctions can give high background noise. These reads are a small (<2%) amount of the total reads but they are highly present in control samples in which the DRB were not washed out. Therefore, we suspect those reads are mainly non-specific enrichment. We removed these reads for the downstream analysis.
module purge
conda activate ngs     ### this is to activate python (pysam is needed)
split_reads_script="${script_dir}split_unspliced_spliced_reads.py"
for sample in ${sample_names[@]}
do
  echo "split unspliced and spliced reads for ${sample}"
  ### bam files
  if $UMI; then
    bam_name=${bam_dir}${sample}.unimappers.deduped
  else
    bam_name=${bam_dir}${sample}.unimappers
  fi
  ~/.conda/envs/ngs/bin/python ${split_reads_script} ${bam_name}.bam ${bam_name}.clean.bam ${bam_name}.splice.junction.bam
done
conda deativate
###############################################################
### 7. Split reads transcribed from forward and reverse strands
###############################################################
module purge
module load samtools/1.15.1
for sample in ${sample_names[@]}
do
  bam_name=${bam_dir}${sample}.clean
  samtools index -@ ${cores} ${bam_name}.bam
  echo "Split forward and reverse strand $sample"
  if [[ "${sequencing_type}" == "SE" && "${strandedness}" == "reverse-stranded" ]]; then
    ### Forward strand reads -- R1 mapped to reverse strand (include read reverse strand 16)
    samtools view -b -f 16 --threads ${cores} ${bam_name}.bam -o ${bam_name}.fwd.bam
    samtools index -@ ${cores} ${bam_name}.fwd.bam
    ### Reverse strand reads -- mapped to forward strand (exclude read reverse strand 16)
    samtools view -b -F 16 --threads ${cores} ${bam_name}.bam -o ${bam_name}.rev.bam
    samtools index -@ ${cores} ${bam_name}.rev.bam
  fi
  if [[ "${sequencing_type}" == "PE" && "${strandedness}" == "reverse-stranded" ]]; then
    ### Forward strand reads -- R1 mapped to reverse strand -- R2 mapped to forward strand
    ### PE1 -- include first in pair 64 + read reverse strand 16
    samtools view -b -f 80 --threads ${cores} ${bam_name}.bam -o ${bam_name}.R1.fwd.bam
    ### PE2 -- include second in pair 128 + exclude read reverse strand 16
    samtools view -b -f 128 -F 16 --threads ${cores} ${bam_name}.bam -o ${bam_name}.R2.fwd.bam
    ### combine PE1 & PE2
    samtools merge -l 9 --threads ${cores} -f ${bam_name}.fwd.bam ${bam_name}.R1.fwd.bam ${bam_name}.R2.fwd.bam
    samtools index -@ ${cores} ${bam_name}.fwd.bam

    ### Reverse strand reads -- R1 mapped to forward strand -- R2 mapped to reverse strand
    ### PE1 -- include first in pair 64 + exclude read reverse strand 16
    ### PE2 -- inlcude second in pair 128 + include read read reverse strand 16
    samtools view -b -f 64 -F 16 --threads ${cores} ${bam_name}.bam -o ${bam_name}.R1.rev.bam
    samtools view -b -f 144 --threads ${cores} ${bam_name}.bam -o ${bam_name}.R2.rev.bam
    ### combine PE1 & PE2
    samtools merge -l 9 --threads ${cores} ${bam_name}.rev.bam ${bam_name}.R1.rev.bam ${bam_name}.R2.rev.bam
    samtools index -@ ${cores} ${bam_name}.rev.bam
    ### remove intermediated files
    # rm ${bam_name}.R1.fwd.bam
    # rm ${bam_name}.R2.fwd.bam
    # rm ${bam_name}.R1.rev.bam
    # rm ${bam_name}.R2.rev.bam
  fi
done
#########################################
### 8. Count mapped reads (library size)
#########################################
module purge
module load samtools/1.15.1
reads_number="${work_dir}/analysis/reads_number.txt"
# printf "sample_name\t\total\t\dedup\t\unimapped\t\human\t\yeast\n" > ${reads_number}
if $UMI; then
  printf '%s\t' 'sample_name' 'unimapper' 'deduped' 'clean' 'human' > ${reads_number}
  printf '\n' >> ${reads_number}
else
  printf '%s\t' 'sample_name' 'unimapper' 'clean' 'human' > ${reads_number}
  printf '\n' >> ${reads_number}
fi

for sample in ${sample_names[@]}
do
  echo $sample
  n_unimapper=$(samtools view --threads ${cores} -c ${bam_dir}${sample}.unimappers.bam)
  if $UMI; then
    bam_name=${bam_dir}${sample}.unimappers.deduped
    n_deduped=$(samtools view --threads ${cores} -c ${bam_name}.bam)
  else
    bam_name=${bam_dir}${sample}.unimappers
  fi
  ### count the total reads in the exon-intron-exon reads removed bam files
  n_clean=$(samtools view --threads ${cores} -c ${bam_name}.clean.bam)
  ### only count human reads in the cleaned bam files
  n_human_clean=$(samtools idxstats ${bam_name}.clean.bam | awk '/^chr/ {s+=$3}END{print s}')
  ### save the numbers to file
  if $UMI; then
    echo $sample | awk -v OFS="\t" -v unimapper="$n_unimapper" -v deduped="$n_deduped" -v clean="$n_clean" -v human_clean="$n_human_clean" '{print $0, unimapper, deduped, clean, human_clean}' >> ${reads_number}
  else
    echo $sample | awk -v OFS="\t" -v unimapper="$n_unimapper" -v clean="$n_clean" -v human_clean="$n_human_clean" '{print $0, unimapper, clean, human_clean}' >> ${reads_number}
  fi
done
#########################################
### 9. Quantify gene counts
#########################################
module purge
module load subread/2.0.3
for sample in ${sample_names[@]}
do
  if $UMI; then
    bam_name=${bam_dir}${sample}.unimappers.deduped
  else
    bam_name=${bam_dir}${sample}.unimappers
  fi
  ### strandedness
  if [[ "${strandedness}" == "unstranded" ]]; then
    strand_type=0
  elif [[ "${strandedness}" == "stranded" ]]; then
    strand_type=1
  else
    strand_type=2
  fi
  echo "featureCount for $sample"
  if [[ "${sequencing_type}" == "SE" ]]; then
    featureCounts -T ${cores} -s ${strand_type} -t gene -g gene_id -a ${gtf} -o ${featureCounts_dir}${sample}.gene.featureCounts.txt ${bam_name}.bam
  fi
  if [[ "${sequencing_type}" == "PE" ]]; then
    featureCounts -T ${cores} -p --countReadPairs -s ${strand_type} -t gene -g gene_id -a ${gtf} -o ${featureCounts_dir}${sample}.gene.featureCounts.txt ${bam_name}.bam
  fi
done
########################################
### 10. Calculate spike-in size factors
########################################
module purge
module load R/4.2.1
echo "Calculate size factors using DESeq2"
Rscript ${script_dir}size_factor.R $work_dir $sample_sheet $gtf
##################################
### 11. Convert bam to bigwig
##################################
module purge
module load samtools/1.15.1
module load bedtools/2.30.0
module load GenomeToolset
## Load the sample size factor file calcaulted using DESeq2 and the reads number file 
size_factors=${work_dir}"analysis/size_factors_deseq2.txt"
reads_number="${work_dir}analysis/reads_number.txt"
[ ! -f ${scale_factors} ] && echo "Size factor file DOES NOT exists!"
[ ! -f ${reads_number} ] && echo "Reads number file DOES NOT exists!"
for sample in ${sample_names[@]}
do
  # echo ${sample}
  if $UMI; then
    bam_name=${bam_dir}${sample}.unimappers.deduped
  else
    bam_name=${bam_dir}${sample}.unimappers
  fi
  ###############################
  ### normalized to spike-ins
  ### bam -> bedgraph -> bigwig 
  ###############################
  scale_factor=$( cat ${size_factors} | awk -v s=$sample '{ if($1==s) {print $2^-1} }' )
  echo "Convert bam to bigwig for sample: ${sample} and scale it by ${scale_factor}"
  ### forward strand
  bedtools genomecov -ibam ${bam_name}.fwd.bam -bg -split -strand - -scale ${scale_factor} | sort --parallel=${cores} -k1,1 -k2,2n > ${bigwig_dir}${sample}.spikein.normalized.fwd.bedgraph
  bedGraphToBigWig ${bigwig_dir}${sample}.spikein.normalized.fwd.bedgraph ${chrsize} ${bigwig_dir}${sample}.spikein.normalized.fwd.bw
  ### strand
  bedtools genomecov -ibam ${bam_name}.rev.bam -bg -split -strand + -scale ${scale_factor} | sort --parallel=${cores} -k1,1 -k2,2n > ${bigwig_dir}${sample}.spikein.normalized.rev.bedgraph
  bedGraphToBigWig ${bigwig_dir}${sample}.spikein.normalized.rev.bedgraph ${chrsize} ${bigwig_dir}${sample}.spikein.normalized.rev.bw
  ### rm intermediate files
  rm ${bigwig_dir}${sample}.spikein.normalized.fwd.bedgraph
  rm ${bigwig_dir}${sample}.spikein.normalized.rev.bedgraph
  ###############################
  ### normalized to library size
  ### bam -> bedgraph -> bigwig 
  ################################
  scale_factor=$( cat ${reads_numbers} | awk -v s=$sample '{ if($1==s) {print 10^6/$2} }')  ### uniquely mapped reads
  echo "Convert bam to bigwig for sample: ${sample} and scale it by ${scale_factor}"
  ### forward strand
  bedtools genomecov -ibam ${bam_name}.fwd.bam -bg -scale ${scale_factor} | sort --parallel=${cores} -k1,1 -k2,2n > ${bigwig_dir}${sample}.unimappers.libsize.normalized.fwd.bedgraph
  bedGraphToBigWig ${bigwig_dir}${sample}.unimappers.libsize.normalized.fwd.bedgraph ${chrsize} ${bigwig_dir}${sample}.libsize.normalized.fwd.bw
  ### reverse strand
  bedtools genomecov -ibam ${bam_name}.rev.bam -bg -scale ${scale_factor} | sort --parallel=${cores} -k1,1 -k2,2n > ${bigwig_dir}${sample}.unimappers.libsize.normalized.rev.bedgraph
  bedGraphToBigWig ${bigwig_dir}${sample}.unimappers.libsize.normalized.rev.bedgraph ${chrsize} ${bigwig_dir}${sample}.libsize.normalized.rev.bw
  rm ${bigwig_dir}${sample}.unimappers.libsize.normalized.fwd.bedgraph
  rm ${bigwig_dir}${sample}.unimappers.libsize.normalized.rev.bedgraph
done
##############################
### end
##############################

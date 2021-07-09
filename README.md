# cfSNV

## INTRODUCTION 

cfSNV is an ultra-sensitive and accurate somatic SNV caller designed for cfDNA sequencing. Taking advantage of modern statistical models and machine learning approaches, cfSNV provides hierarchical mutation profiling and multi-layer error suppression, including error suppression in read mates, site-level error filtration and read-level error filtration. cfSNV can be freely used for educational and research purposes by non-profit institutions and U.S. government agencies only under the UCLA Academic Software License. For information on the use for a commercial purpose or by a commercial or for-profit entity, please contact Prof. Xiangong Jasmine Zhou (https://zhoulab.dgsom.ucla.edu/).



## PREREQUISITE PACKAGES

cfSNV was developed on UCLA Hoffman2 cluster. All the analyses in this study were done there. The system information of UCLA Hoffman2 cluster is:
Linux n6015 2.6.32-754.14.2.el6.x86_64 #1 SMP Tue May 14 19:35:42 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux
There is no non-standard hardware required.

1) java 1.8.0_111 (Just need to be compatible to GATK and picard tools)
2) GATK 3.8.0: https://github.com/broadgsa/gatk/releases/tag/3.8
3) FLASh2: https://github.com/dstreett/FLASH2
4) bwa 0.7.17: http://bio-bwa.sourceforge.net/
5) samtools 1.9: https://sourceforge.net/projects/samtools/files/samtools/
6) picard 2.18.4: https://github.com/broadinstitute/picard/releases/tag/2.18.4
7) bedtools 2.26.0: https://github.com/arq5x/bedtools2/releases/tag/v2.26.0
8) python2 version 2.4.2+
9) numpy 1.13.3
10)   pandas 0.20.3
11)   scipy 1.1.0
12)   sklearn 0.19.1
13)   decimal 1.70


## INSTALLATION

cfSNV algorithm is implemented using python scripts for now. There is no specific installation needed, if all required packages and tools listed in PREREQUISITE PACKAGES are successfully installed. All the source codes (python codes) are saved under scripts folder.

The use of python codes are wrapped using bash scripts (".sh" scripts). These scripts are started with PBS (portable batch system, the job scheduling system on UCLA Hoffman2) configurations. The configuration parameters basically set the resources required to run the bash scripts. If the users' computing environment is under different scheduling system, the user could change the configuration part according to their needs.

The paths to the required tools are set for UCLA Hoffman2 cluster. To run cfSNV algorithm in other environment, absolute paths to the required tools in PREREQUISITE PACKAGES section in ".sh" scripts need to be changed to where these packages are installed in users' computing environment.


## USAGE

The cfSNV algorithm is packed using bash scripts. In each bash script, the code starts with configuration settings: LIST OF INPUT FILES and REQUIRED TOOLS.
   LIST OF INPUT FILES is a input table which contains paths to the required input files. Each line in this table is a complete set of required input files for the script. The bash script could be set to run for different lines of inputs in this table in parallel.
   REQUIRED TOOLS is a list of absolute paths to required tools.

To run cfSNV, the users could run the bash scripts in the following order:

1. getbam.align.jobarray.sh
   This bash script basically maps raw reads from fastq (gzip or not) files to the reference genome. The required inputs are shown below:<br>
      path: working directory (also where the bam files are written)<br>
      fastq1: full path to fastq file (gzip or not) for read 1<br>
      fastq2: full path to fastq file (gzip or not) for read 2<br>
      file_id: file prefix for the output<br>
      SNP_database: full path to common SNP databases (in vcf format)<br>
      reference: full path to reference genome fasta file (indexed by samtools and bwa)<br>
   For a plasma-WBC sample pair from a cancer patient, both plasma sample and WBC sample need to be processed with this bash script. The output bam file from WBC sample will be used as matched normal control in variant calling.

2. getbam.align_after_merge.jobarray.sh
   This bash script is to first merge the overlapping read mates in cfDNA sequencing data and then map raw reads from fastq (gzip or not) files to the reference genome. The required inputs are shown below:<br>
      path: working directory (also where the bam files are written)<br>
      fastq1: full path to fastq file (gzip or not) for read 1<br>
      fastq2: full path to fastq file (gzip or not) for read 2<br>
      file_id: file prefix for the output<br>
      SNP_database: full path to common SNP databases (in vcf format)<br>
      reference: full path to reference genome fasta file (indexed by samtools and bwa)<br>
   For a plasma-WBC sample pair from a cancer patient, ONLY plasma sample needs to be processed with this bash script. There are two output bam files: ${file_id}.extendedFrags.recal.bam and ${file_id}.notCombined.recal.bam, which are alignment results for overlapping read mates (after merging) and non-overlapping read mates (cannot be merged) respectively.

3. dataprep0.get_depth.jobarray.sh
   This bash script calculates an approximate depth of coverage in the targeted regions. The resulting depth of coverage will serve as a parameter for the downstream analysis. The required inputs are shown below:<br>
      path: working directory (also where the temporary files and the outputs are written)<br>
      plasma_unmerged_file: full path to plasma bam file from 1. getbam.align.jobarray.sh<br>
      full_target_bed: full path to target regions in bed format (e.g. capture bait in WES)<br>
      sample_id: file prefix for the output<br>

4. dataprep1.pileup_rebase.quick.jobarray.sh
   This bash script pre-processes the plasma and the WBC bam files to pre-select potential mutation candidates. The required inputs are shown below:<br>
      path: working directory (also where the temporary files and the outputs are written)<br>
      plasma_unmerged_file: full path to plasma bam file from 1. getbam.align.jobarray.sh<br>
      normal_file: full path to WBC bam file from 1. getbam.align.jobarray.sh<br>
      plasma_merged_extendedFrags_file: full path to plasma extendedFrags bam file from 2. getbam.align_after_merge.jobarray.sh<br>
      plasma_merged_notCombined_file: full path to plasma notCombined bam file from 2. getbam.align_after_merge.jobarray.sh<br>
      target_bed: full path to target regions in bed format<br>
      reference: full path to reference genome fasta file (indexed by samtools and bwa)<br>
      sample_id: file prefix for the output<br>
   For large datasets, e.g. WES, this step is going to take a quite long time. In our analyses, we manually parallelized this step by splitting the full target bait into chromosomes, i.e. a big bed file to some smaller bed files (n=24). Then for a single pair of plasma-WES samples, this step could be parallelized into 24 smaller sub-tasks, and thus shorten the running time. Basically, we preprocessed samples at the chromosome level using different nodes (one node per chromosome). The outputs from sub-tasks would be automatically merged together in 5. iterate.estimate_call_filter.quick.jobarray.sh. Somatic variant calling only focuses on the regions included in the targeted bed.

5. iterate.estimate_call_filter.quick.jobarray.sh
   This bash script implemented the iterative and intellegent search of mutation clusters in cfDNA. The required inputs are shown below:<br>
      path: working directory (also where the temporary files and the outputs are written)<br>
      sample_id: file prefix for the output<br>
      SNP_database: full path to common SNP databases (in vcf format)<br>

6. classify.machine_learning_classification_reads.jobarray.sh
   This bash script implemented the read-level filtration of sequencing errors and some site-level filtration. The required inputs are shown below:<br>
      path: working directory (also where the temporary files and the outputs are written)<br>
      sample_id: file prefix for the output<br>
      unmerged_plasma_bam: full path to plasma bam file from 1. getbam.align.jobarray.sh<br>
      reference: full path to reference genome fasta file (indexed by samtools and bwa)<br>
      reference_dict: full path to reference genome dict file (generated by picard CreateSequenceDictionary)<br>
   This script outputs the final results into the working directory, with file name started with ${sample_id}.

Please note that 1 & 2 are basically pre-processing of the raw reads. The results of these two files are saved under the working directory specified. 3-6 are the main statistical and machine learning part of the algorithm. All the outputs from 3-6 are intermediate results and put into a temporary folder under the working directory. The final output files are under the working directory, but all the intermediate files will be removed at the end (could be disabled by comment out the last command in 6. classify.machine_learning_classification_reads.jobarray.sh). Some steps in this workflow could be combined together to further save computing resources, but we haven't done it for now for easier debugging and specfic requirements at UCLA Hoffman2.

## OUTPUT

There are two output files from cfSNV. These files could be found under the working directory (path).

1) estimated tumor fraction from Jenks natural breaks optimization<br>
   <plasma sample name>.jenks_estimate<br>
   This file contains three numbers. The first one is the estimated tumor fraction. The second one is the number of somatic mutations used in the estimation. The last one is the number of mutation groups involved in this estimation.

2) variant report table<br>
   <plasma sample name>.variant_report.txt<br>
   This table is a VCF-like report. It contains seven columns:<br>
         a. CHROM: chromosome name<br>
         b. POSITION: chromosome coordinate<br>
         c. ID: variant ID (usually as ".")<br>
         d. REF: wildtype base (according to reference)<br>
         e. VAR: variant base<br>
         f. SCORE: log-likelihood ratio between somatic variant and other circumstances<br>
         g. VAF: variant allele frequency




## DEMO: EXAMPLE DATA


A demo is provided along with the scripts. Example data and required reference files are saved under demo folder. Due to the upload file limit of github, some files are compressed. Please decompress them before testing: chr22.fa.bwt.7z, chr22.fa.zip, and chr22_dbSNP.zip.
1) chr22.fa: reference genome (from UCSC genome browser hg19)
2) chr22.fa.*: genome indices created for bwa, samtools, GATK and picard tools
3) example_output folder: expected output files
   a) plasma.jenks_estimate: estimated tumor fraction
   b) plasma.variant_report.txt: variant report table
4) example_target_regions.bed: targeted regions (e.g. capture regions of WES)
5) normal_raw.fq1.gz: normal control (white blood cell) read 1
6) normal_raw.fq2.gz: normal control (white blood cell) read 2
7) plasma_raw.fq1.gz: plasma cfDNA read 1
8) plasma_raw.fq2.gz: plasma cfDNA read 2
9) input_list_getbam_align.txt: LIST OF INPUT FILES for 1. getbam.align.jobarray.sh
   This file contains two lines for plasma sample and normal (WBC) sample respectively.
10) input_list_getbam_align_after_merge.txt: LIST OF INPUT FILES for 2. getbam.align_after_merge.jobarray.sh
   This file contains one line for plasma sample ONLY.
11) input_list_dataprep0.txt: LIST OF INPUT FILES for 3. dataprep0.get_depth.jobarray.sh
   This file contains one line for the plasma-WBC pair.
12) input_list_dataprep1.txt: LIST OF INPUT FILES for 4. dataprep1.pileup_rebase.quick.jobarray.sh
   This file contains one line for the plasma-WBC pair.
13) input_list_iterate.txt: LIST OF INPUT FILES for 5. iterate.estimate_call_filter.quick.jobarray.sh
   This file contains one line for the plasma-WBC pair.
14) input_list_classify.txt: LIST OF INPUT FILES for 6. classify.machine_learning_classification_reads.jobarray.sh
   This file contains one line for the plasma-WBC pair.

 This demo was expected to run less than 20min. Please change the paths to LIST OF INPUT FILES and PREREQUISITE TOOLS in each .sh script. Then the algorithm could be run simply through the bash scripts.

**open the code package and run .sh scripts in order**
1. getbam.align.jobarray.sh
2. getbam.align_after_merge.jobarray.sh
3. dataprep0.get_depth.jobarray.sh
4. dataprep1.pileup_rebase.quick.jobarray.sh
5. iterate.estimate_call_filter.quick.jobarray.sh
6. classify.machine_learning_classification_reads.jobarray.sh



## Citation

Li, Shuo, Zorawar Noor, Weihua Zeng, Xiaohui Ni, Zuyang Yuan, Frank Alber, Wenyuan Li, Edward B. Garon, and Xianghong Zhou. "Sensitive detection of tumor mutations from blood and its application to immunotherapy prognosis." medRxiv (2020): 2019-12.


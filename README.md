# NGS_HCV

Manual version 1.0  

*Sebastiano Montante - April 2018*


# Summary

1. Installation  
2. Algorithm flow  
3. Instructions  
    3.1 Running first script  
    3.2 Mandatory arguments  
    3.3 Optional arguments  
    3.4 Examples
4. Output files


# 1. Installation

This software is developed for 64bit Linux system and was tested on Ubuntu.

The following packages must be installed on your system to execute the pipeline.

**Python modules:**
- Biopython (Bio package)
- Numpy
- Pandas

**Bash packages:**
- Bwa
- Mothur
- Ncbi-blast+

The software automatically checks the presence of all required packages. A warning message is shown on the terminal indicating the missing package.


# 2. Algorithm flow

Here is a briefly description of the flow, along with the pipeline.
The first step of the pipeline includes the check of the reads folder and the extraction of input
files if compressed.
During the second step, the reads that will be processed are selected according to user's indications. At this point there are two alternatives:

**1. If the genotype is unknown, a blast against a pool of filtered HCV sequences is performed
to find the genotype of the input reads.**

**2. If the genotype is known the blast is not executed and the pipeline goes directly to the next
step.**

The next step is mapping. Input reads are mapped against a specific reference sequence, related to their gene and genotype, located in the reads folder. If the user provides a personal reference sequence, the mapping is performed against that sequence. The output of this step is a SAM file.

The last step is the SAM processing, the most time-consuming stage of the pipeline. SAM processing briefly consists in the extraction of the CIGAR information for each reads mapped correctly. Based on the CIGAR information, pair-wise alignment between each read and the reference sequence is reconstructed, also taking in consideration the presence of insertions and
deletions. Then each aligned pair is translated into amino-acids sequence. Reads sequence is analyzed and compared to a rule set where are reported the resistance-associated mutations and their relevant positions. The prevalence of the amino-acids, including insertions and deletions, based on the reference position is reported on two tables. The prevalence is calculated upon the total number of reads mapped. Finally, if a resistance-associated mutation reported in the rule is found whithin reads analyzed, their prevalence, drugs influenced by these mutations and other relevant informations are reported on
another table.


# 3. Instructions

**3.1 Execution of the first script**

The first script is *NGS_HCV_analysis.py* and is located inside the software folder. It must be executed with the system's Python interpreter in order to start running the pipeline. This script is responsable for the execution of all other scripts. The user must indicate absolute or relative path to this script.

Example:

`python path/software_folder/NGS_HCV_analysis.py`

**3.2 Mandatory arguments**

These arguments are necessary for the correct execution of the pipeline but someone changes based on the type of input. The -R argument is mandatory regardless the input. Next to the -R argument the user must specify the correct directory where reads to be analyzed are located.

`- R directory_path_input_reads`

Other mandatory arguments depend on the type of input:

*1. First case:*

The input is a single file of reads derived from single-end NGS or from pair-end NGS where the forward and reverse reads are in the same file. In this case, in addition to the -R argument, user must specify only the -r argument indicating input name located inside the reads_folder.

`-r name_single_reads_file_to_be_processed`

*2. Second case:*

The input consists of two reads that have to be merged. The resulting merged reads are processed. In this case, in addition to -R argument, user must specify -rf and -rrv arguments, indicating respectively the name of the file containing forward and reverse reads sequences.

`-rf name_file_reads_forward`

`-rrv name_file_reads_reversed`

*3. Third situation:*

The input consists of two reads derived from a paired-end NGS in which both forward and reverse reads have to be processed, without merging. In this case, in addition to the -R argument user must specifiy -rf and -rrv arguments, indicating respectively the name of the file containing forward and reverse reads sequences. Also, user must specify the -p argument to avoid merging.

- `-rf name_file_reads_forward-rrv name_file_reads_reversed`

- `-p`

Note: input files must be in .fastq format. The last mandatory argument is -gn. Next to this argument user must specify the gene of the amplicon.

- `gn {NS3,NS5B,NS5A}`

Valid input genes are indicated within the braces.

**3.3. Optional arguments**

- `wd new_working_directory`

With this argument user can change working directory in which the software processes temporary files and where output files are reported. Default working directory is the folder where the terminal has been opened (user has to remember that relative reads' folder path changes once the working directory changes)

- `btch n` 

User can analyze a batch of “n” reads casually extracted from the input files. Useful to perform fast testing runs.

- `tv n`

User can specify a threshold_value and can consider only relevant mutations with a prevalence above “n”.

- `prc n`

This argument manages number of threads to execute BLAST, BWA mapping and mothur merging.

- `gt {"1a","1b","2","3","4","unknown"}`

If known, user can indicate reads' genotype. Default is unknown. If unknown, a blast is performed. BLAST aligns all the reads against a pool of thousands HCV sequences. This pool covers main HCV genotypes. The software considers genotype with the highest prevalence as the genotype of the reads. Mapping of the reads will always be executed against a specific reference based on the gene and the genotype of the reads. Gene's and genotype's informations are also used during the rule set comparison, because each relevant resistance mutation and associated drug changes based on gene and genotype. 

- `usrf path_user_ref_folder/user_ref.fasta`

User can also provides a personal reference, specifing the correct path of the reference file next to -usrf argument. Arguments (both optional and mandatory) can be inserted in any order.

**3.4 Examples**

*1) Processing a SAM file derived by a single-end mapping. No genotype specified.*

In this case a blast is executed to find genotype of the input reads.

`python path/software_folder/NGS_HCV_analysis.py -R path/reads_folder -r reads_1.fastq -gn
NS5A`

*2) Processing a SAM file derived by a single-end mapping* 

In this case user already knows the reads' genotype, so blast is not executed. Only 100 reads will be analyzed.

`python path/software_folder/NGS_HCV_analysis.py -R path/reads_folder -r reads_1.fastq -gn
NS5A -gt 1a -btch 100`

*3) Processing a SAM file derived by a paired-end mapping.*

User has indicated the name of forward and reversed reads and -p argument.

`python path/software_folder/NGS_HCV_analysis.py -R path/reads_folder -p -rf reads_1.fastq -rrv
reads_2.fastq -gn NS5A`

*4) Processing a SAM file derived by a single-end mapping.*

The user has indicated a personal reference. In this case the reads will be mapped against this reference. Ony 100 reads will be analyzed. No blast executed because the genotype has been indicated.

`python path/software_folder/NGS_HCV_analysis.py -R path/reads_folder -r reads_1.fastq -gn
NS5A -gt 1a -btch 100 -usrf user_ref_path/user_ref.fasta`

*5) Processing a SAM file derived by a single-end mapping.*

Input consists of two reads that are merged into a single file because the -p argument has not been specified. Blast is executed.

`python path/software_folder/NGS_HCV_analysis.py -R ./reads_folder -rf reads_1.fastq -rrv
reads_2.fastq -gn NS5A4`


# 4. Output files

There are 7 output files:

- **1) df_all_ref_pos_all_aa_prevalence.csv:**

An excel file reporting the amino-acids prevalences compared to the total number of mapped reads based on the reference positions. Table also reports coverage of each position and the prevalence of deletions.

- **2) df_all_ref_pos_aa_ins_prevalence.csv**

An excel file reporting only the prevalence of insertions based on reference positions. Table also reports percentage of reads that have insertions in a position compared to the total number of mapped reads and their absolute frequency.

- **3) df_rule_set**

Portion of the rule set reporting the relevant resistance-associated mutations founded among the reads (based on gene and genotype), related drugs, effect (possible resistant or resistant) and the reference from scientific literature.

- **4) df_rules_found_vs_prevalences**

This file is a table reporting relevant mutations founded among their associated prevalences. If -tv argument was specified, mutations whithin the chosen threshold are not considered.

- **5) drugs_resistant_and_susceptible**

A text file reporting list of drugs that virus can be resistant or susceptible to.

- **6) The raw SAM file before processing.**

- **7) log_file:**

A text file reporting details about operations performed by the pipeline. It reports total number of reads in the SAM file, number of reads correctly mapped and analyzed. It reports all the possible errors that can occur while running the pipeline. It also reports various warnings that can informs user about biological aspects emerged during the processing. For example, it informs user about a low mapping percentage when the number of mapped reads is below 50%. If no relevant mutations are founded, output only includes a warning message advising the user and output files are 1,2,6,7. The other files remain empty.

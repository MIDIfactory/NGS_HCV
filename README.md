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

The next step is the mapping phase. The input reads are mapped against a specific reference
sequence,based on their gene and genotype, located in the reads folder. If the user provides a
personal reference sequence,the mapping of the reads is performed against this user reference.
The output of this step is a SAM file.
The last step is the SAM processing,the most time-consuming phase of the pipeline.
The SAM processing briefly consists in the extraction of the CIGAR information for each reads
mapped correctly. Based on the CIGAR information, the pair-wise alignment between each read and
the reference sequence is reconstructed, taking also in consideration the presence of insertions and
deletions respect to the reference. Then each aligned pair is translated into their amino-acids
sequence. The reads sequence is analyzed and compared to a rule set where are reported the
resistance-associated mutations and their relevant positions.
The prevalence of the amino-acids, including insertions and deletions, based on the reference
position is reported in two tables. The prevalence is calculated on the total number of reads
mapped.
Finally, If a resistance-associated mutation reported in the rule is found in the reads analyzed,their
prevalence, the drugs influenced by these mutations and other relevant informations,are reported in
another table.


# 3. Running instructions

**3.1 Execution of the start script**

The start script must be executed with the Python interpreter of the system in order to start the
execution of the pipeline.
The name of the start script is NGS_HCV_analysis.py located inside the software folder.
The start script is responsible for the execution of all the pipeline’s scripts.
The user must indicate the correct absolute or relative path to the start script located inside the
software folder.

Example:

`python path/software_folder/NGS_HCV_analysis.py`

**3.2 mandatory arguments**

These arguments are necessary to the correct execution of the pipeline and some mandatory
arguments changes based on the type of input.
The -R argument is mandatory regardless the type of input.
Next to the - R argument the user must specify the correct directory where the file of the reads to be
analyzed is located.

`- R directory_path_input_reads`

The other mandatory arguments depend on the type of input to be analyzed:

*1. first situation:*

The input is a single file of reads derived from single end next generation sequencing
experiment.
The input is a single file of reads derived from pair end next generation sequencing
experiment where the forward and reverse reads sequences are in the same file.
In this case, in addition to the -R argument, the user must specify only the -r argument
indicating the name of the input reads file,located inside the reads_folder.

`-r name_single_reads_file_to_be_processed`

*2. Second situation:*

◦ The input consists of two reads files that must be merged and the resulting merged reads
are processed
in this case in addition to -R argument the user must specify the -rf and -rrv
arguments,indicating respectively the name of file that contains the forward reads sequences
and the name of the file that contains the reverse reads sequences.

`-rf nome_file_reads_forward`

`-rrv nome_file_reads_reversed`

*3. Third situation:*

The input consists of two reads files derived from a pair end next generation sequencing
experiment in which both the forward and reverse reads must be processed,without
merging
In this case in additon to the -R argument the user must specify the -rf and -rrv
arguments,indicating respectively the name of file that contains the forward reads sequences
and the name of the file that contains the reverse reads sequences. In addition the user must
specify the -p argument to avoid the execution of the merging

- `-rf name_file_reads_forward-rrv name_file_reads_reversed`

- `-p`

Note: the input reads files must be in .fastq format.
The last mandatory argument is the -gn argument. Next to this argument the user must specify the
gene of the amplicon.

- `gn {NS3,NS5B,NS5A}`

The valid input genes are indicated within the braces.3.3 Optional arguments
These arguments are not necessary to the correct running of the pipeline, but they can be useful.

- `wd new_working_directory`

With this argument the user can change the working directory in which the software processes the
temporary files and where output files are reported. The default working directory is the folder
where the terminal has been opened (the user must remember that the relative reads folder path
changes once the working directory changes)

- `btch n`

The user can choose to analyze a batch of “n” reads casually extracted by the input files. Useful to
perform rapid testing runs.

- `tv n`

The user can specify a threshold_value. The user can consider only relevant mutations with a
prevalence above “n”.

- `prc n`

The user can select the number of processors. This argument controls the number of processors to
execute the blast, the bwa mapping, and the mothur merging.

- `gt {"1a","1b","2","3","4","unknown"}`

The user can indicate the genotype of the reads to be analyed if known. The default is unknown.
If unknown , a blast is performed. The blast align all the reads against a pool of thousands HCV
sequences.
This pool covers the main HCV genotypes.
The software considers the genotype with the highest prevalence as the genotype of the reads.
The mapping of the reads will be always executed against a specific reference based on the gene
and genotype of the reads. The gene and genotype informations are also used during the rule set
comparison phase.
Because each relevant resistance mutation and associated drug changes based on gene and
genotype.

- `usrf path_user_ref_folder/user_ref.fasta`

The user can also provide a personal reference, specifing the correct path of the reference file next
to the -usrf argument.The arguments (both optional and mandatory ) can be inserted in any order.

**3.4 Examples of correct runs**

*1) Processing of the SAM file derived by a single end mapping. No genotype specified.*

In this case a blast is executed to find the genotype of the input reads.

`python path/software_folder/NGS_HCV_analysis.py -R path/reads_folder -r reads_1.fastq -gn
NS5A`

*2) Processing of the SAM file derived by a single end mapping* 

In this case the user already knows the reads genotype so the blast is not executed Only 100 reads will be analyzed.

`python path/software_folder/NGS_HCV_analysis.py -R path/reads_folder -r reads_1.fastq -gn
NS5A -gt 1a -btch 100`

*3) Processingt of the SAM file derived by a pair end mapping.*

The user has indicated the name of the reads forward and reversed and the -p argument.

`python path/software_folder/NGS_HCV_analysis.py -R path/reads_folder -p -rf reads_1.fastq -rrv
reads_2.fastq -gn NS5A`

*4) Processing of the SAM file derived by a single end mapping.*

The user has indicated a personal
reference. In this case the reads will be mapped against this reference. Ony 100 reads will be
analyzed. No blast executed because the genotype has been indicated.

`python path/software_folder/NGS_HCV_analysis.py -R path/reads_folder -r reads_1.fastq -gn
NS5A -gt 1a -btch 100 -usrf user_ref_path/user_ref.fasta`

*5) Processing of the SAM file derived by a single end mapping.*

The input are two reads files that
are merged into a single merged reads file because the -p argument has not been specified. Blast is
executed.

`python path/software_folder/NGS_HCV_analysis.py -R ./reads_folder -rf reads_1.fastq -rrv
reads_2.fastq -gn NS5A4`


# 4. Output files

There is 7 output files:

- **df_all_ref_pos_all_aa_prevalence.csv:**

An excel file that reports the amino-acids prevalences
respect to the total number of reads mapped based on reference position. The table reports
also the coverage of each position and the prevalence of deletions.

- **df_all_ref_pos_aa_ins_prevalence.csv**

an excel file that reports only the prevalence of
insertions based on reference positions.
The table reports also the the total percentage of reads that have insertions in a position
respect to the total number of reads mapped and their absolute frequency.

- **df_rule_set**

portion of the rule set that reports all the relevant resistance-associated
mutations founded in the reads (based on gene and genotype), the related drugs, the effect
(possible resistant or resistant) and the reference from scientific literature.

- **df_rules_found_vs_prevalences**

This file is a table that reports the relevant mutations
founded with their associated prevalences. If -tv argument has been specified, the mutations
under the chosen threshold are not considered.

- **drugs_resistant_and_susceptible**

a text file that reports the list of drugs for which the virus
can be resistant or susceptible.

- **The raw SAM file before processing.**

- **log_file:**

a text file that reports details about the operations performed by the pipeline.
It reports the total number of reads in the SAM file, the number of reads correctly mapped
and analyzed.
It reports all the possible errors that can occur during the running of the pipeline. It reports
also the different warnings that can inform the user about biological aspects emerged during
the processing.
For example it informs the user about a low mapping percentage , when the number of reads
mapped is below the 50%.
If no relevant mutations is founded, the output includes only a warning message that advises the
user and the output files 1,2,6,7. The other files remains empty.

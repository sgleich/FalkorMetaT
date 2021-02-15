**Eukaryotic MetaT Bioinformatic Pipeline**  
**By: Samantha Gleich**  
**Most of the code here is taken or modified from Sarah K. Hu (https://github.com/shu251/SPOT_metatranscriptome)**  
**Last modified: 2/15/21**

## Trim Sequences - Trimmomatic
```
java -jar /path/to/directory/Trimmomatic-0.32/trimmomatic-0.32.jar PE [sample-laneNum-readDir.fastq] ILLUMINACLIP:/path/to/directory/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:10 TRAILING:10 SLIDINGWINDOW:25:10 MINLEN:50
```
## Remove ERCC spike-in mix
Know the directory that your ERCC fasta sequences are in (i.e. /path/to/directory/ERCC92b.fa). Here, laneNum represents the lane number (each sample was sequenced on two lanes) and R1/R2 represents the read direction (PE sequences)  

First, align ERCC reads to fastq files.
```
/path/to/directory/trinity-2.1.1/util/align_and_estimate_abundance.pl --seqType fq --left [sample-laneNum-R1-trimmed.fastq] --right [sample-laneNum-R2-trimmed.fastq] --transcripts /path/to/directory/ERCC92b.fa --est_method RSEM --aln_method bowtie --trinity_mode --output_dir sample-laneNum
```
Using the custom perl script written by Sarah Hu, remove ERCC reads.
```
samtools view bowtie.bam | perl /path/to/directory/remove_ERCC_reads_PE.pl /path/to/directory/trimmed/[sample-laneNum-R1-trimmed.fastq] /path/to/directory/trimmed/[sample-laneNum-R2-trimmed.fastq] /path/to/directory/trimmed/[sample-laneNum-R1-trimmed2.fastq] /path/to/directory/trimmed/[sample-laneNum-R2-trimmed2.fastq]
```
## Concatenate reads from different lanes
If reads were sequenced on more than one lane, concatenate sequences from same sample (sample) that have the same read direction (R1/R2).
```
cat [sample-L2-R1-trimmed2.fastq] [sample-L3-R1-trimmed2.fastq] > [sample-R1-trimmed2.fastq]
```
## Separate rRNA from mRNA - SortMeRNA
First, we will merge the PE reads.
```
/path/to/directory/sortmerna-2.1b/scripts/merge-paired-reads.sh /path/to/directory/trimmed/[sample-R1-trimmed2.fastq] /path/to/directory/trimmed/[sample-R2-trimmed2.fastq] [sample-merged.fastq]
```
Then, we will separate sequences that match the rRNA database from those that do not.
```
sortmerna --ref /path/to/directory/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,/path/to/directory/sortmerna-2.1b/index/silva-euk-18s-db --reads [sample-merged.fastq] --sam --fastx --aligned [sample-rrna] --other [sample-norrna] --log -v --paired_in --best 1
```
Then, we can can unmerge PE reads that didn't hit the rRNA database
```
/path/to/directory/sortmerna-2.1b/scripts/unmerge-paired-reads.sh [sample-norrna.fastq] [sample-R1-unmerged.fastq] [sample-R2-unmerged.fastq]
```
Lastly, we can umerge PE reads that did hit the rRNA database. These can be used to look at rRNA-based taxonomy later on (miTag analysis)
```
/path/to/directory/sortmerna-2.1b/scripts/unmerge-paired-reads.sh [sample-rrna.fastq] [sample-R1-unmerged-rrna.fastq] [sample-R2-unmerged-rrna.fastq]
```
## Assemble reads - MEGAHIT
Concatenate all sequences that will be used in one assembly together. Here, I will concatenate all replicates from the same eddy/depth. Leave R1 and R2 reads separate from one another.
```
cat [sample1-R1-unmerged.fastq] [sample2-R1-unmerged.fastq] [sample3-R1-unmerged.fastq] > [sample-R1.fastq]
```
Once sequences that will be used in each coassembly have been concatenated into single R1 and R2 fastq files, we can run an individual coassembly: 
```
megahit -t 16 -1 [sample-R1.fastq] -2 [sample-R2.fastq] -o [sample_mhit_out]
```
## Estimate transcript abundances - salmon
Fix quality strings
```
Add code
```
Make index using assembly output
```
salmon quant -i all_depths_eddies_index -l A -1 /path/to/directory/[sample-R1-unmerged.fastq] -2 /path/to/directory/[sample-R2-unmerged.fastq] -o [sample.quant]
```

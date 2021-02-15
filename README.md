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
Add assembly tag to the beginning of each contig name. This will be important for downstream purposes. Then we will concatenate all of our assemblies together so we have one big list of contigs. In the example below, we will pretend we have 3 assemblies: assembly_ID1, assembly_ID2, and assembly_ID3.
```
sed 's/^>/>assembly_ID1_/' [final.contigs.fa] > [assembly_ID1_final.contigs2.fa]
cat [assembly_ID1_final.contigs2.fa] [assembly_ID2_final.contigs2.fa] [assembly_ID3_final.contigs2.fa] > [all_contigs.fa]
```
## Estimate transcript abundances - salmon
Salmon will be used to align reads to the assembled contigs. This will allow us to estimate the abundance of a specific transcript. 

If you have any sequences with truncated quality strings, these sequences must be removed before running salmon otherwise, you'll get an error. Here, sequences with truncated quality strings will be removed using a custom python script (remove_bad_quality_strings.py). 
```
python ./remove_bad_quality_strings.py
```
Make index using concatenated assembly output (all_contigs.fa). Make sure contig names have sample info added to them (as done with the sed command above)
```
salmon index -t [all_contigs.fa] -i salmon_index --type quasi -k 31
```
Now execute salmon. You will have one salmon output file per sample. 
```
salmon quant -i all_depths_eddies_index -l A -1 /path/to/directory/[sample-R1-unmerged.fastq] -2 /path/to/directory/[sample-R2-unmerged.fastq] -o [sample.quant]
```
## Taxonomic classification - blastx
We will now assign taxonomy to the contigs we obtained from each of our assemblies. 
```
diamond blastx -d /path/to/custom/db/EukZoo.dmnd -q [assemblyID1_final.contigs2.fa] –sensitive -e 0.01 -o [assemblyID1_blastx_out]
```
Now we will BLAH
```
MORE CODE
```
## Predict proteins - GeneMarkS
Now we can predict proteins from the contigs we've obtained from each of our assemblies. 
```
GeneMarkS-T/gmst.pl --fnn -faa [assemblyID1_final.contigs2.fa]
```
GeneMarkS outputs can be used to for functional annotation via GhostKoala for KEGG annotation (https://www.kegg.jp/ghostkoala/) and/or eggNOG-mapper 

## Functional annotation - eggNOG
We will take the predicted protein .faa files obtained from GeneMarkS, to run eggNOG mapper for functional annotation. 
```
emapper.py -i /path/to/genemark/output/[assemblyID1_final.contigs2.faa] --output assemblyID1_eggnog -m diamond
```
To compile eggNOG, GhostKoala, blastx, and salmon data and to do transcript count normalization, see Falkor_MetaT_DataProcessing_Normalization.R

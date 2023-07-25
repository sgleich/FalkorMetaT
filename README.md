**Falkor MetaT Bioinformatic Pipeline**  
**By: Samantha Gleich**  
**Most of the code here is taken from or modified from Dr. Sarah K. Hu (https://github.com/shu251/SPOT_metatranscriptome)**  
**Last modified: 7/23/23**
\
![](static/Contour.pdf)\

## Trim Sequences - Trimmomatic
Trim sequences using trimmomatic version 0.38
```
java -jar /path/to/directory/Trimmomatic-0.38/trimmomatic-0.38.jar PE [sample-laneNum-readDir.fastq] ILLUMINACLIP:/path/to/directory/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:10 TRAILING:10 SLIDINGWINDOW:25:10 MINLEN:50
```
## Remove ERCC spike-in mix
Make trinity reference for removing ERCC reads
```
/path/to/directory/trinityrnaseq-2.8.4/util/align_and_estimate_abundance.pl --transcripts /path/to/directory/ERCC92b.fa --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir out
```
Know the directory that your ERCC fasta sequences are in (i.e. /path/to/directory/ERCC92b.fa). Here, laneNum represents the lane number (each sample was sequenced on two lanes) and R1/R2 represents the read direction (PE sequences)  

First, align ERCC reads to fastq files.
```
/path/to/directory/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl --seqType fq --left [sample-laneNum-R1-trimmed.fastq] --right [sample-laneNum-R2-trimmed.fastq] --transcripts /path/to/directory/ERCC92b.fa --est_method RSEM --aln_method bowtie --trinity_mode --output_dir sample-laneNum
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
## Separate rRNA from mRNA - SortMeRNA (version 2.0)
First, we will merge the PE reads.
```
/path/to/directory/sortmerna-2.1b/scripts/merge-paired-reads.sh /path/to/directory/trimmed/[sample-R1-trimmed2.fastq] /path/to/directory/trimmed/[sample-R2-trimmed2.fastq] [sample-merged.fastq]
```
Make index for PR2 database through sortmerna (PR2 version 12)
```
indexdb_rna --ref pr2_version_4.12.0_18S_mothur.fasta,pr2_version_4.12.0_18S_mothur.idx
```

Then, we will separate sequences that match the rRNA database from those that do not.
```
sortmerna --ref /path/to/directory/pr2_version_4.12.0_18S_mothur.fasta,/path/to/directory/pr2_version_4.12.0_18S_mothur.idx --reads [sample-merged.fastq] --sam --fastx --aligned [sample-rrna] --other [sample-other] --log -v --paired_in --best 1
```
Then, we can can unmerge PE reads that didn't hit the rRNA database
```
/path/to/directory/sortmerna-2.1b/scripts/unmerge-paired-reads.sh [sample-norrna.fastq] [sample-R1-unmerged.fastq] [sample-R2-unmerged.fastq]
```
Lastly, we can umerge PE reads that did hit the rRNA database. These can be used to look at rRNA-based taxonomy later on (miTag analysis)
```
/path/to/directory/sortmerna-2.1b/scripts/unmerge-paired-reads.sh [sample-rrna.fastq] [sample-R1-unmerged-rrna.fastq] [sample-R2-unmerged-rrna.fastq]
```
## miTag analysis
First we will split the libraries (.fna file; use qiime1 for this analysis)
```
split_libraries_fastq.py -i /path/to/directory/[sample-R1-unmerged-rrna.fastq] --barcode_type 'not-barcoded' --sample_ids sample-R1-rrna -q 15 -n 5 -o sample-rrna
```
Now we will assign taxonomy to our rRNA sequences using qiime1
```
assign_taxonomy.py -i seqs.fna -o sample-rrna-tax -r /path/to/directory/pr2_version_4.12.0_18S_mothur.fastq -t /path/to/directory/pr2_version_4.12.0_18S_mothur.idx --similarity 0.97
```

## Assemble reads - MEGAHIT (version 1.2.8)
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
## Estimate transcript abundances - salmon (version 1.5.1)
Make index using concatenated assembly output (all_contigs.fa). Make sure contig names have sample info added to them (as done with the sed command above)
```
salmon index -t [all_contigs.fa] -i salmon_index
```
Now execute salmon. You will have one salmon output file per sample. 
```
salmon quant -i all_depths_eddies_index -l A -1 /path/to/directory/[sample-R1-unmerged.fastq] -2 /path/to/directory/[sample-R2-unmerged.fastq] -o [sample.quant]
```
## Predict proteins - GeneMarkS
Now we can predict proteins from the contigs we've obtained from each of our assemblies. 
```
GeneMarkS-T/gmst.pl --fnn -faa [assemblyID1_final.contigs2.fa]
```
GeneMarkS outputs can be used to for functional annotation via GhostKoala for KEGG annotation (https://www.kegg.jp/ghostkoala/) and/or eggNOG-mapper 

## Taxonomic classification - EUKulele (version 1.4.0)
We will now assign taxonomy to the contigs we obtained from each of our assemblies. In the example below, the directory "directory" has the .fnn file obtained from GeneMarkS.
```
EUKulele --sample_dir /path/to/directory -m mets --n_ext fnn
```
## Functional annotation - eggNOG (version 5.0)
We will take the predicted protein .faa files obtained from GeneMarkS to run eggNOG mapper for functional annotation. 
```
emapper.py -i /path/to/genemark/output/[assemblyID1_final.contigs2.faa] --output assemblyID1_eggnog -m diamond
```
## Cluster proteins - mmseqs2
We will cluster the predicted proteins. First we will create a database using mmseqs2.
```
mmseqs easy-cluster all_seqs.fnn clusterRes tmp --min-seq-id 0.8 -c 0.95 --cov-mode 1 
```
To compile eggNOG, EUKulele, salmon, and mmseqs2 output data into long and wide format dataframes, see Compile_Data.R. The compiled dataframe can be used to normalize the dataframe using Normalize_Data.R.

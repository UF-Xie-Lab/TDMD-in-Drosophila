# TDMD in Drosophila pipeline
## 1 TDMD identification analysis
### 1.1 Cutadapt （https://cutadapt.readthedocs.io/en/stable/） in FASTQ
cutadapt -a TGGAATTCTCGGGTGCCAAG -A GATCGTCGGACTGTAGAACT -o `test_R1_cut.fastq` -p `test_R2_cut.fastq` `test_R1.fastq test_R2.fastq` --minimum-length 18 -j 10

### 1.2 Pear (https://cme.h-its.org/exelixis/web/software/pear/doc.html) in FASTQ
pear -f `test_R1_cut.fastq` -r `test_R2_cut.fastq` -j 10 -o `test`

### 1.3 Collapse (http://hannonlab.cshl.edu/fastx_toolkit/) PCR duplicated reads
fastx_collapser  -i `test.assembled.fastq` -o `test_collapsed.fasta`

### 1.4 Remove UMI sequences from 5' and 3' of FASTA
cutadapt -u 4 -u -4 -m 18 `test_collapsed.fasta` -o `test_cutN.fasta` -j 10              

### 1.5 Run hyb(https://github.com/gkudla/hyb) software
module load hyb

module load unafold

hyb analyse in=`test_cutN.fasta` db=`20220221_dm6_unique` type=mim pref=mim format=fasta

### 1.6 Potential TDMD miRNA-target RNA hybrids identification
python3 CLASH.py Viennad_to_Table -i `test_cutN_comp_20220221_dm6_unique_hybrids_ua` -c `dm6_all_conservation_score.txt` -t `20220221_dm6_unique.fasta` -n `martquery_1109155317_9_name.txt`

```ruby
20220221_dm6_unique.fasta file size is 95MB. The dataset folder has a small size file, called 20220221_dm6_unique_small.fasta.
  
dm6_all_conservation_score.txt file size is 1.4G. The dataset folder has a small size file, called dm6_all_conservation_score_small.txt.
```

  

## 2 miRNA analysis
### 2.1 Deduplicate clean-miRNA-seq reads from BGI
  python3 CLASH.py deduplicate_BGI -i `test_miRNA_BGI.fq`
### 2.2 miRNA abundance calculation
python3 CLASH.py miRNA_abundance -i `test_miRNA_BGI_deduplicated.fa` -d `20220221_dm6_miRNA_database.fasta`

### 2.3 Differential expression level analysis
run `Deseq.R` code

### 2.4 miRNA length distribution (isoform) count
python3 CLASH.py miRNA_length_distribution -i `test_miRNA_BGI_deduplicated.fa` -d `20220221_dm6_miRNA_database.fasta` -m `[--miRNA_sequence]`

e.g. miR-999 length distribution count

*python3 CLASH.py miRNA_length_distribution -i `test_miRNA_BGI_deduplicated.fa` -d `20220221_dm6_miRNA_database.fasta` -m TGTTAACTGTAAGACTGTGTCT*

## 3 poly-A RNA-seq analysis
### 3.1 Buind Drosophila genome database index by Hisat2
hisat2-build `GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna drosphila`

### 3.2 Cutadapt in FASTQ
cutadapt -a AGATCGGAAGAGCACACGTCTGAACT -A AGATCGGAAGAGCGTCGTGTAGGGA -o `test_R1_cut.fq` -p `test_R2_cut.fq` -m 20 -j 12 `test_R1.fastq `test_R2.fastq` 

### 3.3 Mapping
 hisat2 -x drosphila -1 `test_R1_cut.fq` -2 `test_R2_cut.fq` -S test.sam -p 12
 
### 3.4 transfer SAM file to BAM file
samtools sort -@ 12 -O BAM -o `test.bam test.sam`

### 3.5 Count gene abundance
htseq-count `test.bam` `GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff` -f bam -s no -m union -i gene --additional-attr=Parent --additional-attr=Dbxref --additional-attr=gbkey --additional-attr=transcript_id > `test.count`

```ruby
GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff and GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna

can be downloaded from NCBI genome database.
```

### 3.6 Differential expression level analysis
run `Deseq.R` code

### 3.7 Cumulative Fraction Curve analysis

python3 CLASH.py Cumulative_fraction_curve_targetScan -i [--DEseq_file] -a [--all_targets] -c [--conserved_targets] -b [--baseMean] 

or

python3 CLASH.py Cumulative_fraction_curve_targetScan_CLASH -i [--DEseq_file] -a [--all_targets] -c [--conserved_targets] -l [--clash_targets] -t [--clash_targetScan_interacted] -b [--baseMean]

e.g. AGO1 knockout cumulative fraction curve analysis

*python3 CLASH.py Cumulative_fraction_curve_targetScan_CLASH -i `output_DEseq.csv` -a `TargetScan7.2__miR-999-3p.all_predicted_targets.txt` -c `TargetScan7.2__miR-999-3p.conserved_predicted_targets.txt` -l `Exp170_miR_999_targets_unique.txt` -t `Exp265_CLASH_interact_targetScan_All_miR_999.txt` -b 100*

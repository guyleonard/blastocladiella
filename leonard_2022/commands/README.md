# Assembly
## fastqc
```bash
fastqc -t 56 demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.filt.fastq.gz
```

## HiFiAdapterFilt
```bash
pbadapterfilt.sh -t 56 demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads
```

## LJA
```bash
lja -o lja_ass_filt --reads demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.filt.fastq.gz -t 56
```

## purge_dups
```bash
pd_config.py hifiasm_filt.p_ctg.fa blem_pacbio_filt.txt -s illumina_gaiix.txt -l purged
minimap2 -t 112 -xmap-hifi assembly.fasta demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.filt.fastq.gz | gzip -c - > assembly_mapped_filt.paf.g
pbcstat assembly_mapped_filt.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log
split_fa assembly.fasta >assembly.split.fasta
minimap2 -t 112 -xasm5 assembly.split.fasta assembly.split.fasta | gzip -c - > assembly_mapped_self.split.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov assembly_mapped_self.split.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed assembly.fasta
```

## pilon
```bash
bwa aln -t 112 purged.fa s_1_1_trimmed_paired.fastq.gz > purged_vs_s_1_1.sai
bwa aln -t 112 purged.fa s_1_2_trimmed_paired.fastq.gz > purged_vs_s_1_1.sai
bwa sampe purged.fa purged_vs_s_1_1.sai purged_vs_s_1_2.sai s_1_1_trimmed_paired.fastq.gz s_1_2_trimmed_paired.fastq.gz > purged_aligned_gaiix.sam
samtools view -@ 112 -Sb purged_aligned_gaiix.sam > purged_aligned_gaiix.bam
samtools sort -@ 112 purged_aligned_gaiix.bam -o purged_aligned_gaiix_sorted.bam
samtools index -@ 112 purged_aligned_gaiix_sorted.bam

pilon --genome purged.fa --frags purged_aligned_gaiix_sorted.bam --output pilon1 --changes --vcf --tracks --outdir pilon1 > pilon1.log
```

# Annotation
## funannotate
```bash
funannotate clean -i genome.fasta -o genome_cleaned.fasta --exhaustive 2>&1 | tee 1_clean.out
funannotate sort -i genome_cleaned.fasta -o genome_cleaned_sorted.fasta -b c_be_22 2>&1 | tee 2_sort.out
funannotate mask -i genome_cleaned_sorted.fasta -o genome_cleaned_sorted_mtantan.fasta --cpus 112 | tee 3_mask.out
funannotate predict -i genome_cleaned_sorted_mtantan.fasta -o predictions -s "Blastocladiella emersonii" --strain "ATCC 22665" --busco_db fungi --transcript_evidence blasto_ests.fasta -d /databases/funannotate/1.8.9 --trnascan ../trnascan/trnascan.out  --SeqCenter "Oxford Genomics Centre" --cpus 112 --name "H9P43" 2>&1 | tee 4_predict.out
funannotate iprscan -i predictions -c 112 -m local
funannotate remote -i predictions -m antismash -e guy.leonard@gmail.com
funannotate annotate -i predictions --cpus 112 --sbt template.sbt --antismash antismash/Blastocladiella_emersonii_ATCC_22665.scaffolds.gbk --busco_db fungi -d /databases/funannotate/1.8.9/ -s "Blastocladiella emersonii" --strain "ATCC 22665" --force
```

# Reports
```bash
ASSEMBLY=Blastocladiella_emersonii_ATCC_22665.scaffolds.fa
ESTS=blasto_ests.fa
```

## qualimap
```bash
gffread Blastocladiella_emersonii_ATCC_22665.gff3 -T -o Blastocladiella_emersonii_ATCC_22665.gtf
qualimap bamqc -bam blem_22_final_aligned_gaiix_sorted.bam -gff ./Blastocladiella_emersonii_ATCC_22665.gtf -outdir ./illumina -outformat PDF:HTML
qualimap bamqc -bam blem_22_final_aligned_pacbio_sorted.bam  -gff ./Blastocladiella_emersonii_ATCC_22665.gtf -outdir pacbio -outformat PDF:HTML
```

## blat - ests vs assembly
```bash
blat $ASSEMBLY $ESTS -out=blast8 ests_vs_final_assembly.out
cat ests_vs_final_assembly.out | cut -f 1 | sort | uniq | wc -l
```

## busco
### scaffolds
```bash
busco -i Blastocladiella_emersonii_ATCC_22665.scaffolds.fa -o blem22_euk_meta -l /databases/busco/v5/lineages/eukaryota_odb10 -m geno -c 112
busco -i Blastocladiella_emersonii_ATCC_22665.scaffolds.fa -o blem22_fun_meta -l /databases/busco/v5/lineages/fungi_odb10 -m geno -c 112
busco -i Blastocladiella_emersonii_ATCC_22665.scaffolds.fa -o blem22_bac_meta -l /databases/busco/v5/lineages/bacteria_odb10 -m geno -c 112
```

### proteins
```bash
busco -i Blastocladiella_emersonii_ATCC_22665.proteins.fa -o blem22_prot_euk_meta -l /databases/busco/v5/lineages/eukaryota_odb10 -m prot -c 112
busco -i Blastocladiella_emersonii_ATCC_22665.proteins.fa -o blem22_prot_fun_meta -l /databases/busco/v5/lineages/fungi_odb10 -m prot -c 112
busco -i Blastocladiella_emersonii_ATCC_22665.proteins.fa -o blem22_prot_bac_meta -l /databases/busco/v5/lineages/bacteria_odb10 -m prot -c 112
```

## blobtools
### bwa - illumina gaiix
```bash
bwa index $ASSEMBLY
bwa aln -t 112 $ASSEMBLY s_1_1_trimmed_paired.fastq.gz > blem_22_final_vs_s_1_1.sai
bwa aln -t 112 $ASSEMBLY s_1_2_trimmed_paired.fastq.gz > blem_22_final_vs_s_1_2.sai
bwa sampe $ASSEMBLY blem_22_final_vs_s_1_1.sai blem_22_final_vs_s_1_2.sai s_1_1_trimmed_paired.fastq.gz s_1_2_trimmed_paired.fastq.gz > blem_22_final_aligned_gaiix.sam
samtools view -@ 112 -Sb blem_22_final_aligned_gaiix.sam > blem_22_final_aligned_gaiix.bam
samtools sort -@ 112 blem_22_final_aligned_gaiix.bam -o blem_22_final_aligned_gaiix_sorted.bam
samtools index -@ 112 blem_22_final_aligned_gaiix_sorted.bam
```

### minimap - pacbio hifi
```bash
minimap2 -ax map-hifi -t 112 $ASSEMBLY ../../raw_reads/demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.filt.fastq.gz > blem_22_final_vs_pacbio.sam
samtools view -@ 112 -Sb blem_22_final_vs_pacbio.sam > blem_22_final_aligned_pacbio.bam
samtools sort -@ 112 blem_22_final_aligned_pacbio.bam > blem_22_final_aligned_pacbio_sorted.bam 
samtools index -@ 112 blem_22_final_aligned_pacbio_sorted.bam
```

### searches
```bash
blastn -task blastn -db /databases/ncbi/nt/2020-10-15/nt -query Blastocladiella_emersonii_ATCC_22665.scaffolds.fa -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 5 -max_hsps 1 -evalue 1e-10 -num_threads 56 -out assembly.ncbi.blastn.1e-10.out
diamond blastx --query Blastocladiella_emersonii_ATCC_22665.scaffolds.fa --db /databases/diamond/2.0.15/uniprot/reference_proteomes.dmnd --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 56 > assembly.diamond.uniprot.blastx.out
```

### blobtoolkit
```bash
blobtools create --threads 112 --fasta $ASSEMBLY --meta Blastocladiella_emersonii_ATCC_22665.yaml --taxid 1313167 --taxdump /databases/ncbi/new_taxdump/2022-06-06 blem_22
blobtools add --hits assembly.ncbi.blastn.1e-10.out --hits assembly.diamond.uniprot.blastx.out --taxrule bestsum --taxdump /databases/ncbi/new_taxdump/2022-06-06 --threads 112 blem_22
blobtools add --cov blem_22_final_aligned_pacbio_sorted.bam=pacbio --cov blem_22_final_aligned_gaiix_sorted.bam=illumina --threads 112 blem_22
blobtools add --busco ../busco/blem22_euk_meta/run_eukaryota_odb10/full_table.tsv --busco ../busco/blem22_fun_meta/run_fungi_odb10/full_table.tsv --threads 112 blem_22
```

## telomeres
### quick blast check
```bash
blastn -task blastn-short -db assembly.fasta -query telomeric_repeats.fasta -outfmt '6 std slen' -out telomere_scaffolds.blastn

awk '$9>=1 && $9<=1000{print $0}' telomere_scaffolds.blastn | cut -f 2 | sort | uniq >identified_telomere_starts.list
awk '$9<=$13 && $9>=($13-1000){print $0}' telomere_scaffolds.blastn | cut -f 2 | sort | uniq >identified_telomere_ends.list

comm -12 identified_telomere_starts.list identified_telomere_ends.list >identified_telomeres_complete.list
comm -13 identified_telomere_starts.list identified_telomere_ends.list >identified_telomere_start_only.list
comm -23 identified_telomere_starts.list identified_telomere_ends.list >identified_telomere_end_only.list

faSomeRecords assembly.fasta identified_telomeres_complete.list identified_telomeres_complete.fasta
faSomeRecords assembly.fasta identified_telomere_start_only.list identified_telomere_start_only.fasta
faSomeRecords assembly.fasta identified_telomere_end_only.list identified_telomere_end_only.fasta
```

### tapestry
```bash
weave -a assembly.fasta -r pacbio_reads_filt.fq.gz -t AAACCT TTTGGA -o full -c 112
```

### tidk
```bash
tidk search --fasta Catan2_AssemblyScaffolds_Repeatmasked.fasta --string AAACCT --output Catan2_search_aaacct --dir .
tidk search --fasta Allma1_AssemblyScaffolds_Repeatmasked.fasta --string AACCT --output Allma1_search_aacct --dir .

tidk explore --fasta Allma1_AssemblyScaffolds_Repeatmasked.fasta --output Allma1_explore --dir . --length 6

tidk search --fasta purged.fa --string AAACCT --output purged_search --dir .
```

## tRNAScan
```bash
trnascan -o trnascan_blem_22.out -m trnascan_blem_22.stats $ASSEMBLY -c /packages/tRNAscanSE/2.0.9/tRNAscan-SE.conf
```

## kmers
### kat
```
kat hist -t 112 -o blem22.hist demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.filt.fastq.gz
kat comp -o assembly -t 112 -n Blastocladiella_emersonii_ATCC_22665.scaffolds.fa demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.filt.fastq.gz
kat cold -o cold -t 112 Blastocladiella_emersonii_ATCC_22665.scaffolds.fa demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.filt.fastq.gz
kat sect -o sect -t 112 Blastocladiella_emersonii_ATCC_22665.scaffolds.fa demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.filt.fastq.gz
```

### genomescope
```bash
kmc -k21 -t112 -m500 -ci1 -cs10000 demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.filt.fastq.gz kmc_pacbio tmp
kmc_tools transform kmc_pacbio histogram kmc_pacbio.histo -cx10000
genomescope.R -i kmc_pacbio.histo -o kmc_pacbio/ -k 21
```

### smudgeplot
```bash
L=$(smudgeplot.py cutoff kmc_pacbio.hist L)
U=$(smudgeplot.py cutoff kmc_pacbio.histo U)
kmc_tools transform kmc_pacbio -ci"$L" -cx"$U" reduce kmc_pacbio_L"$L"_U"$U"
kmc_tools transform kmc_pacbio -ci"$L" -cx"$U" dump -s kmc_pacbio_L"$L"_U"$U".dump
smudgeplot.py hetkmers -o kmc_pacbio_L"$L"_U"$U" < kmc_pacbio_L"$L"_U"$U".dump
smudgeplot.py plot kmc_pacbio_L"$L"_U"$U"_coverages.tsv
```

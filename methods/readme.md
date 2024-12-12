### 1.Data collection

```
for i in $(cat sample_id); do
prefetch ${i} -X 60G
done
```



### 2.Quality control

```
#quality control
for i in $(cat metagenomic_id); do
trim_galore -o 01.cleandata --gzip --paired ${i}1.raw.fastq.gz ${i}2.raw.fastq.gz
done

#quality check
for i in $(cat data/pair_end_sample_id ); do
fastqc -t 12 -o QC -f fastq 01.cleandata/${i}_1_val_1.fq.gz
fastqc -t 24 -o 02.QC -f fastq 01.cleandata/${i}_2_val_2.fq.gz
done
```



### 3.Assembly and bining

```
#assembly
for i in $(cat metagenomic_id); do
metaspades.py --tmp-dir tmp/${i}.tmp -t 30 -o Assembly-metaspades/${i}_assembly -1 01.cleandata/${i}.R1.raw_val_1.fq.gz -2 01.cleandata/${i}.R2.raw_val_2.fq.gz
done

#bining
#concoct
for i in $(cat concoct_id1); do
bowtie2-build -f 04.seqtk/${i}_seqtk1500/${i}.fa 05.bining/concoct/bam/${i} --threads 32
bowtie2 -p 30  -x 05.bining/concoct/bam/${i} -1 01.cleandata/${i}.R1.raw_val_1.fq.gz -2 01.cleandata/${i}.R2.raw_val_2.fq.gz -S 05.bining/concoct/sam/${i}.sam
samtools faidx 04.seqtk/${i}_seqtk1500/${i}.fa
samtools view -bt 04.seqtk/${i}_seqtk1500/${i}.fa.fai -S 05.bining/concoct/sam/${i}.sam -o 05.bining/concoct/bam_final/${i}.bam
samtools sort -@ 32 -l 9 -O BAM 05.bining/concoct/bam_final/${i}.bam -o 05.bining/concoct/bam_final/${i}.sorted.bam
samtools index 05.bining/concoct/bam_final/${i}.sorted.bam
cut_up_fasta.py 04.seqtk/${i}_seqtk1500/${i}.fa -c 10000 -o 0 --merge_last -b 05.bining/concoct/bed/${i}_10K.bed > 05.bining/concoct/10K/${i}_10K.fa
concoct_coverage_table.py 05.bining/concoct/bed/${i}_10K.bed 05.bining/concoct/bam_final/${i}.sorted.bam > 05.bining/concoct/tsv/${i}.tsv
concoct --composition_file 05.bining/concoct/10K/${i}_10K.fa --coverage_file 05.bining/concoct/tsv/${i}.tsv -b 05.bining/concoct/concoct_output/${i}.concoct_output --threads 64
merge_cutup_clustering.py 05.bining/concoct/concoct_output/${i}.concoct_output_clustering_gt1000.csv > 05.bining/concoct/concoct_output/${i}.clustering_merged.csv
mkdir bining/${i}bins
extract_fasta_bins.py --output_path 05.bining/concoct/${i}bins 04.seqtk/${i}_seqtk1500/${i}.fa 05.bining/concoct/concoct_output/${i}.clustering_merged.csv
done

#metabat2 
for i in $(cat metabat2_id); do
jgi_summarize_bam_contig_depths --outputDepth 05.bining/metabat2/${i}/${i}_final.depth.txt 05.bining/concoct/bam_final/${i}.sorted.bam
metabat2 -m 1500 -t 0 -i 04.seqtk/${i}_seqtk1500/${i}.fa -a 05.bining/metabat2/${i}/${i}_final.depth.txt -o 05.bining/metabat2/${i}/${i}bins/${i} -v
done

#maxbin2
for i in $(cat metabat2_id); do
genomeCoverageBed -ibam 05.bining/concoct/bam_final/${i}.sorted.bam > 05.bining/maxbin2/${i}.histogram.tab
python ~/calculate-contig-coverage.py 05.bining/maxbin2/${i}.histogram.tab
run_MaxBin.pl -contig 04.seqtk/${i}_seqtk1500/${i}.fa -abund 05.bining/maxbin2/${i}.histogram.tab.coverage.tab -max_iteration 50 -out 05.bining/maxbin2/${i}/bin/${i} -thread 32
done
```



### 4.DAS_TOOL

```
for i in $(cat sample_id); do
Fasta_to_Contig2Bin.sh -i 05.bining/maxbin2/${i}/bin/ -e fasta > 06.dastool/tsv/${i}_maxbin2.Contig2Bin.tsv
Fasta_to_Contig2Bin.sh -i 05.bining/metabat2/${i}/${i}bins/ -e fa > 06.dastool/tsv/${i}_metabat2.Contig2Bin.tsv
Fasta_to_Contig2Bin.sh -i 06.dastool/concoctbins/${i}_concoct/ -e fa > 06.dastool/tsv/${i}_concoct.Contig2Bin.tsv
wait
DAS_Tool -i 06.dastool/tsv/${i}_metabat2.Contig2Bin.tsv,06.dastool/tsv/${i}_maxbin2.Contig2Bin.tsv,06.dastool/tsv/${i}_concoct.Contig2Bin.tsv -l metabat,maxbin,concoct -c 04.seqtk/${i}_seqtk1500/${i}.fa -o 06.dastool_results/${i} --search_engine diamond --write_bins --score_threshold 0 -t 32
wait
done
```



### 5.CheckM and Checkm2

```
#checkm
for i in $(cat sample_id); do
checkm lineage_wf -x fa 06.dastool_results/${i}_DASTool_bins/ 07.checkm/${i} -t 32 --tmpdir bin_checkm.tmp
done

#checkm2
checkm2 predict -i bins -o result --allmodels -x .fa --tmpdir tmp -t 64
```

### 6.dRep

```
## 99% ANI (strain level)
for i in $(cat sample_id); do
dRep dereplicate 10.dRep/drep_${i}/ -g 09.genome/${i}/*.fa -p 2 --ignoreGenomeQuality -pa 0.95 -sa 0.99 --S_algorithm fastANI
done

```

### 7.GTDB_tk

```

for i in $(cat sample_id); do
cp 10.dRep/dereplicated_genomes/*.fa 12.GTDB_tk/00.genome/
done
##Taxon classification
gtdbtk classify_wf --genome_dir 12.GTDB_tk/00.genome/ --out_dir 12.GTDB_tk/01.genome_gtdbtk --cpus 30 --scratch_dir gtdbtk.temp --pplacer_cpus 32 --extension fa
```

### 8. Phylogenetic analysis

```
##Phylogenetic analysis
bmge -i 11.GTDB_tk/01.genome_gtdbtk/align/gtdbtk.bac120.user_msa.fasta -t AA -g 0.5 -h 1 -b 1 -w 1 -of 12.tree/out_bmge_trimmed_2.fasta
FastTree out_bmge_trimmed_2.fasta > nxMAGs.tree
```



### 9.BGCs analysis

```
#antismash
concat () {
antismash 11.GTDB_tk/rename/ACT/${1}.fa --taxon bacteria --output-dir 13.antismash/${1} --genefinding-tool prodigal --cb-knownclusters -c 4 --cc-mibig --cb-general --cb-subclusters --fullhmmer
}
export -f concat
cat binid|parallel -j 4 concat {}

#bigscape
/public/home/wuq8022600160/anaconda3/envs/Bigscape/bin/python3 /public/home/wuq8022600160/biosoft/BiG-SCAPE-master/bigscape.py -i allBGCs -o Bigscape -c 8 --cutoffs 0.5 --include_singletons --mode auto --verbose --mix --pfam_dir /public/home/wuq8022600160/biosoft/BiG-SCAPE-master/hmm > run.log
```



### 10.MGE analysis

```
#Extracting the basic information of each BGC from the output file of the BGC
python touch_new.py #Extracting the information of BGC into an excel table
python excel.py #Remove thousand-bit separators from the start and stop sequences of BGC in the excel table
python excel2BED.py #The excel file was converted into BED file containing the contig, start position and end position of BGC
samtools faidx M1.fa
cut -f 1,2 M1.fa.fai > M1.len #Preparing the length of each contig in the MAG
bedtools flank -i M120.bed -g M120.len -l 5000 -r 0 > M120up.bed
bedtools flank -i M1.bed -g genome.len -l 0 -r 5000 -s > M1down.bed
bedtools getfasta -fi M120.fa -bed M120up.bed -fo M120up.fa -name 
bedtools getfasta -s -fi M1.fa -bed down.bed -fo M1down.fa -name 
#rename
sed 's/>Region />/' allup.fa > cleaned_allup.fa
#prodigal
prodigal -i cleaned_allup.fa -a cleaned_allup.faa -d cleaned_allup.fasta -p meta -f gff > cleaned_allup.gff
#mobileOG
diamond makedb --in mobileOG-db_beatrix-1.6.All.faa -d mobileOG.dmnd #构建数据库
chmod +x /mnt/hpc/home/wuq8022600160/biosoft/mobileOG-db-main/mobileOG-pl/mobileOGs-pl-kyanite.sh
PATH=$PATH:/mnt/hpc/home/wuq8022600160/biosoft/mobileOG-db-main/mobileOG-pl/
diamond blastp -q ./data/allgenecalled.faa --db meta.mini.db/mobileOG-db/mobileOG-db-beatrix-1.X.dmnd --outfmt 6 stitle qtitle pident bitscore slen evalue qlen sstart send qstart qend -k 1 -o ./data/mobileOG/mobileOG.tsv -e 1e-5 --query-cover 60 --id 70
python mobileOGs-pl-kyanite.py --o mobileOG/out --i mobileOG/out/all.tsv -m mobileOG-db-beatrix-1.6-All.csv
python /mnt/hpc/home/wuq8022600160/biosoft/mobileOG-db-main/mobileOG-pl/mobileOGs-pl-kyanite.py --o mobileOG/MAGs --i mobileOG/MAGs.tsv -m /mnt/hpc/home/wuq8022600160/database/mobileOG/mobileOG-db-beatrix-1.6-All.csv
```


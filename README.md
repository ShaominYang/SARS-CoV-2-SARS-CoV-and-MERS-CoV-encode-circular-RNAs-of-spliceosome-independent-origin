
# Core for "SARS-CoV-2, SARS-CoV and MERS-CoV encode a novel type of circular RNA"

Circular RNAs (circRNAs) are a newly recognized component of the transcriptome with critical roles in autoimmune diseases and viral pathogenesis. To address importance of circRNA in RNA viral transcriptome, we systematically identified and characterized circRNAs encoded by the RNA genomes of betacoronaviruses using bioinformatical and experimental approaches. We predicted 351, 224 and 2,764 circRNAs derived from SARS-CoV-2, SARS-CoV and MERS-CoV, respectively, and experimentally identified 75 SARS-CoV-2 circRNAs. 

# Requirements
- 1.Bash (Ubuntu, version 18.04)
- 2.Perl [https://www.perl.org](https://www.perl.org/)
- 3.Java [https://javadl.oracle.com](https://javadl.oracle.com/)
- 3.BWA [ttp://bio-bwa.sourceforge.net](http://bio-bwa.sourceforge.net/)
- 4.CIRI2 [https://sourceforge.net/projects/ciri/files/CIRI2/](https://sourceforge.net/projects/ciri/files/CIRI2/) 
- 5.SAMtools [http://www.htslib.org/](http://www.htslib.org/)
- 6.sratoolkit [https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software/)
- 7.aspera [https://www.ibm.com/products/aspera](https://www.ibm.com/products/aspera/)
- 8.Qualimap [http://qualimap.conesalab.org/](http://qualimap.conesalab.org/)
# Preparing Data 
##  Download RNA_seq data (SARS-CoV-2-infected Vero E6 cells at 24 hours post-infection ,PRJNA644588) from ebi
```Shell
for i in ftp.sra.ebi.ac.uk/vol1/srr/SRR121/098/SRR12164498 ftp.sra.ebi.ac.uk/vol1/srr/SRR121/099/SRR12164499 ftp.sra.ebi.ac.uk/vol1/srr/SRR121/000/SRR12164500
do
echo $i
ascp -QT -l 300m -P33001 \
-i ~/miniconda2/envs/rna/etc/asperaweb_id_dsa.openssh \
era-fasp@fasp.sra.ebi.ac.uk:/${i} .
done
```
## Multi-threading to speed up the extraction of fastq from SRA-accessions

```Shell
fasterq-dump --split-3 SRR* -e 16 -p
```
## Pooling all samples together
```Shell
cat SRR12164498_1.fastq SRR12164499_1.fastq SRR12164500_1.fastq > SARS_CoV_2_Vero_E6_24h_1.fastq
cat SRR12164498_2.fastq SRR12164499_2.fastq SRR12164500_2.fastq > SARS_CoV_2_Vero_E6_24h_2.fastq
```
## Remove all SRR files

```Shell
rm SRR*
```

# Generating genome contained African green monkey and SARS-CoV-2

```Shell
cat ChlSab1.1.101.fa NC_045512.2.fasta > ChlSab1.1.101_NC_045512.2.fa
cat ChlSab1.1.101.gtf NC_045512.2.gtf > ChlSab1.1.101_NC_045512.2.gtf
```

# Build  BWA index

```Shell
bwa index ChlSab1.1.101_NC_045512.2.fa
```
![image](https://github.com/ShaominYang/Core-for-SARS-CoV-2-SARS-CoV-and-MERS-CoV-encode-a-novel-type-of-circular-RNA/blob/main/Illustration.jpg)
### <p align="center"> Illustration of CIRI2-based identification of circRNAs </p>
# Running CIRI2 and circ-full pipeline
  
```Shell
for i in SARS_CoV_2_Vero_E6_24h
do
echo $i
mkdir ${i}_ChlSab1.1.101_SARS_CoV_2_output
bwa mem -t 52 ChlSab1.1.101_NC_045512.2.fa ${i}_1.fastq ${i}_2.fastq >${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}.sam
perl ~/CIRI2/CIRI_v2.0.6/CIRI2.pl -I ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}.sam -O ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}.ciri -F ChlSab1.1.101_NC_045512.2.fa -A ChlSab1.1.101_NC_045512.2.gtf -T 24
## Reconstructed SARS-CoV-2 circRNAs circ-full
perl ~/CIRI2/CIRI_AS/CIRI_AS_v1.2.pl -S ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}.sam -C ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}.ciri -F ChlSab1.1.101_NC_045512.2.fa -A ChlSab1.1.101_NC_045512.2.gtf -O ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i} -D yes
java -jar ~/CIRI2/CIRI-full_v2.0/CIRI-full.jar RO1 -1 ${i}_1.fastq -2 ${i}_2.fastq -o ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}
bwa mem -t 52 ChlSab1.1.101_NC_045512.2.fa ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}_ro1.fq > ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}_ro1.sam
java -jar ~/CIRI2/CIRI-full_v2.0/CIRI-full.jar RO2 -r ChlSab1.1.101_NC_045512.2.fa -s ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}_ro1.sam -l 300 -o ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}RO2
java -jar ~/CIRI2/CIRI-full_v2.0/CIRI-full.jar Merge -c ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}.ciri -as ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}_jav.list -ro ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}RO2_ro2_info.list -a ChlSab1.1.101_NC_045512.2.gtf -r ChlSab1.1.101_NC_045512.2.fa -o ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}
unset DISPLAY
java -jar ~/CIRI2/CIRI_vis/CIRI-vis_v1.4.jar -i ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}_merge_circRNA_detail.anno -l ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}_library_length.list -r ChlSab1.1.101_NC_045512.2.fa -min 1
done
  
```
# Running mapping statistics pipeline

```Shell
for i in SARS_CoV_2_Vero_E6_24h
do
echo $i
samtools view -bS ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}.sam > ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}.bam
rm ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}.sam
samtools sort ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}.bam -o ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}_sorted.bam -@ 42
qualimap bamqc -bam ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}_sorted.bam -oc count.matrix -outdir ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}_bamqc -outformat PDF:HTML --java-mem-size=50G
rm ${i}_ChlSab1.1.101_SARS_CoV_2_output/${i}.bam
done
```



# Citations


>1.  Yang S, Zhou H, Cruz-Cosme R, Liu M, Xu J, Niu X, Li Y, Xiao L, Wang Q, Zhu H, Tang Q. Circular RNA profiling reveals abundant and diverse circRNAs of SARS-CoV-2, SARS-CoV and MERS-CoV origin. bioRxiv [Preprint]. 2020 Dec 8:2020.12.07.415422. doi: 10.1101/2020.12.07.415422. PMID: 33330860; PMCID: PMC7743059.
>2.  Yuan Gao†, Jinfeng Wang† and Fangqing Zhao*. CIRI: an efficient and unbiased algorithm for de novo circular RNA identification. Genome Biology (2015) 16:4.
>3.  Yuan Gao, Jinyang Zhang and Fangqing Zhao*. Circular RNA identification based on multiple seed matching. Briefings in Bioinformatics (2017) DOI: 10.1093/bib/bbx014.



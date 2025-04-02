##########trim_galore##########
trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 4 -j 5 --paired -o rep1_1.fq.gz rep1_2.fq.gz

##########fastqc##########
fastqc -o /fastqc -t 10 ../clean/*.fq.gz

##########bowtie2##########
bowtie2-build GCF_003254395.2_Amel_HAv3.1_genomic.fna ApisMellifera_index
bowtie2 --very-sensitive -X 700 -x $bowtie2_index -1 rep1_1.fq.gz -2 rep1_2.fq.gz -p 15 rep1 > rep1_bowtie2.txt | samtools sort -O bam -@ 10 -o rep1.raw.sort.bam
samtools index rep1.raw.sort.bam rep1.raw.sort.bam.bai
samtools flagstat rep1.raw.sort.bam > rep1.raw.sort.stat
sambamba markdup --overflow-list-size 600000 --tmpdir=./ -t 5 -r rep1.raw.sort.bam rep1.rmdup.bam
samtools sort -@ 5 rep1.rmdup.bam -o rep1.rmdup.sort.bam
samtools index rep1.rmdup.sort.bam rep1.rmdup.sort.bam.bai
samtools flagstat rep1.rmdup.sort.bam > rep1.rmdup.sort.stat
rm rep1.rmdup.bam
rm rep1.rmdup.bam.bai
samtools view -h -b -F 1804 -f 2 rep1.rmdup.sort.bam | samtools sort -O bam -@ 5 -o rep1.filtered.sort.bam
samtools index rep1.filtered.sort.bam rep1.filtered.sort.bam.bai
samtools flagstat rep1.filtered.sort.bam > rep1.filtered.sort.stat
samtools view -h -q 30 rep1.filtered.sort.bam | samtools sort -O bam -@ 5 -o rep1.last.sort.bam
samtools index rep1.last.sort.bam rep1.last.sort.bam.bai
samtools flagstat rep1.last.sort.bam > rep1.last.sort.stat
rm rep1.raw.sort.bam
rm rep1.raw.sort.bam.bai
rm rep1.rmdup.sort.bam
rm rep1.rmdup.sort.bam.bai
rm rep1.filtered.sort.bam
rm rep1.filtered.sort.bam.bai

##########macs2##########
macs2 callpeak -t rep1.last.sort.bam -c IgG.last.sort.bam --outdir ./ -f BAMPE -g 2.239e8 -n rep1 -B -q 0.05

##########TSS##########
bamCoverage -p 20 --normalizeUsing CPM -b rep1.last.sort.bam -o rep1.last.bw
computeMatrix reference-point --referencePoint TSS -p 20 -b 2500 -a 2500 -S rep1.last.bw -R genomic.gtf --skipZeros -o rep1_TSS.gz --outFileSortedRegions rep1_region_genes.bed

##########R##########
library(DiffBind)
rm(list = ls())
setwd("/DiffBind")

SampleID = c("Butyrate-1","Butyrate-2","Butyrate-3","GF-1","GF-2","GF-3")
Condition = c(rep("Treatment", 3), rep("Control", 3))
Replicate = c(1:3, 1:3)
PeakCaller = rep("narrow", 6)
bam_file_path <- "/DiffBind/"
bamReads <- c(
  paste(bam_file_path, "Butyrate-1_last.sort.bam", sep = "/"), 
  paste(bam_file_path, "Butyrate-2_last.sort.bam", sep = "/"),
  paste(bam_file_path, "Butyrate-3_last.sort.bam", sep = "/"),
  paste(bam_file_path, "GF-1_last.sort.bam", sep = "/"),
  paste(bam_file_path, "GF-2_last.sort.bam", sep = "/"),
  paste(bam_file_path, "GF-3_last.sort.bam", sep = "/")
)

peak_file_path <- "/DiffBind/"
Peaks <- c(
  paste(peak_file_path, "Butyrate-1_peaks.narrowPeak", sep = "/"),
  paste(peak_file_path, "Butyrate-2_peaks.narrowPeak", sep = "/"),
  paste(peak_file_path, "Butyrate-3_peaks.narrowPeak", sep = "/"),
  paste(peak_file_path, "GF-1_peaks.narrowPeak", sep = "/"),
  paste(peak_file_path, "GF-2_peaks.narrowPeak", sep = "/"),
  paste(peak_file_path, "GF-3_peaks.narrowPeak", sep = "/")
)
samples <- data.frame(SampleID,Condition,Replicate,bamReads,Peaks,PeakCaller)

H3K27ac <- dba(sampleSheet = samples, minOverlap = 1)
dba.plotHeatmap(H3K27ac)
save(H3K27ac,file="H3K27acOriginal.data")

H3K27ac_count <- dba.count(H3K27ac,minOverlap = 1)
save(H3K27ac_count,file="H3K27acCount.data")
H3K27ac_count_norm <- dba.normalize(H3K27ac_count,normalize=DBA_NORM_LIB)
save(H3K27ac_count_norm,file="h3k27normlized.data")
###
library(edgeR)
library(DiffBind)
library(profileplyr)

H3K27ac_diff_contrast <- dba.contrast(H3K27ac_count_norm,
                      group1=H3K27ac_count_norm$masks$Treatment,
                      group2=H3K27ac_count_norm$masks$Control,
                      name1="Butyrate", 
                      name2="GF",
                      minMembers=2) 

H3K27ac_diff_result <- dba.analyze(H3K27ac_diff_contrast,method=DBA_ALL_METHODS)
dba.show(H3K27ac_diff_result, bContrasts=TRUE)

comp1.deseq2<- dba.report(H3K27ac_diff_result,th = 1,method=DBA_DESEQ2, contrast = 1, bCounts = TRUE, bNormalized = TRUE)
out_of_deseq2 <- as.data.frame(comp1.deseq2)
colnames(out_of_deseq2)[colnames(out_of_deseq2) == "Fold"] <- "log2FC"
write.table(out_of_deseq2,"GF_vs_Butyrate_deseq2.txt",col.names = T,row.names = F,sep = "\t",quote = F)

out_of_deseq2$threshold[out_of_deseq2$p.value >= 0.05 | is.na(out_of_deseq2$p.value)] = "None"
out_of_deseq2$threshold[out_of_deseq2$p.value < 0.05 & (out_of_deseq2$log2FC >= -log2(1.5) & out_of_deseq2$log2FC <= log2(1.5))] = "None"
out_of_deseq2$threshold[out_of_deseq2$p.value < 0.05 & out_of_deseq2$log2FC > log2(1.5)] = "Up"
out_of_deseq2$threshold[out_of_deseq2$p.value < 0.05 & out_of_deseq2$log2FC < -log2(1.5)] = "Down"
write.table(out_of_deseq2,"GF_vs_Butyrate_allgene.txt",col.names = T,row.names = F,sep = "\t",quote = F)

up_down_data <- out_of_deseq2[out_of_deseq2$threshold %in% c("Up", "Down"), ]
write.table(up_down_data,"GF_vs_Butyrate_UpDowngene.txt",col.names = T,row.names = F,sep = "\t",quote = F)
###
library(ChIPseeker)
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("genomic.gtf")
diff_peak <- read.delim("GF_vs_Butyrate_UpDowngene.txt", header = TRUE)
diff_peak_bed <- diff_peak[,]
write.table(diff_peak_bed,"GF_vs_Butyrate_deseq2.bed",sep = "\t",quote = F,col.names = F,row.names = F)
diff_peak_bed <- readPeakFile("GF_vs_Butyrate_deseq2.bed")
diffpeak_Anno <- annotatePeak(diff_peak_bed, 
                              level = "gene", 
                              tssRegion = c(-3000, 3000), 
                              TxDb = txdb)
diffpeak_Anno_df <- as.data.frame(diffpeak_Anno)

write.table(diffpeak_Anno_df, 
            "Chipseeker_GF_vs_Butyrate.txt", 
            sep = "\t", 
            quote = FALSE, 
            col.names = TRUE, 
            row.names = FALSE)
###
library(rtracklayer)
library(dplyr)

diffpeak_Anno <- read.delim("Chipseeker_GF_vs_Butyrate.txt",header=T)
gtf <- rtracklayer::import('genomic.gtf')
gtf <- as.data.frame(gtf)
gtf <- dplyr::select(gtf, c(gene, gene_id)) 
gtf <- unique(gtf) 
diffpeak_Anno <- merge(diffpeak_Anno, gtf, by.x = "geneId", by.y = "gene_id")
amel_gene <- read.csv("Amel_gene.csv", header = TRUE)
merged_data <- merge(diffpeak_Anno, amel_gene,by.x = "gene", by.y = "Symbol", all.x = TRUE)

write.csv(merged_data, "GFvsButyrate_UpDowngeneDescription.csv", row.names = FALSE)
##########fastp##########
fastp -i rep1_1.fq.gz -o rep1_1.clean.fq.gz -I rep1_2.fq.gz -O rep1_2.clean.fq.gz -q 20 -u 10 -w 10 --adapter_sequence AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGCTCTTCCGATCT --adapter_sequence_r2 CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC -h rep1.html -j rep1.json

##########fastqc##########
fastqc rep1_1.clean.fq.gz rep1_2.clean.fq.gz -o ../fastqc -t 10

##########hisat2##########
nohup hisat2-build -p 10 GCF_003254395.2_Amel_HAv3.1_genomic.fna ApisMellifera_index &
hisat2 -p 10 --phred33 --dta-cufflinks --no-mixed --no-discordant --summary-file rep1_summary.log -x ../ApisMellifera_index -1 rep1_1.clean.fq.gz -2 rep1_2.clean.fq.gz -S rep1.sam
samtools view -bS rep1.sam -o rep1.bam
samtools sort rep1.bam -o rep1_sorted.bam -@ 10

##########stringtie#########
mv rep1_sorted.bam rep1.bam
stringtie -A rep1_gene_abund.tab -C rep1_cov_refs.gtf -p 10 -G genomic.gff -o ./rep1.gtf -e -b rep1_Ballgown rep1.bam;

##########R##########
###DESeq2###
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)

count_data = read.csv(file = "gene_count.csv", header = T)
rownames(count_data)=count_data[,3]
count_data_1=count_data[,5:ncol(count_data)]
count_data_2=count_data_1[rowSums(count_data_1) != 0,]
col_data = read.table(file = "0.col_data.txt", header = T, sep = "\t")
colnames(count_data_2) == col_data$id 
gene_decs = read.csv(file = "Amel_gene.csv",header = T)
colnames(gene_decs) = c("geneid","Symbol","Description")

col_data$condition=as.factor(col_data$condition)
dds = DESeqDataSetFromMatrix(countData = count_data_2, 
                             colData = col_data,
                             design = ~ condition)
colData(dds)
dds = DESeq(dds)
res <- results(dds,contrast = c("condition","Butyrate","GF"),alpha = 0.05)

sink("GF-Butyrate_fatbody_DESeq2_summary.txt")
mcols(res, use.names=TRUE)
summary(res)
sink()

res <- res[order(res$padj),]
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.table(normalized_counts, file="GF-Butyrate_fatbody_DESeq2.normalized.xls",
            quote=F, sep="\t", row.names=T, col.names=T)

resdata_total <- merge (as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
colnames(resdata_total)[1] = "geneid"
resdata_total_desc = merge(as.data.frame(resdata_total),as.data.frame(gene_decs),by.x="geneid",by.y="geneid",all.x=T,sort=FALSE)
new_order <- c("geneid", "Symbol", "Description", setdiff(colnames(resdata_total_desc), c("geneid", "Symbol", "Description")))
merged_data <- resdata_total_desc[, new_order]
merged_data %>% 
  mutate(group = case_when(
    log2FoldChange >= 1 & padj <= 0.05 ~ "Up",
    log2FoldChange <= -1 & padj <= 0.05 ~ "Down",
    TRUE ~ "NotChange"
  )) -> merged_data_1
write.csv(merged_data_1, file="GF-Butyrate_fatbody_total.csv",row.names = F)
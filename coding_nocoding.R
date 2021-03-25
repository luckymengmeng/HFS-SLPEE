rm(list = ls())
# setwd("/Users/mac/Downloads/exp9.0/Py_R_data_process/luad_data/1_learn/") #文件路径
# 以下代码将近6万维的转录组数据，分为编码蛋白质基因、长分编码RNA，miRNA和剩余的非编码基因。
gene <- read.table("../dataset/bio_contrast/TCGA_BRCA/TCGA_BRCA_logFC3_14/raw_data/gene_symbol_preprocess.txt",
                   header=T,sep="\t",check.names=F) #文件夹名字
names <-gene$id
# names
miRNA <- grep("\\|miRNA$", names)
mRNA <- grep("\\|protein_coding$", names)
lncRNA <- grep("\\|sense_overlapping|\\|lincRNA|\\|3prime_overlapping_ncRNA|\\|processed_transcript|antisense|\\|sense_intronic$", names)
other <- grep("\\|protein_coding|\\|sense_overlapping|\\|lincRNA|3prime_overlapping_ncRNA|\\|processed_transcript|\\|antisense|sense_intronic|\\|miRNA$", names)

miRNA_result <- subset(gene,grepl("\\|miRNA$", names))
mRNA_result <- subset(gene,grepl("\\|protein_coding$", names))
lncRNA_result <- subset(gene,grepl("\\|sense_overlapping|\\|lincRNA|\\|3prime_overlapping_ncRNA|\\|processed_transcript|antisense|\\|sense_intronic$", names))
other_result <- subset(gene,!grepl("\\|protein_coding|\\|sense_overlapping|\\|lincRNA|3prime_overlapping_ncRNA|\\|processed_transcript|\\|antisense|sense_intronic|\\|miRNA$", names))

write.table(miRNA_result,file="../dataset/bio_contrast/TCGA_BRCA/TCGA_BRCA_logFC3_14/raw_data/miRNA.txt",sep="\t",quote=F,row.names = F)
write.table(other_result,file="../dataset/bio_contrast/TCGA_BRCA/TCGA_BRCA_logFC3_14/raw_data/other.txt",sep="\t",quote=F,row.names = F)
write.table(mRNA_result,file="../dataset/bio_contrast/TCGA_BRCA/TCGA_BRCA_logFC3_14/raw_data/mRNA.txt",sep="\t",quote=F,row.names = F)
write.table(lncRNA_result,file="../dataset/bio_contrast/TCGA_BRCA/TCGA_BRCA_logFC3_14/raw_data/lncRNA.txt",sep="\t",quote=F,row.names = F)


# # # 以下是将近6万维的数据分为编码RNA和非编码RNA
# rm(list = ls())
# # setwd("/Users/mac/Downloads/exp/test/datasets/luad") #文件路径
# gene <- read.table("/Users/xujunlin/Downloads/bio_data/dataset/bio_contrast/TCGA_KIRC/TCGA_KIRC_logFC3_14/ncRNA/gene_symbol_preprocess.txt",header=T,sep="\t",check.names=F) #文件夹名字
# # 改下id,根据自己数据的第一行第一个字符串是什么就改成什么
# names <-gene$id
# names
# grep("\\protein_coding$", names)
# grepl("\\protein_coding$", names)
# coding_result <- subset(gene,grepl("\\protein_coding$", names))
# nocoding_result <- subset(gene,!grepl("\\protein_coding$", names))
# write.table(coding_result,file="/Users/xujunlin/Downloads/bio_data/dataset/bio_contrast/TCGA_KIRC/TCGA_KIRC_logFC3_14/ncRNA/coding_RNA.txt",sep="\t",quote=F,row.names = F)
# write.table(nocoding_result,file="/Users/xujunlin/Downloads/bio_data/dataset/bio_contrast/TCGA_KIRC/TCGA_KIRC_logFC3_14/ncRNA/nocoding_RNA.txt",sep="\t",quote=F,row.names = F)

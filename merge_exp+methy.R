rm(list = ls())
# 设置工作目录
setwd("/Users/xujunlin/Downloads/bio_data/dataset/")
# 甲基化文件
methyFile="../dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/Cpg/all_CpgDiff.txt"      
# 表达文件
expFile="../dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneExp/normalizeExp_diff.txt"                         

methy = read.table(methyFile, row.names=1 ,header=T,sep="\t",check.names=F)
RNA = read.table(expFile, row.names=1 ,header=T,sep="\t",check.names=F)
colnames(methy)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(methy))
colnames(RNA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(RNA))
rownames(methy)=paste(rownames(methy),"methy",sep="|")
rownames(RNA)=paste(rownames(RNA),"exp",sep="|")
sameSample=intersect(colnames(methy),colnames(RNA))
merge=rbind(id=sameSample,methy[,sameSample],RNA[,sameSample])
write.table(merge,file="../dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneCpg/merge_CpgExp.txt",
            sep="\t",quote=F,col.names=F)
name_merge <- union(rownames(methy),rownames(RNA))

write.table(name_merge, file = '../dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneCpg/merge_CpgExp_name.txt',
            sep = '\t',row.names = F,col.names = F,quote = F)

updata_methy = rbind(id=sameSample,methy[,sameSample])
updata_RNA = rbind(id=sameSample,RNA[,sameSample])
# 存储数据
write.table(updata_methy,file="../dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneCpg/Cpg_.txt",
            sep="\t",quote=F,col.names=F)
write.table(updata_RNA,file="../dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneCpg/exp_.txt",
            sep="\t",quote=F,col.names=F)
# 获取标签
# 根据样本名获取类别标签
# 01-09表示肿瘤，10-29表示正常对照组
sameSample_label <- grepl("\\-01|-02|-03|-04|-05|-06|-07|-08|-09",sameSample)
# 查看样本的类别标签
sameSample_label
# 将类别标签转为字符型 
sameSample_label <- as.character(sameSample_label) #标签
write.table(sameSample_label,file="../dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneCpg/sameSample_label.txt",
            sep="\t",quote=F,row.names=F)

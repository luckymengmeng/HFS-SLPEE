rm(list = ls())
options(stringsAsFactors = F)
timestart <- Sys.time()
# # 设置工作路径
setwd("/Users/xujunlin/Downloads/bio_data/dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneMethy/")
# 输入文件 行是属性 列是样本
inputFile="geneMethy_process.txt" 
# 正常样品的数目
normalNum=96
# 癌症样品的数目
tumorNum=792
fdrFilter=0.05    
logFCfilter=0.5                          

library(limma)
library(data.table)
outTab=data.frame()
grade=c(rep(1,normalNum),rep(2,tumorNum))
Type=c(rep("Normal",normalNum),rep("Tumor",tumorNum))
rt <- fread(inputFile,sep="\t",header=T)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data) # 处理基因有重复名字
data=data[rowMeans(data)>0,] # 去掉基因在所有样本中都为0的

# 矫正数据
data=normalizeBetweenArrays(data)
normalData=cbind(id=row.names(data),data)
fwrite(normalData,file = "normalizeCpg0.txt",
       sep = "\t",row.names = F)
# write.table(normalData,file="/Users/xujunlin/Downloads/bio_data/dataset/bio_contrast/TCGA_BRCA/TCGA_BRCA_logFC3_14/cpg/normalizeMethy.txt",
#             sep="\t",row.names=F,quote=F) # 矫正后的文件

# 差异分析
for(i in row.names(data)){
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)  # 差异分析

  normalGeneMeans=mean(data[i,1:normalNum])
  tumorGeneMeans=mean(data[i,(normalNum+1):ncol(data)])
  logFC=log2(tumorGeneMeans)-log2(normalGeneMeans)  
  pvalue=wilcoxTest$p.value
  normalMed=median(data[i,1:normalNum])
  tumorMed=median(data[i,(normalNum+1):ncol(data)])
  diffMed=tumorMed-normalMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		  outTab=rbind(outTab,
		               cbind(gene=i,
		               normalMean=normalGeneMeans,
		               TumorMean=tumorGeneMeans,
		               logFC=logFC,
		               pValue=pvalue))
	 }
}

# 对P值进行矫正
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

# 输出所有基因的甲基化差异情况
# write.table(outTab,file="/Users/xujunlin/Downloads/bio_data/dataset/bio_contrast/TCGA_BRCA/TCGA_BRCA_logFC3_14/cpg/allGene.xls",
#             sep="\t",row.names=F,quote=F)

# 输出差异甲基化的基因
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
# write.table(outDiff,file="/Users/xujunlin/Downloads/bio_data/dataset/bio_contrast/TCGA_BRCA/TCGA_BRCA_logFC3_14/cpg/diff.xls",
#             sep="\t",row.names=F,quote=F)
fwrite(outDiff,file = "diffCpg0.xls",
       sep = "\t",row.names = F)

timeend <- Sys.time()
runtime <- timeend - timestart
print(runtime)





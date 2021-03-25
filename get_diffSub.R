# dt：行为属性，列为样本的数据 
# sl：一列属性
# re：根据该列，得到子矩阵数据
# 注：保证得到的子集和diff_logFC1_fdr0.05_name.txt中的名字次序一致。
# clean workspace
rm(list = ls())
options(stringsAsFactors = F)
nor_exp <- fread('/Users/xujunlin/Downloads/bio_data/dataset/bio_contrast/TCGA_BRCA/TCGA_BRCA_logFC3_14/cpg/Cpg7/normalizeCpg7.txt',
            sep="\t",header=T)
diffName_exp <- read.table('/Users/xujunlin/Downloads/bio_data/dataset/bio_contrast/TCGA_BRCA/TCGA_BRCA_logFC3_14/cpg/Cpg7/diffCpg7.txt',
                 header = F,sep = '\t',check.names = F)
diff_nor_exp <- nor_exp[(nor_exp$id %in% diffName_exp$V1),]
# ##如果re的行数和sl的行数不一致，一般是因为有3-mar这样的基因名
# ##查看有多少个不一致的，并且显示行号
# which(sl$V1 %in% re$id == FALSE)
write.table(diff_nor_exp, file = '/Users/xujunlin/Downloads/bio_data/dataset/bio_contrast/TCGA_BRCA/TCGA_BRCA_logFC3_14/cpg/Cpg7/normalizeExp_diff7.txt',
            sep = '\t',row.names = F,col.names = T,quote = F)

# dt：行为样本，列为属性的数据 
# sl：属性，为一行
# re：根据该行，得到子矩阵数据
# 最后虽得到子集，但不保证子集和特征顺序是一致的
#clean workspace
# rm(list = ls())
# options(stringsAsFactors = F)
# # ##read data
# # ##导入数据,数据框
# dt <- read.table('../dataset/TCGA_BRCA/methyGeneMerge_final/merge_methyExp_T_Nor.txt',
#                  header = T,sep = '\t',check.names = F)
# sl <- read.table('../dataset/TCGA_BRCA/methyGeneMerge_final/merge_100.txt',
#                  header = T,sep = '\t',check.names = F)
# d <- colnames(dt)
# s <- colnames(sl)
# 
# re <- dt[,(d %in% s)]
# id <- dt$id
# re_f = cbind(id,re)
# # ##如果re的行数和sl的行数不一致，一般是因为有3-mar这样的基因名
# # ##查看有多少个不一致的，并且显示行号
# # which(sl$V1 %in% re$id == FALSE)
# write.table(re_f, file = '../dataset/TCGA_BRCA/methyGeneMerge_final/fea_100.txt',
#             sep = '\t',row.names = F,col.names = T,quote = F)
# 
# ?match

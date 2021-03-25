rm(list = ls())
options(stringsAsFactors = F)
# 读入每折生成的特征index,即最优结果对应的特征lindex
# 每折生成的特征索引是从0开始
feaname_index = read.table("/Users/xujunlin/Downloads/bio_data/dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneCpg/merge_mrmr21/feature.txt",
                      sep="",check.names=F,comment.char = "]")
# 去掉第一列所有字符的首字符"["
v <- c(1:10)
for(i in v){
   f = substr(feaname_index$V1[i], start = 2, stop = 8)
   feaname_index$V1[i] = f
}
# 每行组合起来
one_fold <- feaname_index[1,]
two_fold <- feaname_index[2,]
three_fold <- feaname_index[3,]
four_fold <- feaname_index[4,]
five_fold <- feaname_index[5,]
six_fold <- feaname_index[6,]
seven_fold <- feaname_index[7,]
eight_fold <- feaname_index[8,]
nine_fold <- feaname_index[9,]
ten_fold <- feaname_index[10,]
# 按列组合
merge2 <- cbind(one_fold,two_fold)
merge3 <- cbind(merge2,three_fold)
merge4 <- cbind(merge3,four_fold)
merge5 <- cbind(merge4,five_fold)
merge6 <- cbind(merge5,six_fold)
merge7 <- cbind(merge6,seven_fold)
merge8 <- cbind(merge7,eight_fold)
merge9 <- cbind(merge8,nine_fold)
merge10 <- cbind(merge9,ten_fold)
# 将每折的特征索引组合为一行
fea_merge = as.character(merge10[1,])
table(duplicated(fea_merge))
# 不重复的索引分布
!duplicated(fea_merge)
# 去掉重复的索引
fea_quchong = merge10[,!duplicated(fea_merge)]
# 读入特征名字
feaname_name = read.table("/Users/xujunlin/Downloads/bio_data/dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneCpg/merge_mrmr21/merge_CpgExp_name.txt",
                           sep="",check.names=F)
# 交集转换为数字
fea_quchong_int <- as.numeric(fea_quchong)
# 每个数字加1
fea_quchong_update <- fea_quchong_int + 1
# 存储交集的名字
intersect_feaName <- feaname_name[fea_quchong_update,]
# 全集转换为数字
merge10_int <- as.numeric(merge10)
# 每个数字加1
merge10_update <- merge10_int + 1
# 存储交集的名字
union_feaName <- feaname_name[merge10_update,]

# 保存
write.table(union_feaName, file = '/Users/xujunlin/Downloads/bio_data/dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneCpg/merge_mrmr21/union_feaName.txt',
            sep = '\t',row.names = F,col.names = F,quote = F)
write.table(intersect_feaName, file = '/Users/xujunlin/Downloads/bio_data/dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneCpg/merge_mrmr21/intersect_feaName.txt',
            sep = '\t',row.names = F,col.names = F,quote = F)

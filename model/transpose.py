# transpose the data
def T(fileData,fileData_T):

    filename1 = fileData
    f = open(filename1)

    temp = f.readlines()
    temp_1 = []

    for i in range(0, len(temp)):
        data = temp[i].strip().split('\t')
        temp_1.append(data)
    print(len(temp_1))
    print(len(temp_1[0]))
    le = len(temp_1)
    le1 = len(temp_1[0])

    filename2 = fileData_T
    ff = open(filename2, 'w')

    for i in range(0, len(temp_1[0])):
        for j in range(0, len(temp_1)):
            if(j < le-1):
                ff.write(temp_1[j][i])
                ff.write('\t')
            else:
                ff.write(temp_1[j][i])
        if i < le1 - 1:
            ff.write("\n")
    ff.close()
    return filename2

file_Dir = '/Users/xujunlin/Downloads/bio_data/dataset/TCGA_BRCA/TCGA_BRCA_logFC3_logFC0.5_14/geneCpg/'
# file_Dir = '../../dataset/bio_contrast/TCGA_LUAD/TCGA_LUAD_logFC3_20/geneMethy/'
# file_Dir = '../../dataset/bio_contrast/TCGA_KIRC/TCGA_KIRC_logFC3_21/geneMethy/'
file_data = file_Dir + 'merge_CpgExpT_nor.txt'
file_dataT = file_Dir + 'merge_CpgExp_nor.txt'
T(file_data,file_dataT)

def loadDataset(filename, Data_X=[], Data_Y=[]):
    f = open(filename)
    lines = f.readlines()

    dataMat = []
    N = len(lines)
    for i in range(0, N):
        data = lines[i].strip().split("\t")
        dataMat.append(data)

    for i in range(1, N):
        Y = dataMat[i][0].strip().split("-")
        if (Y[3][0] == '1' and Y[3][1] == '1'):
            Data_Y.append(-1)
        elif (Y[3][0] == '0' and Y[3][1] == '1'):
            Data_Y.append(1)
        else:
            Data_Y.append(1)
        for index in range(1, len(dataMat[i])):
            dataMat[i][index] = float(dataMat[i][index])
        del dataMat[i][0]
        Data_X.append(dataMat[i])

    normal = 0
    tumor = 0
    for i in range(0, N - 1):
        if (Data_Y[i] == -1):
            normal = normal + 1
        if (Data_Y[i] == 1):
            tumor = tumor + 1

    print("样本总数:"+str(N-1))
    print("正常样本数:"+str(normal))
    print("癌症样本数:"+str(tumor))
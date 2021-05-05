# -*- coding: UTF-8 -*- 
'''
@Time       :   2020/2/10 5:39 下午
@Author     :   Meng Yajie
@FileName   :   1_main_KIRC.py
@Software   :   PyCharm
'''

import load
import os
import numpy as np
#import transpose
import text_save
from sklearn.model_selection import StratifiedKFold
from sklearn import preprocessing
from skfeature.function.information_theoretical_based import MRMR
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from mlxtend.classifier import StackingCVClassifier
from sklearn.metrics import confusion_matrix, classification_report
from sklearn import tree
import xgboost as xgb

# the dataset
# 初始化各维度的实验结果衡量指标
accuracy_file = []
precision_file = []
recall_file = []
f1_score_file = []
specificity_file = []

# 设置目录
data_Dir = "../../dataset/TCGA_KIRC/TCGA_KIRC_logFC3_logFC0.5_14/methyGeneMerge/"
# 数据，离散格式，行是样本，列为特征
source_file_disc = data_Dir + "merge_methyExp_disc.txt"
# 原始融合数据
# 行是属性，列是样本
source_file = data_Dir + "merge_methyExp.txt"
# 转置后，行是样本，列是属性
source_fileT = data_Dir + "merge_methyExpT.txt"
# 转置原始数据，行是样本，列是特征
# transpose.T(source_file,source_fileT)

# 设置维度的变化范围，步长设置为1
for num_fea in range(16,17,1):
    print(num_fea)
    # # 设置目录
    # data_Dir = "../../dataset/TCGA_LUAD/TCGA_LUAD_logFC1.5_logFC0.5/methyGeneMerge2/"

    # 新建各个维度目录
    nDim_DIR = data_Dir + "/merge_mrmr" + str(num_fea) + "/"
    if (os.path.exists(nDim_DIR) != True):
        os.mkdir(nDim_DIR)

    # 导入离散数据
    # get the feature data and the label
    X = []
    y = []
    load.loadDataset(source_file_disc, X, y)
    X = np.array(X)
    y = np.array(y)

    # 导入原始数据
    # get the feature data and the label
    X_raw = []
    y_raw = []
    load.loadDataset(source_fileT, X_raw, y_raw)
    X_raw = np.array(X_raw)
    y_raw = np.array(y_raw)

    # ten fold cross validation
    # BRCA: random_state = 14
    # KIRC: random_state = 14
    # LUAD: random_state = 20
    skf = StratifiedKFold(n_splits=10,random_state=14,shuffle=True)
    # get the number of fold
    n = skf.get_n_splits(X, y)

    # 初始化各指标
    TN = np.zeros(n,dtype=int)
    FP = np.zeros(n,dtype=int)
    FN = np.zeros(n,dtype=int)
    TP = np.zeros(n,dtype=int)

    for (train_index, test_index),foldNumber in zip(skf.split(X, y),range(n)):
        # print("TRAIN:", train_index, "TEST:", test_index)

        # 离散数据
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        print('训练集数量:', X_train.shape, '测试集数量:', X_test.shape)
        # print(y_train)
        # print(y_test)

        # 原始数据
        X_train_raw, X_test_raw = X_raw[train_index], X_raw[test_index]
        y_train_raw, y_test_raw = y_raw[train_index], y_raw[test_index]
        print('训练集数量:', X_train_raw.shape, '测试集数量:', X_test_raw.shape)

        # obtain the index of each feature on the training set which is selected
        idx, _, _ = MRMR.mrmr(X_train, y_train, n_selected_features=num_fea)
        print("the index of selected feature is: " + str(idx))

        # save the index of selected feature
        feature_file = nDim_DIR + "特征.txt"
        f0 = open(feature_file, 'a+')  # 基因
        f0.write(str(idx))
        f0.write('\n')
        f0.close()

        # obtain the dataset on the selected features
        # update the discrete train data and the discrete test data
        X_train_select_disc = X_train[:, idx[0:num_fea]]
        X_test_select_disc = X_test[:, idx[0:num_fea]]

        # obtain the dataset on the selected features
        # update the raw train data and the raw test data
        X_train_select_raw = X_train_raw[:, idx[0:num_fea]]
        X_test_select_raw = X_test_raw[:, idx[0:num_fea]]

        # 标准化数据方法一
        # min_max_scaler = preprocessing.MinMaxScaler()
        # X_train_minmax = min_max_scaler.fit_transform(X_train)
        # X_test_minmax = min_max_scaler.fit_transform(X_test)

        # 标准化数据方法二
        scaler = preprocessing.StandardScaler().fit(X_train_select_raw)
        # print(scaler.mean_)
        # print(scaler.scale_)
        X_train_scale_raw = scaler.transform(X_train_select_raw)
        X_test_scale_raw = scaler.transform(X_test_select_raw)

        # 设置随机种子，使结果可以复现
        RANDOM_SEED = 200
        # Stacked CV Classification and GridSearch
        # Initializing models
        # clf1 = KNeighborsClassifier(n_neighbors=5)
        # clf1 = GaussianNB()
        clf1 = tree.DecisionTreeClassifier(random_state=RANDOM_SEED)
        clf2 = RandomForestClassifier(n_estimators=50,
                                      oob_score=True,
                                      random_state=RANDOM_SEED)
        clf3 = SVC(kernel='rbf', probability=True, gamma='scale',random_state=RANDOM_SEED)
        # clf4 = GaussianNB()
        clf4 = AdaBoostClassifier(n_estimators=50,
                                  random_state=RANDOM_SEED)
        # clf4 = xgb.XGBClassifier(max_depth=6, n_estimators=100, num_round=5)

        # lr = LogisticRegression(C=0.1,max_iter=100,solver='lbfgs',random_state=RANDOM_SEED)
        # lr = SVC(kernel='rbf', C=0.1, probability=True, gamma='scale')
        lr = xgb.XGBClassifier(max_depth=6, n_estimators=100, num_round=5,random_state=RANDOM_SEED)

        # 设置use_probas = True，第一层单个模型预测结果是类别概率，而不是label
        sclf = StackingCVClassifier(classifiers=[clf1, clf2, clf3, clf3, clf4],
                                    use_probas=True,
                                    meta_classifier=lr,
                                    random_state=79)
        # 参数优化
        params = {
            # 错7个
            # 'kneighborsclassifier__n_neighbors': [5,7],
            'randomforestclassifier__n_estimators': [50, 100],
            # 'gradientboostingclassifier_learning_rate':[0.1,0.2,0.3],
            # 'adaboostclassifier__n_estimators': [100,250,500],
            # 'svmclassifier-1__C': [0.1, 0.25],
            'svc-1__C': [0.001,0.01,0.1],
            # 'svc-1__C': [x / 0.1 for x in range(1, 10)],
            'svc-2__gamma': [1.0,10.0,100.0],
            # 'meta-logisticregression__C': [0.01, 0.1, 1.0],
            # 'svc-1__C': [x / 0.1 for x in range(1, 10)],
            # 'svc-3__kernel': ['rbf','linear'],
            # 'meta_xgbclassifier__n_estimators': [100, 200]
            'meta_classifier__n_estimators':[100, 200, 300]
            # 'meta-logisticregression__C': [0.01,0.1,1.0]
        }

        # The scoring parameter: defining model evaluation rules
        # accuracy,average_precision,balanced_accuracy,f1,f1_micro,f1_macro,f1_weighted,recall,roc_auc
        # 优先选择f1,average_precision,roc_auc
        # 参数寻优
        grid = GridSearchCV(estimator=sclf,
                            param_grid=params,
                            cv=3,
                            refit=True,
                            verbose=1,
                            n_jobs=-1,
                            scoring="roc_auc")
        # 模型训练
        print("正在训练第" + str(foldNumber) + "折： ")
        grid.fit(X_train_scale_raw, y_train_raw)
        cv_keys = ('mean_test_score', 'std_test_score', 'params')
        # 打印训练过程中验证集的平均精度，最高精度
        # for r, _ in enumerate(grid.cv_results_['mean_test_score']):
        #     print("%0.3f +/- %0.2f %r"
        #           % (grid.cv_results_[cv_keys[0]][r],
        #              grid.cv_results_[cv_keys[1]][r] / 2.0,
        #              grid.cv_results_[cv_keys[2]][r]))

        print("正在保存第" + str(foldNumber) + "折最优模型的参数： ")
        print('Best parameters: %s' % grid.best_params_)
        print('Accuracy: %.2f' % grid.best_score_)

        # params = grid.best_params_
        # 存储每折最优模型参数在当前文件夹中
        optimalModel = nDim_DIR + "bestParams.txt"
        f3 = open(optimalModel, 'a+')
        f3.write(str(grid.best_params_))
        f3.write('\n')
        f3.close()

        # test the data
        print("正在测试第" + str(foldNumber) + "折： ")
        y_pred = grid.predict(X_test_scale_raw)

        # 得到相应每折的TN，FP，FN，TP
        conf = confusion_matrix(y_test_raw, y_pred)
        TN[foldNumber] = conf[0][0]
        FP[foldNumber] = conf[0][1]
        FN[foldNumber] = conf[1][0]
        TP[foldNumber] = conf[1][1]
        print(classification_report(y_test_raw, y_pred))
        print()
        print(confusion_matrix(y_test_raw, y_pred))
        print()
        print("第" + str(foldNumber) + "折结束： ")

        # 保存每折的预测结果
        final_pred = nDim_DIR + "最终预测结果.txt"
        f1 = open(final_pred, 'a+')
        f1.write(str(y_pred))
        f1.write('\n')
        f1.close()

        # 保存每折的混淆矩阵
        con_file = nDim_DIR + "混淆矩阵.txt"
        f2 = open(con_file, 'a+')  # 基因
        f2.write(str(confusion_matrix(y_test_raw, y_pred)))
        f2.write('\n')
        f2.close()
    # 打印十折的结果
    print("每折的TN是：")
    print(TN)
    print("每折的FP是：")
    print(FP)
    print("每折的FN是：")
    print(FN)
    print("每折的TP是：")
    print(TP)
    TN_sum = TN.sum()
    FP_sum = FP.sum()
    FN_sum = FN.sum()
    TP_sum = TP.sum()
    # print(TN_sum,FP_sum,FN_sum,TP_sum)
    # 最终的大矩阵
    conf_final = np.zeros(4,dtype=int).reshape(2,2)
    conf_final[0][0] = TN_sum
    conf_final[0][1] = FP_sum
    conf_final[1][0] = FN_sum
    conf_final[1][1] = TP_sum
    print("模型训练完毕，最终的混淆矩阵是：")
    print()
    print(str(conf_final))

    # 评价指标
    accu = (TN_sum + TP_sum)/(TN_sum + FP_sum + FN_sum + TP_sum)
    precision = TP_sum / (TP_sum + FP_sum)
    recall = TP_sum / (TP_sum + FN_sum)
    f1_score = (2*TP_sum) / (2*TP_sum + FP_sum + FN_sum)
    specificity = TN_sum / (TN_sum + FP_sum)
    print("accuracy: " + str(accu))
    print("precision: " + str(precision))
    print("recall: " + str(recall))
    print("f1_score: " + str(f1_score))
    print("specificity: " + str(specificity))

    # 保存十折的总混淆矩阵
    conf_final_file = nDim_DIR + "总混淆矩阵.txt"
    f4 = open(conf_final_file, 'a+')  # 基因
    f4.write(str(conf_final))
    f4.write('\n')
    f4.close()

    # 保存十折的结果衡量指标值
    result_final_file = nDim_DIR + "各个衡量指标值.txt"
    f5 = open(result_final_file, 'a+')  # 基因
    f5.write("accuracy: " + str(accu))
    f5.write('\n')
    f5.write("precision: " + str(precision))
    f5.write('\n')
    f5.write("recall : " + str(recall))
    f5.write('\n')
    f5.write("f1_score: " + str(f1_score))
    f5.write('\n')
    f5.write("specificity: " + str(specificity))
    f5.write('\n')
    f5.close()
    # 存储每折的各个指标值
    accuracy_file.append(accu)
    precision_file.append(precision)
    recall_file.append(recall)
    f1_score_file.append(f1_score)
    specificity_file.append(specificity)
# 合并
temp1 = np.vstack((accuracy_file, precision_file))
temp2 = np.vstack((temp1, recall_file))
temp3 = np.vstack((temp2, f1_score_file))
re = np.vstack((temp3, specificity_file))
# 存储
re_fileName = data_Dir + "各个维度的衡量指标值汇总表.txt"
text_save.text_save(filename=re_fileName,data=re)
f6 = open(re_fileName, 'a+')  # 基因
f6.write("第一行是各个维度的accuracy")
f6.write('\n')
f6.write("第二行是各个维度的precision")
f6.write('\n')
f6.write("第三行是各个维度的recall")
f6.write('\n')
f6.write("第四行是各个维度的f1_score")
f6.write('\n')
f6.write("第五行是各个维度的specificity")
f6.write('\n')
f6.close()

# 保存各个维度的accuracy
accuracy_result = data_Dir + "accuracy.txt"
f7 = open(accuracy_result,'w')
f7.write(str(accuracy_file))
f7.close()

# 保存各个维度的precision
precision_result = data_Dir + "precision.txt"
f8 = open(precision_result,'w')
f8.write(str(precision_file))
f8.close()

# 保存各个维度的recall
recall_result = data_Dir + "recall.txt"
f9 = open(recall_result,'w')
f9.write(str(recall_file))
f9.close()

# 保存各个维度的f1
f1_result = data_Dir + "f1-score.txt"
f10 = open(f1_result,'w')
f10.write(str(f1_score_file))
f10.close()

# 保存各个维度的specificity
specificity_result = data_Dir + "specificity.txt"
f11 = open(specificity_result,'w')
f11.write(str(specificity_file))
f11.close()

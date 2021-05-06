import matplotlib.pyplot as plt
# 说明字体
font2 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 15}
fea = [5,6,7,8,9,10,11,12,13,14,15,
       16,17,18,19,20,21,22,23,24,
       25,26,27,28,29,30]
acc =[0.98490128,0.99070848,0.98954704,0.99186992,0.99303136,0.99186992,
 0.98954704,0.99070848,0.99186992,0.99070848,0.99070848,0.99186992,
 0.99186992,0.99070848,0.98373984,0.9941928,0.99651568,0.9941928,
 0.99186992,0.99303136,0.98722416,0.99186992,0.99303136,0.99186992,
 0.99303136,0.98954704]
ax = plt.subplots(figsize=(7,6))
# 显示指定点
# plt.scatter(21,0.9965)
# # 显示坐标点横线、竖线
plt.vlines(21,min(acc),0.9965,colors="black",linestyles="dashed")
plt.hlines(0.9965,5,21,colors="black",linestyles="dashed")
# 显示坐标点坐标值
plt.text(21,0.9965,(21,float('%.4f'% 0.9965)),ha='left',va='baseline',fontsize = 15)
plt.text(21,0.9965,(21,float('%.4f'% 0.9965)),ha='left',va='baseline',fontsize = 15)
# plot
plt.plot(fea,acc,lw=3)
# 横坐标起止
plt.xlim(5,30)
plt.xlabel("Feature number",font2)
plt.ylabel("Accuracy",font2)
plt.title("Feature number vs. Accuracy of BRCA",font2)
plt.tick_params(labelsize = 12)
plt.savefig('../../dataset/plot_optimal_matrics/exp_result/brca_acc2.tiff', dpi=300)
# plt.show()

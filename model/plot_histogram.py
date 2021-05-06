import numpy as np
import matplotlib.pyplot as plt
# 图例字体
font1 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 27}
# 说明字体
font2 = {'family': 'Times New Roman',
         'weight': 'normal',
         'size': 38}
labels = ['BRCA_21', 'LUAD_12', 'KIRC_16']
# 画specificity
ncRNA = [90.3,98.3,94.4]
proteinGene = [96.5,88.1,97.2]
lncRNA = [93.8,91.5,94.4]
miRNA_low = [83.2,81.4,90.3]
merged = [100.0,100.0,100.0]
gene = [92.0,96.6,97.2]
methyl = [94.8,96.9,99.4]
proteinGene_methy = [98.80,90.48,95.83]

x = np.arange(len(labels))  # the label locations
width = 0.1  # the width of the bars

fig, ax = plt.subplots(figsize=(18,16)) # 18，16

# 颜色方案
rects1 = ax.bar(x, proteinGene, width,color="#FC84FC",label='mRNA',edgecolor='black',lw=0)
rects2 = ax.bar(x + width, miRNA_low, width,color="#FC0404",label='miRNA',edgecolor='black',lw=0)
rects3 = ax.bar(x + 2*width, lncRNA, width,color="#A00434",label='lncRNA',edgecolor='black',lw=0)
rects4 = ax.bar(x + 3*width, ncRNA, width,color="#F8CCB4",label='ncRNA',edgecolor='black',lw=0)
rects5 = ax.bar(x + 4*width, methyl, width,color="#FCF404",label='Methylation',edgecolor='black',lw=0)
rects6 = ax.bar(x + 5*width, gene, width,color="#7094BC",label='Transcriptomic',edgecolor='black',lw=0)
rects7 = ax.bar(x + 6*width, proteinGene_methy, width,color="#A0DCEC",label='mRNA+DM',edgecolor='black',lw=0)
rects8 = ax.bar(x + 7*width, merged, width,color="#28B44C",label='TDS',edgecolor='black',lw=0)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Specificity(%)',font2)
# ax.set_title('Specificity vs. different type of data on three cancer dataset',font2)
# 设置横坐标每个区域中的中心线
ax.set_xticks(x + 3.5*width)
ax.set_xticklabels(labels,font2)
# 画图例
# ax.legend()
plt.ylim(top=101, bottom=80)
plt.tick_params(width = 5,length = 10, labelsize = 30)

# def autolabel(rects):
#     """Attach a text label above each bar in *rects*, displaying its height."""
#     for rect in rects:
#         height = rect.get_height()
#         ax.annotate('{}'.format(height),
#                     xy=(rect.get_x() + rect.get_width() / 2, height),
#                     xytext=(0, 3),  # 3 points vertical offset
#                     textcoords="offset points",
#                     ha='center', va='bottom')
# autolabel(rects1)
# autolabel(rects2)
# autolabel(rects3)
# autolabel(rects4)
# autolabel(rects5)
# autolabel(rects6)
# autolabel(rects7)
# fig.tight_layout()
# 设置图例格式
# plt.legend(loc='upper left', fancybox=True, ncol=3,prop= font1)
plt.savefig('../../dataset/plot_optimal_matrics/exp_case1_histogram/Specificity.tiff',dpi=300)
plt.show()
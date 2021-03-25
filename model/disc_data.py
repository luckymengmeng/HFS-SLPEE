# -*- coding: UTF-8 -*- 
'''
@Time       :   2020/1/13 5:33 下午
@Author     :   Meng Yajie
@FileName   :   disc_data.py
@Software   :   PyCharm
'''

'''function disc = disc_matrix(matrix)

% matrix is a value matrix
% each row is a sample and each column is a variable/attribute/feature
%% discretize by using mean+/-alpha*std. Here let alpha=0.5 for simplify
m = mean(matrix);
s = std(matrix);
disc = zeros(size(matrix));
for i = 1: size(matrix,2)
  disc(matrix(:, i)>=m(i)+s(i)/2, i) = 2;
  disc(matrix(:, i)<=m(i)-s(i)/2, i) = -2;
end
'''

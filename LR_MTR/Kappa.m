
% testlabel标准图拉成一维向量
% predict_val即为你自己的分类结果图，也拉成一维向量
% nlabel为类别个数，如果是三类，nlabel=3
% 
% testlabel的标签和predict_val的标签一致，例如在testlabel中农田标记为1，则，在predict_val中分类正确的点也要标记为1，其他类似。
% 
% 你在用的时候，把选择的点，依次排成一个向量就行了。


function [Kappavalue] = Kappa(nlabel,testlabel,predict_val)

%% 说明 
% input:
% nlabel 类别数
% testlabel 测试样本标签 （标准图的标签）
% predict_val 预测测试样本标签
% output: Kappavalue Kappa系数值
% labels=[0 102 204];
% term=0;
% for i=1:nlabel
%     term=term+length(find(testlabel==labels(i)))*length(find(predict_val==labels(i)));
% end
% Kappavalue=(length(testlabel)*(length(find(predict_val==testlabel)))-term)/(length(testlabel)^2-term);
[n m]=size(testlabel);
pe=0;
labels=[1 2 3 ];
for i=1:nlabel
    pe=pe+length(find(testlabel==labels(i)))*length(find(predict_val==labels(i)));
end
pe=pe/(n*n);
pa=length(find(testlabel==predict_val));
pa=pa/n;
Kappavalue=(pa-pe)/(1-pe);
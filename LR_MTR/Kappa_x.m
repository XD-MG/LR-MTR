function [Kappavalue] = Kappa_x(nlabel,testlabel,predict_val)

%% ˵�� 
% input:
% nlabel �����
% testlabel ����������ǩ ����׼ͼ�ı�ǩ��
% predict_val Ԥ�����������ǩ
% output: Kappavalue Kappaϵ��ֵ
% labels=[0 102 204];
% term=0;
% for i=1:nlabel
%     term=term+length(find(testlabel==labels(i)))*length(find(predict_val==labels(i)));
% end
% Kappavalue=(length(testlabel)*(length(find(predict_val==testlabel)))-term)/(length(testlabel)^2-term);
[n m]=size(testlabel);
pe=0;

labels=[];
for i=1:1:nlabel
   labels=[labels,i] ;
end
% labels=[1 2 3 4 5 6 7 8 9];
for i=1:nlabel
    pe=pe+length(find(testlabel==labels(i)))*length(find(predict_val==labels(i)));
end
pe=pe/(n*n);
pa=length(find(testlabel==predict_val));
pa=pa/n;
Kappavalue=(pa-pe)/(1-pe);
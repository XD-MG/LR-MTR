function [rate, rate1, rate2, rate3, rate4,rate5, testlabel] = KNN_Classfier(train, trainlabel, test,true_testlabel, K, NUM)
%K Nearest Neighbor Classfier
%train, test: each column is a data
%trainlabel, true_testlabel: row vectors containing the labels
n1=NUM(1);n2=NUM(2);n3=NUM(3);n4=NUM(4);n5=NUM(5);
aa=sum(train.*train,1); bb=sum(test.*test,1); ab=train'*test;
dist = sqrt(repmat(aa',[ 1 numel(bb)]) + repmat(bb, [numel(aa) 1]) - 2*ab);

% make sure result is all real
dist = real(dist);

[B,IX] = sort(dist);

testlabel = zeros(1,size(test,2));

for i  = 1:size(test,2)
    minindex = IX(1:K,i);
    
    neighborlabels = trainlabel(minindex);%��ȡ�����������ǩ
    
    class_ids = unique(neighborlabels);%��ȡ�����ǩ
    C = numel(class_ids);%������ĸ���
    neighbor_nums_per_class = zeros(1,C);
    for j = 1:C
        neighbor_nums_per_class(j) = sum(neighborlabels == class_ids(j));%����ÿ����Ľ���������
    end
    [temp,id_id] = max(neighbor_nums_per_class);%������������Ǹ���
    if sum(neighbor_nums_per_class==temp)==1
        testlabel(i) = class_ids(id_id);%�õ���������������,���ֻ��һ����������K����ѵ��������
    else
        %               fprintf('�˵�����������ͬ�������������')
        %               i
        testlabel(i)=trainlabel(IX(1,i)); %���ĳ������������ͬ�Ľ��������㣬�������ڵ����ڵ��Ǹ����
    end
end
rate=(sum(testlabel==true_testlabel)/length(true_testlabel))*100;
testlabel=testlabel;
num=0; num1=0; num2=0; num3=0; num4=0;
for i  = 1:size(test,2)
    AA = testlabel(i);
    BB = true_testlabel(i);
    if AA == BB && BB == 1
        num = num+1;
    end
    
    if AA == BB && BB == 2
        num1 = num1+1;
    end
    
    if AA == BB && BB == 3
        num2 = num2+1;
    end
    
    if AA == BB && BB == 4
        num3 = num3+1;
    end
    
    if AA == BB && BB == 5
        num4 = num4+1;
    end
end

rate1=(num/n1)*100;
rate2=(num1/n2)*100;
rate3=(num2/n3)*100;
rate4=(num3/n4)*100;
rate5=(num4/n5)*100;

end






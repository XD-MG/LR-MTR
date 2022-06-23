clear;
clc;

load PolSAR_feature_rectangle.mat;
disp('load data done.');
groundtruth=imread('Got_label_picture.bmp');
[rows, cols, ~] = size(groundtruth);


T_vector=PolSAR_feature_rectangle';%�����T_vector������涨Ϊÿһ��Ϊһ��������ÿһ�д���һά������
orgin_vector=T_vector(1:9,:);%ԭʼ9ά���� %85
% orignDB_vector=T_vector(10:19,:);%ԭʼ9ά����ȡģȡ��λ ������ܲ� %21
Freeman_vector=T_vector(20:22,:);%Freeman 3���� %85.4
Vanzyl_vector=T_vector(23:25,:);%Vanzyl 3���� %82.29
Yamaguchi_vector=T_vector(26:29,:);% Yamaguchi 4���� %82.91
% T_vector=T_vector(30:33,:);% Neumann 4����  �ܲ� %21
Krogager_vector=T_vector(34:36,:);% Krogager 3���� %87 ���
% TSVM_vector=T_vector(37:52,:);% TSVM 16����   �ܲ� %28
% H_alpha_vector=T_vector(53:93,:);% H/alpha 41����  �ܲ� %28
% T_vector=orgin_vector;
% T_vector=[Freeman_vector;Vanzyl_vector;Yamaguchi_vector;Krogager_vector];
T_vector=[orgin_vector;Freeman_vector;Vanzyl_vector;Yamaguchi_vector;Krogager_vector];
T_vector=T_vector';






[sample_number,deep] = size(T_vector);
PolSAR_feature_mat=reshape(T_vector,rows,cols,deep );

M=9;
t=floor(M/2);
PolSAR_feature_mat_sp=padarray(PolSAR_feature_mat,[t t],'symmetric','both');
PolSAR_feature_tensor=ones(M,M,deep,sample_number);

num=1;
%% **********************************************************************************************
%^^^^^^^^Ϊ����reshape ����ƥ�䣬����rows*colsתΪsample_numberʱ��Ҫ���а������У�һ��һ�е�������������ʾ^^^^^^^^^
%% **********************************************************************************************

tic
for j=1+t:cols+t
    for i=1+t:rows+t   
        block=PolSAR_feature_mat_sp(i-t:i+t,j-t:j+t,:);
        PolSAR_feature_tensor(:,:,:,num)=block; 
        num=num+1;
    end
    disp(['��',num2str(j-t),'��'])
end


save('PolSAR_feature_tensor.mat','PolSAR_feature_tensor','-v7.3');
save('T_vector.mat');
toc










T_vector1=PolSAR_feature_rectangle';%�����T_vector������涨Ϊÿһ��Ϊһ��������ÿһ�д���һά������
orgin_vector=T_vector1(1:9,:);%ԭʼ9ά���� %85
% orignDB_vector=T_vector(10:19,:);%ԭʼ9ά����ȡģȡ��λ ������ܲ� %21
Freeman_vector=T_vector1(20:22,:);%Freeman 3���� %85.4
Vanzyl_vector=T_vector1(23:25,:);%Vanzyl 3���� %82.29
Yamaguchi_vector=T_vector1(26:29,:);% Yamaguchi 4���� %82.91
% T_vector=T_vector(30:33,:);% Neumann 4����  �ܲ� %21
Krogager_vector=T_vector1(34:36,:);% Krogager 3���� %87 ���
% TSVM_vector=T_vector(37:52,:);% TSVM 16����   �ܲ� %28
% H_alpha_vector=T_vector(53:93,:);% H/alpha 41����  �ܲ� %28
% T_vector=orgin_vector;
T_vector1=[Freeman_vector;Vanzyl_vector;Yamaguchi_vector;Krogager_vector];
% T_vector1=[orgin_vector;Freeman_vector;Vanzyl_vector;Yamaguchi_vector;Krogager_vector];
T_vector1=T_vector1';






[sample_number,deep] = size(T_vector1);
PolSAR_feature_mat1=reshape(T_vector1,rows,cols,deep );

M=9;
t=floor(M/2);
PolSAR_feature_mat_sp1=padarray(PolSAR_feature_mat1,[t t],'symmetric','both');
PolSAR_feature_tensor1=ones(M,M,deep,sample_number);

num=1;
%% **********************************************************************************************
%^^^^^^^^Ϊ����reshape ����ƥ�䣬����rows*colsתΪsample_numberʱ��Ҫ���а������У�һ��һ�е�������������ʾ^^^^^^^^^
%% **********************************************************************************************

tic
for j=1+t:cols+t
    for i=1+t:rows+t   
        block1=PolSAR_feature_mat_sp1(i-t:i+t,j-t:j+t,:);
        PolSAR_feature_tensor1(:,:,:,num)=block1; 
        num=num+1;
    end
    disp(['��',num2str(j-t),'��'])
end


save('PolSAR_feature_tensor1.mat','PolSAR_feature_tensor1','-v7.3');
save('T_vector1.mat');
toc

T_vector2=PolSAR_feature_rectangle';%�����T_vector������涨Ϊÿһ��Ϊһ��������ÿһ�д���һά������
orgin_vector=T_vector2(1:9,:);%ԭʼ9ά���� %85
% orignDB_vector=T_vector(10:19,:);%ԭʼ9ά����ȡģȡ��λ ������ܲ� %21
Freeman_vector=T_vector2(20:22,:);%Freeman 3���� %85.4
Vanzyl_vector=T_vector2(23:25,:);%Vanzyl 3���� %82.29
Yamaguchi_vector=T_vector2(26:29,:);% Yamaguchi 4���� %82.91
% T_vector=T_vector(30:33,:);% Neumann 4����  �ܲ� %21
Krogager_vector=T_vector2(34:36,:);% Krogager 3���� %87 ���
% TSVM_vector=T_vector(37:52,:);% TSVM 16����   �ܲ� %28
% H_alpha_vector=T_vector(53:93,:);% H/alpha 41����  �ܲ� %28
% T_vector=orgin_vector;
% T_vector=[Freeman_vector;Vanzyl_vector;Yamaguchi_vector;Krogager_vector];
T_vector2=[orgin_vector;Freeman_vector;Yamaguchi_vector;Krogager_vector];
T_vector2=T_vector2';






[sample_number,deep] = size(T_vector2);
PolSAR_feature_mat2=reshape(T_vector2,rows,cols,deep );

M=9;
t=floor(M/2);
PolSAR_feature_mat_sp2=padarray(PolSAR_feature_mat2,[t t],'symmetric','both');
PolSAR_feature_tensor2=ones(M,M,deep,sample_number);

num=1;
%% **********************************************************************************************
%^^^^^^^^Ϊ����reshape ����ƥ�䣬����rows*colsתΪsample_numberʱ��Ҫ���а������У�һ��һ�е�������������ʾ^^^^^^^^^
%% **********************************************************************************************

tic
for j=1+t:cols+t
    for i=1+t:rows+t   
        block2=PolSAR_feature_mat_sp2(i-t:i+t,j-t:j+t,:);
        PolSAR_feature_tensor2(:,:,:,num)=block2; 
        num=num+1;
    end
    disp(['��',num2str(j-t),'��'])
end


save('PolSAR_feature_tensor2.mat','PolSAR_feature_tensor2','-v7.3');
save('T_vector2.mat');
toc

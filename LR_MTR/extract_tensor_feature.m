clear;
clc;

load PolSAR_feature_rectangle.mat;
disp('load data done.');
groundtruth=imread('Got_label_picture.bmp');
[rows, cols, ~] = size(groundtruth);


T_vector=PolSAR_feature_rectangle';%这里的T_vector在这里规定为每一列为一个样本，每一行代表一维的特征
orgin_vector=T_vector(1:9,:);%原始9维数据 %85
% orignDB_vector=T_vector(10:19,:);%原始9维数据取模取相位 用这个很差 %21
Freeman_vector=T_vector(20:22,:);%Freeman 3分量 %85.4
Vanzyl_vector=T_vector(23:25,:);%Vanzyl 3分量 %82.29
Yamaguchi_vector=T_vector(26:29,:);% Yamaguchi 4分量 %82.91
% T_vector=T_vector(30:33,:);% Neumann 4分量  很差 %21
Krogager_vector=T_vector(34:36,:);% Krogager 3分量 %87 最好
% TSVM_vector=T_vector(37:52,:);% TSVM 16分量   很差 %28
% H_alpha_vector=T_vector(53:93,:);% H/alpha 41分量  很差 %28
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
%^^^^^^^^为了与reshape 进行匹配，当将rows*cols转为sample_number时需要进行案列排列，一列一列叠起来，如下所示^^^^^^^^^
%% **********************************************************************************************

tic
for j=1+t:cols+t
    for i=1+t:rows+t   
        block=PolSAR_feature_mat_sp(i-t:i+t,j-t:j+t,:);
        PolSAR_feature_tensor(:,:,:,num)=block; 
        num=num+1;
    end
    disp(['第',num2str(j-t),'列'])
end


save('PolSAR_feature_tensor.mat','PolSAR_feature_tensor','-v7.3');
save('T_vector.mat');
toc










T_vector1=PolSAR_feature_rectangle';%这里的T_vector在这里规定为每一列为一个样本，每一行代表一维的特征
orgin_vector=T_vector1(1:9,:);%原始9维数据 %85
% orignDB_vector=T_vector(10:19,:);%原始9维数据取模取相位 用这个很差 %21
Freeman_vector=T_vector1(20:22,:);%Freeman 3分量 %85.4
Vanzyl_vector=T_vector1(23:25,:);%Vanzyl 3分量 %82.29
Yamaguchi_vector=T_vector1(26:29,:);% Yamaguchi 4分量 %82.91
% T_vector=T_vector(30:33,:);% Neumann 4分量  很差 %21
Krogager_vector=T_vector1(34:36,:);% Krogager 3分量 %87 最好
% TSVM_vector=T_vector(37:52,:);% TSVM 16分量   很差 %28
% H_alpha_vector=T_vector(53:93,:);% H/alpha 41分量  很差 %28
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
%^^^^^^^^为了与reshape 进行匹配，当将rows*cols转为sample_number时需要进行案列排列，一列一列叠起来，如下所示^^^^^^^^^
%% **********************************************************************************************

tic
for j=1+t:cols+t
    for i=1+t:rows+t   
        block1=PolSAR_feature_mat_sp1(i-t:i+t,j-t:j+t,:);
        PolSAR_feature_tensor1(:,:,:,num)=block1; 
        num=num+1;
    end
    disp(['第',num2str(j-t),'列'])
end


save('PolSAR_feature_tensor1.mat','PolSAR_feature_tensor1','-v7.3');
save('T_vector1.mat');
toc

T_vector2=PolSAR_feature_rectangle';%这里的T_vector在这里规定为每一列为一个样本，每一行代表一维的特征
orgin_vector=T_vector2(1:9,:);%原始9维数据 %85
% orignDB_vector=T_vector(10:19,:);%原始9维数据取模取相位 用这个很差 %21
Freeman_vector=T_vector2(20:22,:);%Freeman 3分量 %85.4
Vanzyl_vector=T_vector2(23:25,:);%Vanzyl 3分量 %82.29
Yamaguchi_vector=T_vector2(26:29,:);% Yamaguchi 4分量 %82.91
% T_vector=T_vector(30:33,:);% Neumann 4分量  很差 %21
Krogager_vector=T_vector2(34:36,:);% Krogager 3分量 %87 最好
% TSVM_vector=T_vector(37:52,:);% TSVM 16分量   很差 %28
% H_alpha_vector=T_vector(53:93,:);% H/alpha 41分量  很差 %28
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
%^^^^^^^^为了与reshape 进行匹配，当将rows*cols转为sample_number时需要进行案列排列，一列一列叠起来，如下所示^^^^^^^^^
%% **********************************************************************************************

tic
for j=1+t:cols+t
    for i=1+t:rows+t   
        block2=PolSAR_feature_mat_sp2(i-t:i+t,j-t:j+t,:);
        PolSAR_feature_tensor2(:,:,:,num)=block2; 
        num=num+1;
    end
    disp(['第',num2str(j-t),'列'])
end


save('PolSAR_feature_tensor2.mat','PolSAR_feature_tensor2','-v7.3');
save('T_vector2.mat');
toc

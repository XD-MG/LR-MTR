clear all;clc;
close all;
%% Ëã·¨Ö´ÐÐ

% load('ShangHai_LBP.mat');
% load('ShangHai_GLCM.mat');
load('Tokyo_GLCM.mat');
% load('Tokyo_LBP.mat');
[acc,nmi,f,ri] = lr_mtr(X, gt, lambda, lambda2, X_test, gt_train);
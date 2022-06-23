function [acc,nmi,f,ri] = lr_mtr(X, gt, lambda, lambda2, X_test, gt_train)

% Tokyo data set
NUM(1) = 855;   NUM(2) = 56; NUM(3) = 697; NUM(4) = 205; NUM(5) = 187;
% Shanghai data set
%         NUM(1) = 262;   NUM(2) = 346; NUM(3) = 378; NUM(4) = 417; NUM(5) = 197;

V = length(X);
cls_num = length(unique(gt));
train_num=50;%每一类选二百个做训练样本
%% Note: each column is an sample (same as in LRR)
%%
K = length(X); N = size(X{1},2); %sample number
ModCount = 5; %unfold_mode_number
big_loop=1;
for FT =15% 降维维度
    FT
    % FFT=(1:2:15);% %降维降到的维度
    % NN=size(FFT,2);% 降维需要循环次数
    XX=[X{1};X{2};X{3}];
    Xtest=[X_test{1};X_test{2};X_test{3}];
    dim1 = N;dim2 = N;dim3 = K;
    myNorm = 'tSVD_1';
    sX = [N, N, K];
    Isconverg = 0;
    epson = 1e-5;
    ModCount = 3; %unfold_mode_number
    for v=1:ModCount
        %     para_ten{v} = lambda;
        para_ten{v} = 0.1;
    end
    %  for  ii=1:1:NN
    %      FT=FFT(ii);
    [UUU01,DX,DV]=svd(XX,'econ');
    UUU0=rand(size(UUU01(:,1:FT)));
    EE = zeros(FT,N);
    HH = zeros(N,N);
    YY1 = zeros(FT,N);
    for k=1:K
        Z{k} = zeros(N,N); %Z{2} = zeros(N,N);
        W{k} = zeros(N,N);
        G{k} = zeros(N,N);
        SS{k} = zeros(N,N);
        RR{k} = zeros(N,N);
        E{k} = zeros(size(X{k},1),N); %E{2} = zeros(size(X{k},1),N);
        Y{k} = zeros(size(X{k},1),N); %Y{2} = zeros(size(X{k},1),N);
    end
    w = zeros(N*N*K,1);
    g = zeros(N*N*K,1);
    mu = 10e-10; max_mu = 10e1000; pho_mu = 2;
    rho = 10e-10; max_rho = 10e1000; pho_rho = 2;
    Z_tensor = cat(3, Z{:,:});
    G_tensor = cat(3, G{:,:});
    W_tensor = cat(3, W{:,:});
    for i=1:ModCount
        WT{i} = W_tensor;
    end
    
    for lambda_loop=2
        %     for lambda_loop=2
        lambda = lambda_loop/10;
        plot_x(lambda_loop) = lambda;
        iter_out=1
        close all;
        if iter_out==1
            UUU=UUU0;
        end
        iter_out;
        if iter_out >1
            YY1= zeros(FT,N);
            EE = zeros(FT,N);
        end
        for iter=1:300
            plot_iter(iter) = iter;
            if mod(iter,10) == 0 && iter>=100
                flag = 1;
            end
            for k=1:K
                tmp = (X{k}'*Y{k} + mu*X{k}'*X{k} - mu*X{k}'*E{k} + mu*HH - mu*SS{k} - RR{k}- W{k})./rho +  G{k};
                Z{k}=inv(2*eye(N,N)+ (mu/rho)*X{k}'*X{k})*tmp;
                F = [];
                for v=1:K
                    F = [F;X{v}-X{v}*Z{v}+Y{v}/mu];
                end
                [Econcat] = solve_l1l2(F,lambda/mu);
                beg_ind = 0;
                end_ind = 0;
                for v=1:K
                    if(v>1)
                        beg_ind = beg_ind+size(X{v-1},1);
                    else
                        beg_ind = 1;
                    end
                    end_ind = end_ind+size(X{v},1);
                    E{v} = Econcat(beg_ind:end_ind,:);
                end
                Y{k} = Y{k} + mu*(X{k}-X{k}*Z{k}-E{k});
                tmpp = (mu*XX'*UUU *UUU' *XX + XX'*UUU *YY1 - mu*XX'*UUU *EE + mu*SS{k} + RR{k}+ mu*Z{k} )./rho;
                HH=inv(eye(N,N)+ (mu/rho)*XX'*UUU *UUU' *XX)*tmpp;
                RR{k} = RR{k} + mu*(SS{k}- HH +Z{k});
                stmp = ( mu*HH - mu*Z{k} - RR{k})./rho;
                SS{k}=inv(2*eye(N,N))*stmp;
            end
            stmmp = ( mu*UUU' *XX - mu*UUU' *XX *HH + YY1)./rho;
            EE=inv(2*eye(FT,FT))*stmmp;
            YY1 = YY1 + mu*(UUU' *XX-UUU' *XX *HH-EE);
            Z_tensor = cat(3, Z{:,:});
            W_tensor = cat(3, W{:,:});
            z = Z_tensor(:);
            w = W_tensor(:);
            for umod=1:ModCount
                G_tensor = updateG_tensor(WT{umod},Z,sX,mu,para_ten,V,umod);
                WT{umod} = WT{umod}+mu*(Z_tensor-G_tensor);
            end
            Isconverg = 1;
            for k=1:K
                if (norm(X{k}-X{k}*Z{k}-E{k},inf)>epson)
                    history.norm_Z = norm(X{k}-X{k}*Z{k}-E{k},inf);
                    %                     fprintf('    norm_Z %7.10f    ', history.norm_Z);
                    Isconverg = 0;
                end
                G{k} = G_tensor(:,:,k);
                W_tensor = cat(3, WT{:,:});
                W{k} = W_tensor(:,:,k);
                if (norm(Z{k}-G{k},inf)>epson)
                    history.norm_Z_G = norm(Z{k}-G{k},inf);
                    %                     fprintf('norm_Z_G %7.10f    \n', history.norm_Z_G);
                    Isconverg = 0;
                end
                if (norm(SS{k}-HH+Z{k},inf)>epson)
                    history.norm_Z_S = norm(SS{k}-HH+Z{k},inf);
                    Isconverg = 0;
                end
            end
            BB{iter}= Z{1};
            if iter>1
                BBB=norm( BB{iter},inf);
            end
            BB1{iter}= G{1};
            if iter>1
                BBB1=norm( BB1{iter},inf);
            end
            BB2{iter}= SS{1};
            if iter>1
                BBB2=norm( BB2{iter},inf);
            end
            BB3{iter}= HH;
            if iter>1
                BBB3=norm( BB3{iter},inf);
            end
            Error1(iter) = history.norm_Z;
            Error2(iter) = history.norm_Z_G;
            Error3(iter) = history.norm_Z_S;

            if (norm(UUU' *XX-UUU' *XX *HH-EE,inf)>epson)
                history.norm_Z_S_s = norm(UUU' *XX-UUU' *XX *HH-EE,inf);
                Error4(iter) = history.norm_Z_S_s;
                Isconverg = 0;
            end

            if Isconverg == 1
                Isconverg
            else
                mu = min(mu*pho_mu, max_mu);
                rho = min(rho*pho_rho, max_rho);
                %     iter = iter+1;
            end
        end
        ttemp=UUU'*(XX-XX *HH);
        DIAG=diag(0.5./sqrt(sum(ttemp.*ttemp))+eps);
        SO=(XX-XX *HH)*DIAG*(XX-XX *HH)';
        SO=(SO+SO')/2;
        [Pall,DS]=eig(SO);
        [ds,indd]=sort(diag(DS),'ascend');
        Pall=Pall(:,indd);
        [ind2]=find(ds>10^-10);
        clear UUU;
        ind2=ind2(1:FT,:);
        UUU=Pall(:,ind2);
        
        % 东京数据
        eigvector=fliplr(UUU);%max to min
        projection=UUU;
        YYtrain=projection'*XX;
        YYtest=projection'*Xtest;
        Xtrain_label=gt';
        Xtest_label=gt_train';
        
        %% KNN
        [rateYale_LRPPmaxd1_1(big_loop), rate1(big_loop), rate2(big_loop), rate3(big_loop), rate4(big_loop),rate5(big_loop),testlabel] = KNN_Classfier(YYtrain, Xtrain_label, YYtest,Xtest_label, 5, NUM);%Kneighbor,stop,
        plot_y_KNN(lambda_loop) = rateYale_LRPPmaxd1_1;
        KK=5;
        actual_label_constract = Xtest_label;
        result = testlabel;
        Kappa_KNN = Kappa_x(KK,actual_label_constract,result);
        [A nmi avgent] = compute_nmi(Xtest_label,testlabel);
        [f,p,r] = compute_f(Xtest_label,testlabel);
        [ar,ri,MI,HI]=RandIndex(Xtest_label,testlabel);
        [confmatrix] = cfmatrix(Xtest_label,testlabel);
    end
    temp_18((FT+1)/2) = rateYale_LRPPmaxd1_1(big_loop);
end
%% SVM
[~,bestc,bestg] = SVMcg(gt,YYtrain');
cmd = ['-c ',num2str(bestc),' -g ',num2str(bestg)];
model = svmtrain(gt,YYtrain',cmd);
[predict_label,acc,~] = svmpredict(Xtest_label',YYtest', model);
[rate11(big_loop), rate22(big_loop), rate33(big_loop), rate44(big_loop), rate55(big_loop)] = RATE_Classfier(predict_label', YYtest, Xtest_label, NUM);


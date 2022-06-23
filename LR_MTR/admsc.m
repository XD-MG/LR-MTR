function [acc,nmi,f,ri] = lt_msc(X, gt, lambda)
V = length(X);
cls_num = length(unique(gt));
%% Note: each column is an sample (same as in LRR)
%% 
K = length(X); N = size(X{1},2); %sample number

EE = zeros(size(X{1},1),N);
HH = zeros(N,N);
UU{1} =  zeros(size(X{1},1),size(X{1},1));
UU{2} =  zeros(size(X{1},1),size(X{1},1));
UU{3} =  zeros(size(X{1},1),size(X{1},1));
UUU=[UU{1};UU{2};UU{3}];
YY1 = zeros(size(X{1},1),N);
XX=[X{1};X{2};X{3}];
DD = zeros(size(X{1},1),N);
CC = zeros(size(X{1},1),size(X{1},1));
II = zeros(size(X{1},1),size(X{1},1));


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
dim1 = N;dim2 = N;dim3 = K;
myNorm = 'tSVD_1';
sX = [N, N, K];
%set Default
parOP         =    false;
ABSTOL        =    1e-6;%ø…ƒ‹ «¶Ã
RELTOL        =    1e-4;


Isconverg = 0;epson = 1e-7;
ModCount = 3; %unfold_mode_number
for v=1:ModCount
    para_ten{v} = lambda;
end
iter = 0;
mu = 10e-5; max_mu = 10e10; pho_mu = 2;
rho = 10e-5; max_rho = 10e12; pho_rho = 2;
tic;

Z_tensor = cat(3, Z{:,:});
G_tensor = cat(3, G{:,:});
W_tensor = cat(3, W{:,:});

for i=1:ModCount
    WT{i} = W_tensor;
end

maxIter = 10;
[UUU01,DX,DV]=svd(XX,'econ');
        UUU0=rand(size(UUU01(:,1:8)));

    while iter<maxIter
    iter = iter + 1;
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
        RR{k} = RR{k} + mu*(SS{k}- HH +Z{k});
         FFF = [];
        for v=1:K
             FFF = [FFF;HH+RR{v}-Z{v}/mu];
        end

        [Sconcat] = solve_l1l2(FFF,lambda/mu);
        
        bbeg_ind = 0;
        eend_ind = 0;
        for v=1:K
            if(v>1)
                bbeg_ind = bbeg_ind+size(RR{v-1},1);
            else
                bbeg_ind = 1;
            end
            eend_ind = eend_ind+size(RR{v},1);
            SS{v} = Sconcat(bbeg_ind:eend_ind,:);
        end
        tmpp = (mu*XX'*UUU *UUU' *XX + XX'*UUU *YY1 - mu*XX'*UUU *EE + mu*SS{k} + RR{k}+ Z{k} )./rho;
        HH=inv(eye(N,N)+ (mu/rho)*XX'*UUU *UUU' *XX)*tmpp;
        F11 = [];
        F11 = [F11;UUU' *XX-UUU' *XX *HH+YY1/mu];
        [E11concat] = solve_l1l2(F11,lambda/mu);
        EE = E11concat;
               
        YY1 = YY1 + mu*(UUU' *XX-UUU' *XX *HH-EE);
    end
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
           % fprintf('    norm_Z %7.10f    ', history.norm_Z);
            Isconverg = 0;
        end
        
        G{k} = G_tensor(:,:,k);
        W_tensor = cat(3, WT{:,:});
        W{k} = W_tensor(:,:,k);
        if (norm(Z{k}-G{k},inf)>epson)
            history.norm_Z_G = norm(Z{k}-G{k},inf);
            %fprintf('norm_Z_G %7.10f    \n', history.norm_Z_G);
            Isconverg = 0;
        end
  
        if (norm(SS{k}-HH+Z{k},inf)>epson)
            history.norm_Z_S = norm(SS{k}-HH+Z{k},inf);
           % fprintf('    norm_Z %7.10f    ', history.norm_Z);
            Isconverg = 0;
        end
        
    end
   
     if (norm(UUU' *XX-UUU' *XX *HH-EE,inf)>epson)
            history.norm_Z_S_s = norm(UUU' *XX-UUU' *XX *HH-EE,inf);
           % fprintf('    norm_Z %7.10f    ', history.norm_Z);
            Isconverg = 0;
     end
      
    if Isconverg == 0
        break;
    else
    mu = min(mu*pho_mu, max_mu);
    rho = min(rho*pho_rho, max_rho); 
    end
 end


    







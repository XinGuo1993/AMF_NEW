%% Demo Codes for YaleB face 
%%
clear;
clc;

warning('off')
addpath('code');
addpath('code/SkewNorm');
addpath('code/AsymLaplace_modify');
addpath('code/AMF');
addpath('code/BRMF');
addpath('code/BRMF/Utilities');
addpath('code/BRMF/mex');
addpath('code/PCP');
addpath('code/Mog-RPCA');
addpath('code/CWM');
addpath('code/VBLR');
addpath('code/BRPCA');
addpath('code/PRMF');

%% Prepare the data

m = 100;
n = 100;                     % Data size
R = 4;                       % Ture Rank
c =1;                        % 实验次数
 
for i = 1:c
%% normalization  
    %X_Ori = normalize_std(X_Ori);
    RU = randn(m,R);
    RV = randn(R,n);
    X_Ori = RU * RV;
    W =  rand(size(X_Ori))>0.20;
%   W = ones(size(X_Ori));
   
    X_Noi = X_Ori;
    Ind = randperm(m*n);
    p1 = floor(m*n*0.15);
    p2 = floor(m*n*0.2);
    X_Noi(Ind(1:p1))    =   X_Noi(Ind(1:p1))+randn(1,p1)* 0.25;  %Add Gaussian noise
    X_Noi(Ind(p1+1:p1+p2)) = X_Noi(Ind(p1+1:p1+p2)) + (rand(1,p2)*10) -5; %Add uniform noise
    p3 = floor(m*n*0.2);
    X_Noi(Ind(p1+p2+1:p1+p2+p3)) = X_Noi(Ind(p1+p2+1:p1+p2+p3)) + (rand(1,p3)*4) -2; %Add uniform noise
% %     X_Noi(Ind(p1+p2+p3+1:end)) = X_Noi(Ind(p1+p2+p3+1:end)) + randn(1,m*n-p1-p2-p3)* 0.1; 
% %   X_Noi = W.*X_Noi;       %Add missing components
%     X_Noi = X_Noi + ones(m,n)*0.05 + randn(m,n)* 0.25;


%%
%%%%%%%%%%%%%%%%%%%%%%%% PCP (Robust PCA) %%%%%%%%%%%%%%%%%%%%%%%%
r = R*2;

tic;
[B,E] = PCP(X_Noi, 1/sqrt(max(m,n)), 1e-3, r);
time_PCP = toc;

E1_PCP(i) = norm((X_Ori - B),'fro')/norm(X_Ori,'fro');

%%
%%%%%%%%%%%%%%%%%%%%%%%%  BRMF %%%%%%%%%%%%%%%%%%%%%%%%  
r = R*2;
opts.maxIter = 100;
opts.burnin = 50;
opts.invW_0 = 1000 * eye(r);
opts.beta_0 = 2;
opts.nu_0 = r;
opts.a = 1e-4;
opts.b = 1e0;
% We set the maximum rank to be twice of the ground truth.
opts.r = r;
opts.alpha = 0.5;

tic;
[U, V, Tau, p] = BRMF(X_Noi, opts);  %p是U*V
time_BRMF = toc;

E1_BRMP(i) = norm((X_Ori - p),'fro')/norm(X_Ori,'fro');

%% 
%%%%%%%%%%%%%%%%%%%%%%%%  Markov BRMF  %%%%%%%%%%%%%%%%%%%%%%%%
 % We use the results of BRMF as initialization. We may also use random
 % initialization.
opts.maxIter = 200;
opts.burnin = 100;
tic;
[U, V, Tau, p] = MBRMF(X_Noi, opts, U, V);
time_MBRMF = toc;

E1_MBRMP(i) = norm((X_Ori - p),'fro')/norm(X_Ori,'fro');


%%
%%%%%%%%%%%%%%%%%%%%%%%% MRPCA/MoG %%%%%%%%%%%%%%%%%%%%%%%%
r= 2*R;
param.mog_k = 3;
param.lr_init = 'SVD';
param.maxiter = 100;
param.initial_rank = r;
param.tol = 1e-3;
    
lr_prior.a0 = 1e-6;
lr_prior.b0 = 1e-6;

mog_prior.mu0 = 0;
mog_prior.c0 = 1e-3;
mog_prior.d0 = 1e-3;
mog_prior.alpha0 = 1e-3;
mog_prior.beta0 = 1e-3;

tic;
[lr_model, mog_model, r] = mog_rpca(X_Noi, param, lr_prior, mog_prior);
time_MRPCA = toc;

L = lr_model.U*lr_model.V';
E1_MRPCA(i) = norm(L-X_Ori,'fro')/norm(X_Ori,'fro');

%%
%%%%%%%%%%%%%%%%%%%%%%%% CWM %%%%%%%%%%%%%%%%%%%%%%%%

r = R*2;

% Initialize U and V
    MedValX = median(abs(X_Noi(:)));
    MedValX = sqrt(MedValX/r);
    param.IniU = rand(m,r)*MedValX*2-MedValX;
    param.IniV = rand(n,r)*MedValX*2-MedValX;
    tol = 1e-10;
% CWM implementation
    tic;
    [U,V] = CWM(X_Noi,r,param);
    time_CWM = toc;
    
    E1_CWM(i) = norm(U*V'-X_Ori,'fro')/norm(X_Ori,'fro');

%%
%%%%%%%%%%%%%%%%%%%%%%%% VBLR %%%%%%%%%%%%%%%%%%%%%%%%
    
options.verbose = 1;
options.initial_rank = 'auto'; % This sets to the maximum possible rank
options.inf_flag = 2; % inference flag for the sparse component
options.MAXITER = 200;
options.UPDATE_BETA = 1; 
options.mode = 'VB';

tic
[X_hat, A_hat, B_hat, E_hat] = VBRPCA(X_Noi,options);
time_VBLR = toc;

E1_VBLR(i) = norm( X_hat - X_Ori, 'fro' ) / norm( X_Ori, 'fro');    
    
%%
%%%%%%%%%%%%%%%%%%%%%%%% BRPCA %%%%%%%%%%%%%%%%%%%%%%%%
r = R*2;
Theta0 = InitialPara_random(X_Noi,r);
tic;
Output = Bayesian_RPCAmcmc(X_Noi,Theta0);
time_BRPCA = toc;
E1_BRPCA(i) = norm(X_Ori-Output.Lowrank_mean,'fro')/norm(X_Ori,'fro'); 

%%
%%%%%%%%%%%%%%%%%%%%%%%% PRMF %%%%%%%%%%%%%%%%%%%%%%%%
r = R*2;
X_Noi_N = normalize(X_Noi);
tic;
[P, Q] = RPMF(X_Noi, r, 1, 1, 1e-2);
time_PRMF = toc;
E1_PRMF(i) = norm(X_Ori- P * Q,'fro')/norm(X_Ori,'fro');

%%
%%%%%%%%%%%%%%%%%%%%%%%% AMF method %%%%%%%%%%%%%%%%%%%%%%%%
K=64;   %initialize cluster number
r=R*2; %initial rank
Prune =0; %decide whether prune Rank

        opt.maxIter =300;
        opt.r  =r;
        opt.a0 =1e-4;
        opt.b0 =1e-4;
        opt.a1 =1e-1;
        opt.b1 =1e-1;
        opt.alpha0 = 1;
        opt.K      =K;
        opt.tol = 1e-5;
        opt.Prune  =Prune;
        opt.verbose =1;
        opt.minIter =50;
tic;
[U V X_sb Phi_sb opt] = DemoVarInference(X_Noi, W, opt);
time_AMF = toc;
E1_AMF(i) = norm((X_Ori - X_sb),'fro')/norm(X_Ori,'fro');

%%
%%%%%%%%%%%%%%%%%%%%%%%% SkewNormAMF method %%%%%%%%%%%%%%%%%%%%%%%%
K=64;   %initialize cluster number
r=R*2; %initial rank
Prune =1; %decide whether prune Rank
        opt.lamda_s = 0; % 偏移量
        opt.maxIter =300;
        opt.r  =r;
        opt.a0 = 1e-4;
        opt.b0 =1e-4;
        opt.a1 =1e-1;  
        opt.b1 =1e-1;
        opt.alpha0 = 1;
        opt.K      =K;
        opt.tol = 1e-5;
        opt.Prune  =Prune;
        opt.verbose =1;
        opt.minIter =50;
        tic;
        [A_sb B_sb X_sb Phi_sb opt] = SkewNormVarInference(X_Noi, W, opt);
        time_SkewNormAMF = toc;

 E1_SkewAMF(i) = norm((X_Ori - X_sb),'fro')/norm(X_Ori,'fro');       

%%
%%%%%%%%%%%%%%%%%%%%%%%% AsymLaplaceAMF method %%%%%%%%%%%%%%%%%%%%%%%%
K=64;   %initialize cluster number
r=R*2; %initial rank
Prune =1; %decide whether prune Rank

        opt.tao_s = 0.50; % 偏移量
        
        opt.s_alpha  = 1;
        opt.s_beta = 0.5;
        opt.maxIter =300;
        opt.r  =r;
        opt.a0 =1e-4;
        opt.b0 =1e-4;
        opt.a1 =1e-1;
        opt.b1 =1e-1;
        opt.alpha0 = 1;
        opt.K      =K;
        opt.tol = 1e-5;
        opt.Prune  =Prune;
        opt.verbose =1;
        opt.minIter =50;
        tic;
        [A_sb B_sb X_sb Phi_sb opt] = AsymLaplaceVarInference(X_Noi, W, opt);
        time_AsymLaplaceAMF = toc;

E1_AsymAMF(i) = norm((X_Ori - X_sb),'fro')/norm(X_Ori,'fro');
end

   
E2_PCP = mean(E1_PCP);
std_PCP =  std(E1_PCP);

E2_BRMP = mean(E1_BRMP);
std_BRMP =  std(E1_BRMP);

E2_MBRMP = mean(E1_MBRMP);
std_MBRMP =  std(E1_MBRMP);

E2_MRPCA = mean(E1_MRPCA);
std_MRPCA = std(E1_MRPCA);

E2_CWM = mean(E1_CWM);
std_CWM = std(E1_CWM);

E2_VBLR = mean(E1_VBLR);
std_VBLR = std(E1_VBLR);

E2_BRPCA = mean(E1_BRPCA);
std_BRPCA = std(E1_BRPCA);

E2_PRMF = mean(E1_PRMF);
std_PRMF = std(E1_PRMF);

E2_AMF = mean(E1_AMF);
std_AMF =  std(E1_AMF);

E2_SkewAMF = mean(E1_SkewAMF);
std_SkewAMF =  std(E1_SkewAMF);

E2_AsymAMF = mean(E1_AsymAMF);
std_AsymAMF =  std(E1_AsymAMF);

matlabpool close;

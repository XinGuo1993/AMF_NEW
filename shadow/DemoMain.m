%% Demo Codes for test removal 
%%
clear;
clc;

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

% matlabpool open;
%% Prepare the data
R = 10;                       % Ture Rank

pic = imread('160.bmp');
B0 = (double(pic(:, :, 1)) - 128) / 128;  %无文字矩阵  
pic = imread('161.bmp');
D = (double(pic(:, :, 1)) - 128) / 128;   %有文字矩阵Y
[m, n] = size(B0);

S0 = (D == (41 - 128) / 128);  %是文字just a simple test,41是文字的值，S0为0，1矩阵,1在matlab中是白色
outNum = sum(sum(S0)); %有文字的像素和
D(S0) = rand(outNum, 1) * 2 - 0.5; %对原图有文字的地方加噪声,值域【0，1】  噪声与文中不符和
W = ones(m,n); %AMF等算法的mask

%%
%%%%%%%%%%%%%%%%%%%%%%%% PCP (Robust PCA) %%%%%%%%%%%%%%%%%%%%%%%%
r = 2*R;

tic;
[B,E] = PCP(D, 1/sqrt(max(m,n)), 1e-3, r);
time_PCP = toc;
% We choose the threshold to maximize the F measure.
[f S] = findFMeasure(abs(E), S0);

figure;
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (PCP), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(B, [0 1]),colormap gray; axis off; title(['Recovery (PCP), error = ', num2str(norm((B0 - B), 'fro')/norm(B,'fro'))]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%  BRMF  %%%%%%%%%%%%%%%%%%%%%%%%
r = R*3;
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
[U, V, Tau, p] = BRMF(D, opts);  %p是U*V
time_BRMF = toc;

[f S] = findFMeasure(abs(D - p), S0); %f是f测度，S是mask

figure;
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (BRMF), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(p, [0 1]),colormap gray; axis off; title(['Recovery (BRMF), error = ', num2str(norm((B0 - p), 'fro')/norm(p,'fro'))]);

%% Markov BRMF
 % We use the results of BRMF as initialization. We may also use random
 % initialization.
opts.maxIter = 200;
opts.burnin = 100;
tic;
[U, V, Tau, p] = MBRMF(D, opts, U, V);
time_MBRMF = toc;

[f S] = findFMeasure(abs(D - p) .* Tau , S0);

figure;
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (MBRMF), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(p, [0 1]),colormap gray; axis off; title(['Recovery (MBRMF), error = ', num2str(norm((B0 - p), 'fro')/norm(p,'fro'))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%% MRPCA/MoG %%%%%%%%%%%%%%%%%%%%%%%%

r = R*2;
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
[lr_model, mog_model, r] = mog_rpca(D, param, lr_prior, mog_prior);
time_MRPCA = toc;
L = lr_model.U*lr_model.V'; % L

[f S] = findFMeasure(abs(D - p) , S0);

figure;
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (MRPCA), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(p, [0 1]),colormap gray; axis off; title(['Recovery (MRPCA), error = ', num2str(norm((B0 - p), 'fro')/norm(p,'fro'))]);

%%
%%%%%%%%%%%%%%%%%%%%%%%% CWM %%%%%%%%%%%%%%%%%%%%%%%%

r = R*2;

% Initialize U and V
    MedValX = median(abs(D(:)));
    MedValX = sqrt(MedValX/r);
    param.IniU = rand(m,r)*MedValX*2-MedValX;
    param.IniV = rand(n,r)*MedValX*2-MedValX;
    tol = 1e-10;
% CWM implementation
    tic;
    [U,V] = CWM(D,r,param);
    time_CWM = toc;
    p = U*V';
    
[f S] = findFMeasure(abs(D - p), S0);

figure;
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (CWM), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(p, [0 1]),colormap gray; axis off; title(['Recovery (CWM), error = ', num2str(norm((B0 - p), 'fro')/norm(p,'fro'))]);

%%
%%%%%%%%%%%%%%%%%%%%%%%% VBLR %%%%%%%%%%%%%%%%%%%%%%%%
    
options.verbose = 1;
options.initial_rank = 'auto'; % This sets to the maximum possible rank
options.inf_flag = 2; % inference flag for the sparse component
options.MAXITER = 200;
options.UPDATE_BETA = 1; 
options.mode = 'VB';

tic
[p, A_hat, B_hat, E_hat] = VBRPCA(D,options);
time_VBLR = toc;

[f S] = findFMeasure(abs(D - p) .* Tau , S0);

figure;
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (VBLR), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(p, [0 1]),colormap gray; axis off; title(['Recovery (VBLR), error = ', num2str(norm((B0 - p), 'fro')/norm(p,'fro'))]);  
    
%%
%%%%%%%%%%%%%%%%%%%%%%%% BRPCA %%%%%%%%%%%%%%%%%%%%%%%%

r = R*2;
Theta0 = InitialPara_random(D,r);
tic;
Output = Bayesian_RPCAmcmc(D,Theta0);
time_BRPCA = toc;
p = Output.Lowrank_mean;

[f S] = findFMeasure(abs(D - p), S0);

figure;
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (BRPCA), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(p, [0 1]),colormap gray; axis off; title(['Recovery (BRPCA), error = ', num2str(norm((B0 - p), 'fro')/norm(p,'fro'))]);

%%
%%%%%%%%%%%%%%%%%%%%%%%% PRMF %%%%%%%%%%%%%%%%%%%%%%%%
r = R*2;
tic;
[P, Q] = RPMF(D, r, 1, 1, 1e-2);
time_PRMF = toc;
p = P*Q';

[f S] = findFMeasure(abs(D - p), S0);

figure;
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (PRMF), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(p, [0 1]),colormap gray; axis off; title(['Recovery (PRMF), error = ', num2str(norm((B0 - p), 'fro')/norm(p,'fro'))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        [A_sb B_sb X_sb Phi_sb opt] = DemoVarInference(D, W, opt);
        time_AMF = toc;

        [f S] = findFMeasure(abs(D - X_sb), S0); %f是f测度，S是mask
        
figure;       
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (AMF), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(X_sb, [0 1]),colormap gray; axis off; title(['Recovery (AMF), error = ', num2str(norm((B0 - X_sb), 'fro')/norm(X_sb,'fro'))]);

%%
%%%%%%%%%%%%%%%%%%%%%%%% SkewNormAMF method %%%%%%%%%%%%%%%%%%%%%%%%
K=64;   %initialize cluster number
r=R*2; %initial rank
Prune =0; %decide whether prune Rank
        opt.lamda_s = 0; % 偏移量  
        
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
        [U V X_sb Phi_sb opt] = SkewNormVarInference(D, W, opt);
        time_SkewNormAMF = toc;

        [f S] = findFMeasure(abs(D - X_sb), S0); %f是f测度，S是mask
        
figure;       
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (SkewNormAMF), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(X_sb, [0 1]),colormap gray; axis off; title(['Recovery (SkewNormAMF), error = ', num2str(norm((B0 - X_sb), 'fro')/norm(X_sb,'fro'))]);

%%
%%%%%%%%%%%%%%%%%%%%%%%% AsymLaplaceAMF method %%%%%%%%%%%%%%%%%%%%%%%%
K=64;   %initialize cluster number
r=R*2; %initial rank
Prune =0; %decide whether prune Rank
        opt.tao_s = 0.5; % 偏移量    
        
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
        [U V X_sb Phi_sb opt] = AsymLaplaceVarInference(D, W, opt);
        time_AsymLaplaceAMF = toc;

        [f S] = findFMeasure(abs(D - X_sb), S0); %f是f测度，S是mask
        
figure;       
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (AsymLaplaceAMF), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(X_sb, [0 1]),colormap gray; axis off; title(['Recovery (AsymLaplaceAMF), error = ', num2str(norm((B0 - X_sb), 'fro')/norm(X_sb,'fro'))]);

matlabpool close;



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

%% Prepare the data
R = 8;

gx_m = 48;
gx_n = 42;
% path = 'yaleB01/';
path = '01/';
files = dir(fullfile(path,'*.bmp'));
len = size(files,1); 
face  = zeros(gx_m*gx_n,len);
face2  = zeros(gx_m*gx_n,len);
random_gx = 0;
if random_gx == 1
    idex = randperm(len);
else
    idex = 1:len;
end
for i=1:len
    fileName = strcat(path,files(i,1).name); 
    pic = imread(fileName);
    pic2 = imresize(pic,[gx_m,gx_n],'bicubic');
    face(:,idex(i)) = reshape(pic2,[gx_m*gx_n,1]);
    face2(:,idex(i)) = uint8(Normalize(double(face(:,idex(i)))));
end

Recovery0 = uint8(Normalize(face));
org_pic = reshape(Recovery0,[gx_m,gx_n*len]);
% org_pic2 = reshape(face2,[gx_m,gx_n*len]);
% figure(1)
% imshow(org_pic)
pic = face;
%% normalization  
% D = (double(pic(:, :, 1)) - 128) / 128;
D = double(pic(:, :, 1));
[m,n]=size(pic);   
W = ones(m,n);


%%
%%%%%%%%%%%%%%%%%%%%%%%% PCP (Robust PCA) %%%%%%%%%%%%%%%%%%%%%%%%
r = R;

tic;
[B,E] = PCP(D, 1/sqrt(max(m,n)), 1e-3, r);
time_PCP = toc;
% Recovery = uint8(Normalize(p));
% new_pic = reshape(Recovery,[gx_m,gx_n*len]);
Recovery2 = B;
for i = 1:len
    Recovery2(:,i) = uint8(Normalize(B(:,i)));
end
pic_PCP = uint8(reshape(Recovery2,[gx_m,gx_n*len]));

%%
%%%%%%%%%%%%%%%%%%%%%%%%  BRMF initial %%%%%%%%%%%%%%%%%%%%%%%%  
r = R;
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
% Recovery = uint8(Normalize(p));
% new_pic = reshape(Recovery,[gx_m,gx_n*len]);
Recovery2 = p;
for i = 1:len
    Recovery2(:,i) = uint8(Normalize(p(:,i)));
end
pic_BRMF = uint8(reshape(Recovery2,[gx_m,gx_n*len]));

%% 
%%%%%%%%%%%%%%%%%%%%%%%%  Markov BRMF  %%%%%%%%%%%%%%%%%%%%%%%%
 % We use the results of BRMF as initialization. We may also use random
 % initialization.
opts.maxIter = 200;
opts.burnin = 100;
tic;
[U, V, Tau, p] = MBRMF(D, opts, U, V);
time_MBRMF = toc;

% Recovery = uint8(Normalize(p));
% new_pic = reshape(Recovery,[gx_m,gx_n*len]);
Recovery2 = p;
for i = 1:len
    Recovery2(:,i) = uint8(Normalize(p(:,i)));
end
pic_MBRMF = uint8(reshape(Recovery2,[gx_m,gx_n*len]));

%%
%%%%%%%%%%%%%%%%%%%%%%%% MRPCA/MoG %%%%%%%%%%%%%%%%%%%%%%%%
r = R;
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

p = lr_model.U*lr_model.V';
% Recovery = uint8(Normalize(p));
% new_pic = reshape(Recovery,[gx_m,gx_n*len]);
Recovery2 = p;
for i = 1:len
    Recovery2(:,i) = uint8(Normalize(p(:,i)));
end
pic_MRPCA = uint8(reshape(Recovery2,[gx_m,gx_n*len]));

%%
%%%%%%%%%%%%%%%%%%%%%%%% CWM %%%%%%%%%%%%%%%%%%%%%%%%

r = R;

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
    % Recovery = uint8(Normalize(p));
    % new_pic = reshape(Recovery,[gx_m,gx_n*len]);
    Recovery2 = p;
    for i = 1:len
        Recovery2(:,i) = uint8(Normalize(p(:,i)));
    end
    pic_CWM = uint8(reshape(Recovery2,[gx_m,gx_n*len]));

%%
%%%%%%%%%%%%%%%%%%%%%%%% VBLR %%%%%%%%%%%%%%%%%%%%%%%%
    
options.verbose = 1;
options.initial_rank = 'auto'; % This sets to the maximum possible rank
options.inf_flag = 2; % inference flag for the sparse component
options.MAXITER = 200;
options.UPDATE_BETA = 1; 
options.mode = 'VB';

D1 = (D - 128) / 128;   %不正则化会跪
tic
[p, A_hat, B_hat, E_hat] = VBRPCA(D1,options);
time_VBLR = toc;

% Recovery = uint8(Normalize(p));
% new_pic = reshape(Recovery,[gx_m,gx_n*len]);
Recovery2 = p;
for i = 1:len
    Recovery2(:,i) = uint8(Normalize(p(:,i)));
end
pic_VBLR = uint8(reshape(Recovery2,[gx_m,gx_n*len])); 
    
%%
%%%%%%%%%%%%%%%%%%%%%%%% BRPCA %%%%%%%%%%%%%%%%%%%%%%%%
r = R;
Theta0 = InitialPara_random(D,r);
D1 = (D - 128) / 128;   %不正则化会跪
tic;
Output = Bayesian_RPCAmcmc(D1,Theta0);
time_BRPCA = toc;
p = Output.Lowrank_mean;
% Recovery = uint8(Normalize(p));
% new_pic = reshape(Recovery,[gx_m,gx_n*len]);
Recovery2 = p;
for i = 1:len
    Recovery2(:,i) = uint8(Normalize(p(:,i)));
end
pic_BRPCA = uint8(reshape(Recovery2,[gx_m,gx_n*len]));

%%
%%%%%%%%%%%%%%%%%%%%%%%% PRMF %%%%%%%%%%%%%%%%%%%%%%%%
r = R;
X_Noi_N = normalize(D);
tic;
[P, Q] = RPMF(D, r, 1, 1, 1e-2);
time_PRMF = toc;
p = P*Q;
% Recovery = uint8(Normalize(p));
% new_pic = reshape(Recovery,[gx_m,gx_n*len]);
Recovery2 = p;
for i = 1:len
    Recovery2(:,i) = uint8(Normalize(p(:,i)));
end
pic_PRMF = uint8(reshape(Recovery2,[gx_m,gx_n*len]));


%%
%%%%%%%%%%%%%%%%%%%%%%%% AMF method %%%%%%%%%%%%%%%%%%%%%%%%
K=64;   %initialize cluster number
r = R; %initial rank
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
[U V X_sb Phi_sb opt] = DemoVarInference(D, W, opt);
time_AMF = toc;
% Recovery = uint8(Normalize(X_sb));
% new_pic = reshape(Recovery,[gx_m,gx_n*len]);
Recovery2 = X_sb;
for i = 1:len
    Recovery2(:,i) = uint8(Normalize(X_sb(:,i)));
end
pic_AMF = uint8(reshape(Recovery2,[gx_m,gx_n*len]));

%%
%%%%%%%%%%%%%%%%%%%%%%%% SkewNormAMF method %%%%%%%%%%%%%%%%%%%%%%%%
K=64;   %initialize cluster number
r = R; %initial rank
Prune =0; %decide whether prune Rank
        opt.lamda_s = 0; % 偏移量
        opt.maxIter =300;
        opt.r  =r*2;
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
        [A_sb B_sb X_sb Phi_sb opt] = SkewNormVarInference(D, W, opt);
        time_SkewNormAMF = toc;

% Recovery = uint8(Normalize(X_sb));
% new_pic = reshape(Recovery,[gx_m,gx_n*len]);
Recovery2 = X_sb;
for i = 1:len
    Recovery2(:,i) = uint8(Normalize(X_sb(:,i)));
end
pic_SkewNormAMF = uint8(reshape(Recovery2,[gx_m,gx_n*len]));

%%
%%%%%%%%%%%%%%%%%%%%%%%% AsymLaplaceAMF method %%%%%%%%%%%%%%%%%%%%%%%%
K=64;   %initialize cluster number
r = R; %initial rank
Prune =0; %decide whether prune Rank

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
        [A_sb B_sb X_sb Phi_sb opt] = AsymLaplaceVarInference(D, W, opt);
        time_AsymLaplaceAMF = toc;

% Recovery = uint8(Normalize(X_sb));
% new_pic = reshape(Recovery,[gx_m,gx_n*len]);
Recovery2 = X_sb;
for i = 1:len
    Recovery2(:,i) = uint8(Normalize(X_sb(:,i)));
end
pic_AsymLaplaceAMF = uint8(reshape(Recovery2,[gx_m,gx_n*len]));


guoxin_pic = [org_pic;pic_BRMF;pic_MBRMF;pic_MRPCA;pic_CWM;pic_VBLR;pic_BRPCA;pic_PRMF;pic_AMF;pic_SkewNormAMF;pic_AsymLaplaceAMF];
imwrite(guoxin_pic,'result.bmp');


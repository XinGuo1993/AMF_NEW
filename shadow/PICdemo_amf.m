clear all;
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
                 
    
%% normalization  
pic = imread('160.bmp');
B0 = (double(pic(:, :, 1)) - 128) / 128;
pic = imread('161.bmp');
D = (double(pic(:, :, 1)) - 128) / 128;
[m, n] = size(B0);
S0 = (D == (41 - 128) / 128);
outNum = sum(sum(S0));
D(S0) = 0.5+rand(outNum, 1) * 2 - 0.5;
W = ones(m,n);
     
%%%%%%%%%%%%% AMF method %%%%%%%%%%%%%%%%%%%%%%%%
R = 4;                      % Ture Rank
c =3;
K=10;   %initialize cluster number
r=10; %initial rank
Prune =1; %decide whether prune Rank

        opt.maxIter =300;
        opt.r  =r*2;
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

        [f S] = findFMeasure(abs(D - X_sb), S0); %f «f≤‚∂»£¨S «mask
        
figure;       
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (AMF), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(X_sb, [0 1]),colormap gray; axis off; title(['Recovery (AMF), error = ', num2str(norm((B0 - X_sb), 'fro'))]);


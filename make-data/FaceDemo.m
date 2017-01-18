%face demo for brmf
clear all;
clc;
warning('off')
addpath('Utilities');
K=64;   %initialize cluster number
r=8; %initial rank
Prune =0; %decide whether prune Rank

%% normalization  
pic = imread('yaleB01_P00A+110E-20.pgm');
D = (double(pic(:, :, 1)) - 128) / 128;
D = double(pic(:, :, 1));
[m,n]=size(pic);   
W = ones(m,n);
     
%%%%%%%%%%%%% AsymLaplaceAMF method %%%%%%%%%%%%%%%%%%%%%%%%
r = 10;
opts.maxIter = 100;
opts.burnin = 50;
opts.invW_0 = 1000 * eye(r * 2);
opts.beta_0 = 2;
opts.nu_0 = r * 2;
opts.a = 1e-4;
opts.b = 1e0;
% We set the maximum rank to be twice of the ground truth.
opts.r = r * 2;
opts.alpha = 0.5;

[U, V, Tau, X_sb] = BRMF(D, opts);  %p «U*V
        
figure;       
mask = uint8(Normalize(D - X_sb));
Recovery = uint8(Normalize(X_sb));
subplot(2,2,1),imshow(pic),colormap gray; axis off; title('Input');
subplot(2,2,3),imshow(mask),colormap gray; axis off; title(['Mask (AsymLaplaceAMF)']);
subplot(2,2,4),imshow(Recovery),colormap gray; axis off; title(['Recovery (AsymLaplaceAMF) ']);


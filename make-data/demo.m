%% Demo Codes for "Bayesian Robust Matrix Factorization for Image and Video Processing"
%% If you have any questions, please contact winsty@gmail.com

addpath('PCP');
addpath('Utilities');
addpath('mex');
% matlabpool open;
%% Prepare the data
pic = imread('1.bmp');
% pic = imresize(pic,[160,160],'bicubic');
B0 = (double(pic(:, :, 1)) - 128) / 128;  %无文字矩阵  
pic = imread('2.bmp');
% pic = imresize(pic,[160,160],'bicubic');
D = (double(pic(:, :, 1)) - 128) / 128;   %有文字矩阵Y
[m, n] = size(B0);

S0 = (D == (41 - 128) / 128);  %是文字just a simple test,41是文字的值，S0为0，1矩阵,1在matlab中是白色
% imshow(S0,[0,1])
outNum = sum(sum(S0)); %有文字的像素和
D(S0) = rand(outNum, 1) * 2 - 0.5; %对原图有文字的地方加噪声,值域【0，1】

%% Set up the (hyper)parameters. See more details in the paper.
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

%% PCP (Robust PCA)
tic;
[B,E] = PCP(D, 1/sqrt(max(m,n)), 1e-3, r * 2);
toc;
% We choose the threshold to maximize the F measure.
[f S] = findFMeasure(abs(E), S0);

figure;
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (PCP), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(B, [0 1]),colormap gray; axis off; title(['Recovery (PCP), error = ', num2str(norm((B0 - B), 'fro'))]);

%% BPMF
tic;
[U, V, Tau, p] = BRMF(D, opts);  %p是U*V
toc;

[f S] = findFMeasure(abs(D - p), S0); %f是f测度，S是mask

figure;
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (BRMF), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(p, [0 1]),colormap gray; axis off; title(['Recovery (BRMF), error = ', num2str(norm((B0 - p), 'fro'))]);

% %% Markov BRMF
%  % We use the results of BRMF as initialization. We may also use random
%  % initialization.
% opts.maxIter = 200;
% opts.burnin = 100;
% tic;
% [U, V, Tau, p] = MBRMF(D, opts, U, V);
% toc;
% 
% [f S] = findFMeasure(abs(D - p) .* Tau , S0);
% 
% figure;
% subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
% subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
% subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (MBRMF), F measure = ', num2str(f)]);
% subplot(2,2,4),imshow(p, [0 1]),colormap gray; axis off; title(['Recovery (MBRMF), error = ', num2str(norm((B0 - p), 'fro'))]);

matlabpool close;

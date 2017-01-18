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
B0 = (double(pic(:, :, 1)) - 128) / 128;  %�����־���  
pic = imread('161.bmp');
D = (double(pic(:, :, 1)) - 128) / 128;   %�����־���Y
[m, n] = size(B0);

S0 = (D == (41 - 128) / 128);  %������just a simple test,41�����ֵ�ֵ��S0Ϊ0��1����,1��matlab���ǰ�ɫ
outNum = sum(sum(S0)); %�����ֵ����غ�
D(S0) = rand(outNum, 1) * 2 - 0.5; %��ԭͼ�����ֵĵط�������,ֵ��0��1��  ���������в�����
W = ones(m,n); %AMF���㷨��mask

%%
%%%%%%%%%%%%%%%%%%%%%%%% PRMF %%%%%%%%%%%%%%%%%%%%%%%%
r = R*2;
tic;
[P, Q] = RPMF(D, r, 1, 1, 1e-2);
time_PRMF = toc;
p = P*Q;

[f S] = findFMeasure(abs(D - p), S0);

figure;
subplot(2,2,1),imshow(D, [0 1]),colormap gray; axis off; title('Input');
subplot(2,2,2),imshow(B0, [0 1]),colormap gray; axis off; title('Ground Truth');
subplot(2,2,3),imshow(S, [0 1]),colormap gray; axis off; title(['Mask (PRMF), F measure = ', num2str(f)]);
subplot(2,2,4),imshow(p, [0 1]),colormap gray; axis off; title(['Recovery (PRMF), error = ', num2str(norm((B0 - p), 'fro')/norm(p,'fro'))]);

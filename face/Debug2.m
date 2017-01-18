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
        
        D1 = (D - 128) / 128;   %不正则化会跪
        tic;
        [A_sb B_sb X_sb Phi_sb opt] = SkewNormVarInference(D1, W, opt);
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
        
        D1 = (D - 128) / 128;   %不正则化会跪
        tic;
        [A_sb B_sb X_sb Phi_sb opt] = AsymLaplaceVarInference(D1, W, opt);
        time_AsymLaplaceAMF = toc;

% Recovery = uint8(Normalize(X_sb));
% new_pic = reshape(Recovery,[gx_m,gx_n*len]);
Recovery2 = X_sb;
for i = 1:len
    Recovery2(:,i) = uint8(Normalize(X_sb(:,i)));
end
pic_AsymLaplaceAMF = uint8(reshape(Recovery2,[gx_m,gx_n*len]));


guoxin_pic = [org_pic;pic_PCP;pic_SkewNormAMF;pic_AsymLaplaceAMF];
imwrite(guoxin_pic,'result.bmp');
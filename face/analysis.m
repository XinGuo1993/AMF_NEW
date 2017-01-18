clear;
clc;
load('guoxin.mat');
addpath('code');
m = 48;
n = 42;
num = 7;
Norm = pic_SkewNormAMF(:,(num-1)*n+1:num*n);
Norm_fig = zeros(256,1);
for i =1:m*n
    Norm_fig(Norm(i)+1) = Norm_fig(Norm(i)+1)+1;
end

Asym = pic_AsymLaplaceAMF(:,(num-1)*n+1:num*n);
Asym_fig = zeros(256,1);
for i =1:m*n
    Asym_fig(Asym(i)+1) = Asym_fig(Asym(i)+1)+1;
end

tem = 0;
for i=70:256
    tem = tem + Asym_fig(i);
end
    
Asym_new = Asym;
for i =1:m*n
    if Asym(i)>70
        Asym_new(i) = 25; 
    end   
end

Asym_new2 = uint8(Normalize(double(Asym_new)));
for i =1:m*n
    if Asym(i)>70
        Asym_new2(i) = 255; 
    end   
end
imshow(Asym_new2)
% imwrite(Asym_new2,'analysis.bmp');
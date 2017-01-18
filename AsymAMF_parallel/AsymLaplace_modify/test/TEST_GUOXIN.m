clear all;
clc;
warning('off')
m = 100;
n = 100;                     %Data size
R = 4;                      % Ture Rank
% rep = 5;
c =3;
K=64;   %initialize cluster number
r=8; %initial rank
Prune =1; %decide whether prune Rank
% randn('seed',2);
    
 
%% normalization  
    
    RU = randn(m,R);
    RV = randn(R,n);
    
    X_Ori = RU * RV;
    W =  rand(size(X_Ori))>0.20;  
   
    X_Noi = X_Ori;
    Ind = randperm(m*n);
    p1 = floor(m*n*0.15);
    p2 = floor(m*n*0.2);
    X_Noi(Ind(1:p1))    =   X_Noi(Ind(1:p1))+randn(1,p1)* 0.25;  %Add Gaussian noise
    X_Noi(Ind(p1+1:p1+p2)) = X_Noi(Ind(p1+1:p1+p2)) + (rand(1,p2)*10) -5; %Add uniform noise
    p3 = floor(m*n*0.2);
    X_Noi(Ind(p1+p2+1:p1+p2+p3)) = X_Noi(Ind(p1+p2+1:p1+p2+p3)) + (rand(1,p3)*4) -2; %Add uniform noise

    U0 =randn(m,r);
    V0 =randn(n,r);
      
%%%%%%%%%%%%% AsymLaplaceAMF method %%%%%%%%%%%%%%%%%%%%%%%%
        opt.tao_s = 0.5; % Æ«ÒÆÁ¿
   
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
      
result = [0,0,0];  

for i  = 0.1:0.1:2
    for j = 0.1:0.1:1
        i
        j
        opt.s_alpha = i;
        opt.s_beta  = j;
        try
            tic;
            [A_sb B_sb X_sb Phi_sb opt] = AsymLaplaceVarInference(X_Noi, W, opt);
            E1_amf = norm((X_Ori - X_sb),'fro')/norm(X_Ori,'fro');
            result = [result [opt.s_alpha;opt.s_beta;E1_amf]];
            time_AMF = toc
        catch err
            disp(err)
            time_AMF = toc
        end
    end 
end


clear all;
clc;
m = 100;
n = 100;                     %Data size
R = 4;                      % Ture Rank
% rep = 5;
c =1;
K=64;   %initialize cluster number
r=8; %initial rank
Prune =1; %decide whether prune Rank
% randn('seed',2);
    
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

    U0 =randn(m,r);
    V0 =randn(n,r);
   
      
        
%%%%%%%%%%%%% AMF method %%%%%%%%%%%%%%%%%%%%%%%%
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
        [A_sb B_sb X_sb Phi_sb opt] = DemoVarInference(X_Noi, W, opt);
        time_AMF = toc;
        E1_amf(i) = norm((X_Ori - X_sb),'fro')/norm(X_Ori,'fro');
        
       
% %%%%%%%%%%%%% BRMF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         opts.maxIter = 100;
%         opts.burnin = 50;
%         opts.invW_0 = 1000 * eye(r);
%         opts.beta_0 = 2;
%         opts.nu_0 = r;
%         opts.a = 1e-4;
%         opts.b = 1e1;
%         opts.r = r;
%         opts.alpha = 0.5;
%         tic;
%         [A_w1 B_w1 Tau_w1 L1] = BRMF(X_Noi, opts,U0,V0);
%         time_brmf = toc;
%     
%         E1_brmf(i) = norm((X_Ori - L1),'fro')/norm(X_Ori,'fro');
        
        
             
   
        
end

        E2_amf = mean(E1_amf);
        std_amf =  std(E1_amf);
%         E2_brmf = mean(E1_brmf);
%         std_brmf= std(E1_brmf);
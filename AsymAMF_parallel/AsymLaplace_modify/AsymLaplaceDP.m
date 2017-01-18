function [opts,Gamma1,Gamma2, Phi, E_inv_sigma1,E_inv_sigma2,E_s] = AsymLaplaceDP(Y, U, V, Cov_U, Cov_V,Gamma1,Gamma2,Phi,E_inv_sigma1,E_inv_sigma2,E_s,Mask,opts)

%%Initialize 
[M, N] = size(Y);

tao_s     = opts.tao_s;
s_alpha        = opts.s_alpha;
s_beta     = opts.s_beta;

alpha0   = opts.alpha0;
a0       = opts.a0;
b0       = opts.b0;
K        = opts.K;
R        = opts.r;     % 初始rank
Pai      = zeros(1,K);
err      = zeros(M,N);   % E(y-uv) = err
err2     = zeros(M,N);   % E(y-uv)^2 = err2
temp0    = (U*V').*Y;
E_inv_s = ones(K,1);   % 对E（1/s）初始化，，这个初始化不影响后续计算，因为在使用前会被更新
E_log_s = zeros(K,1);   % 对E（log(s)）初始化，这个初始化不影响后续计算，因为在使用前会被更新

% 计算err, err2
err = Y - U*V';
N_1 = reshape(V',R,1,N);
N_2 = reshape(V',1,R,N);
temp_N=  bsxfun(@times,N_1,N_2 )+Cov_V; 

M_1 = reshape(U',R,1,M);
M_2 = reshape(U',1,R,M);
temp_M=  bsxfun(@times,M_1,M_2 )+Cov_U ;

MatU =reshape(temp_M, R*R,M);
MatV =reshape(temp_N, R*R,N);
err2= Y.^2+MatU'*MatV-2*temp0;
guoxin_err2 = sum(sum(err2))
clear temp0 temp_M temp_N;
clear MatU;
clear MatV;
parfor i =1:K
    phi(i) =  sum(sum(Phi(:,:,i)));
end
Sum_to_k = sum(phi);
%% Update
% try  

for iter = 1 : 1
%     for k=1:K
%          Gamma(k,1) = 1 + phi(k);
%         Sum_to_k   = Sum_to_k - phi(k);
%         Gamma(k,2) = alpha0+Sum_to_k;
%         psi_vec(k) = psi(Gamma(k,2))-psi(Gamma(k,1)+Gamma(k,2));
%         pai_vec(k) = Gamma(k,2)/(Gamma(k,1)+Gamma(k,2));
%     end

    for x=1:K
        Gamma1(x) =1+phi(x);
        Sum_to_k = Sum_to_k -phi(x);
        Gamma2(x) = alpha0+Sum_to_k;
        if(x==1)
            psi_vec_sum(x) =0;
            pai_vec_prod(x)=1;
        elseif(x==2)
            psi_vec_sum(x) =psi(Gamma2(x-1))-psi(Gamma1(x-1)+Gamma2(x-1));
            pai_vec_prod(x)=Gamma2(x-1)/(Gamma1(x-1)+Gamma2(x-1));
        else
            psi_vec_sum(x) = psi_vec_sum(x-1)+psi(Gamma2(x-1))-psi(Gamma1(x-1)+Gamma2(x-1));
            pai_vec_prod(x) = pai_vec_prod(x-1)*Gamma2(x-1)/(Gamma2(x-1)+Gamma2(x-1));
        end
    end
    Sum_of_psi = zeros(K,1);
    Prod_of_pai = zeros(K,1);
    E_log_sigma = zeros(K,1);
  
    Gamma2(K) = 0;

    parfor k = 1 : K
        warning('off');
        Sum_of_psi(k) =psi_vec_sum(k);  
        Prod_of_pai(k) =pai_vec_prod(k);
%         Tao(k,1)   = a0 + 0.5*phi(k);
%         Tao(k,2)   = b0 + 0.5*sum(sum(Phi(:,:,k).* err));
        A_sigma = 2*a0+2+phi(k);
        B_sigma = (0.5-tao_s)*sum(sum(Phi(:,:,k).* err));
        C_sigma = b0 + 0.5*tao_s*(1-tao_s)* E_s(k) * sum(sum(Phi(:,:,k).* err2));
        J_sigma = J_plus(A_sigma-3,B_sigma,C_sigma);
        E_inv_sigma1(k) = exp(J_plus(A_sigma+1-3,B_sigma,C_sigma)-J_sigma);
        E_inv_sigma2(k) = exp(J_plus(A_sigma+2-3,B_sigma,C_sigma)-J_sigma);
        E_log_sigma(k) = J_log(A_sigma-3,B_sigma,C_sigma,J_sigma);

       
%         [tem0,tem1] = J_log(A_sigma,B_sigma,C_sigma);
%          if isnan(tem1)
%              E_log_sigma = -exp(tem0-J_sigma);
%          else
%              E_log_sigma = -exp(tem0-J_sigma)+exp(tem1-J_sigma);
%          end
    
        %对s进行估计
        alpha_s = tao_s*(1-tao_s)*E_inv_sigma2(k)*sum(sum(Phi(:,:,k).* err2));
%         if alpha_s>1e4
%             alpha_s=1e4
%         end
        beta_s =  2*s_beta+(0.5-tao_s)^2/(tao_s*(1-tao_s)) * phi(k);
        gamma_s = 0.5 * phi(k) - s_alpha;
        sq_ab = sqrt(alpha_s*beta_s);
        r_bessel = div_bessel(gamma_s,sq_ab);
        E_s(k) = sqrt(beta_s/alpha_s)*r_bessel;
        E_inv_s(k) = sqrt(alpha_s/beta_s)*r_bessel - 2*gamma_s/beta_s ;
        E_log_s(k) = 0.5*log(beta_s/alpha_s) + diff_bessel(gamma_s,sq_ab);


    
%         Phi(:,:,k) = exp((psi(Gamma(k,1)) - psi(Gamma(k,1)+Gamma(k,2))+Sum_of_psi) *ones(M,N) - 0.5*Tao(k,1)/Tao(k,2)* err -0.5*(log(Tao(k,2))-psi(Tao(k,1)))*ones(M,N));
        Phi(:,:,k) = exp((psi(Gamma1(k)) - psi(Gamma1(k)+Gamma2(k))+Sum_of_psi(k)) *ones(M,N) -0.5*(0.5-tao_s)^2/(tao_s*(1-tao_s))*E_inv_s(k)*ones(M,N) + (E_log_sigma(k)+0.5*E_log_s(k))*ones(M,N) - 0.5*tao_s*(1-tao_s)*E_inv_sigma2(k)*E_s(k)*err2 + (0.5-tao_s)*E_inv_sigma1(k)*err);   %%change guoxin，此处需要继续修改
        Pai(k)     = Prod_of_pai(k)*Gamma1(k)/(Gamma1(k)+Gamma2(k));
        

    end 
    
%      if opts.verbose
%          disp(Pai);
%      end
    % catch 
     Pai   = Pai/sum(Pai);
     flag  = (Pai>1e-4);   % 相当于将Pai值太小的类删掉
     Pai            = Pai(flag);
     K              = size(Pai,2);
     opts.K         = K;
     E_s = E_s(flag);
     E_inv_s = E_inv_s(flag);
     E_log_s = E_log_s(flag);
     Gamma1          = Gamma1(flag);
     Gamma2          = Gamma2(flag);
     E_inv_sigma1            = E_inv_sigma1(flag,:);
     E_inv_sigma2            = E_inv_sigma2(flag,:);
     Phi            = Phi(:,:,flag);
     Phi = bsxfun(@rdivide, Phi,sum(Phi,3));
     Phi = bsxfun(@times, Phi, Mask);
end
 
end

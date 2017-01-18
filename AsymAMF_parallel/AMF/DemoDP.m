function [opts,Gamma, Phi, Tao] = DemoDP(Y, U, V, Cov_U, Cov_V,Gamma,Phi,Tao,Mask,opts)


%%Initialize 
[M, N] = size(Y);

alpha0   = opts.alpha0;
a0       = opts.a0;
b0       = opts.b0;
K        = opts.K;
R        = opts.r;
Pai      = zeros(1,K);
err      =zeros(M,N);
temp0    = (U*V').*Y;



N_1 = reshape(V',R,1,N);
N_2 = reshape(V',1,R,N);
temp_N=  bsxfun(@times,N_1,N_2 )+Cov_V ;


M_1 = reshape(U',R,1,M);
M_2 = reshape(U',1,R,M);
temp_M=  bsxfun(@times,M_1,M_2 )+Cov_U ;

MatU =reshape(temp_M, R*R,M);
MatV =reshape(temp_N, R*R,N);
err= Y.^2+MatU'*MatV-2*temp0;
guoxin_err2 = sum(sum(err))
clear temp0 temp_M temp_N;
clear MatU;
clear MatV;
for i =1:K
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
        Gamma(x,1) =1+phi(x);
        Sum_to_k = Sum_to_k -phi(x);
        Gamma(x,2) = alpha0+Sum_to_k;
        if(x==1)
            psi_vec_sum(x) =0;
            pai_vec_prod(x)=1;
        elseif(x==2)
            psi_vec_sum(x) =psi(Gamma(x-1,2))-psi(Gamma(x-1,1)+Gamma(x-1,2));
            pai_vec_prod(x)=Gamma(x-1,2)/(Gamma(x-1,1)+Gamma(x-1,2));
        else
            psi_vec_sum(x) = psi_vec_sum(x-1)+psi(Gamma(x-1,2))-psi(Gamma(x-1,1)+Gamma(x-1,2));
            pai_vec_prod(x) = pai_vec_prod(x-1)*Gamma(x-1,2)/(Gamma(x-1,1)+Gamma(x-1,2));
        end
    end
   
    for k = 1 : K
        Sum_of_psi =psi_vec_sum(k);  
        Prod_of_pai =pai_vec_prod(k);
        Tao(k,1)   = a0 + 0.5*phi(k);
        Tao(k,2)   = b0 + 0.5*sum(sum(Phi(:,:,k).* err));
        if k == K
            Gamma(K,2) = 0;
        end

        Phi(:,:,k) = exp((psi(Gamma(k,1)) - psi(Gamma(k,1)+Gamma(k,2))+Sum_of_psi) *ones(M,N) - 0.5*Tao(k,1)/Tao(k,2)* err -0.5*(log(Tao(k,2))-psi(Tao(k,1)))*ones(M,N));
        Pai(k)     = Prod_of_pai*Gamma(k,1)/(Gamma(k,1)+Gamma(k,2));
    end 
    
%      if opts.verbose
%          disp(Pai);
%      end
    % catch 
     Pai   = Pai/sum(Pai);
     flag  = (Pai>1e-4);
     Pai            = Pai(flag);
     K              = size(Pai,2);
     opts.K         = K;
     Gamma          = Gamma(flag,:);
     Tao            = Tao(flag,:);
     Phi            = Phi(:,:,flag);
     Phi = bsxfun(@rdivide, Phi,sum(Phi,3));
%      Phi = bsxfun(@times, Phi, Mask);
end
 
end

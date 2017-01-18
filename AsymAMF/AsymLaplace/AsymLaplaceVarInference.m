function [U V LowRank Phi opts] = AsymLaplaceVarInference(Y,Mask, opts, U, V)
    
    
    tao_s     = opts.tao_s;   
    R           = opts.r;      % 初始rank
    a0          = opts.a0;    
    b0          = opts.b0;
    a1          = opts.a1;
    b1          = opts.b1;
    alpha0      = opts.alpha0;
    K           = opts.K;      % 初始化类中心数目
    [M N]       = size(Y);  
    Gamma    = ones(K,2);
    Gamma(:,2)=alpha0*ones(K,1);
    E_inv_sigma      = a0*ones(K,2); %待修改
    E_inv_sigma(:,2) = b0*ones(K,1);  % 待修改
    Phi      = 1/K *ones(M,N,K);
    E_s      = 4*ones(K,1);   % 对E（s）初始化，因为一开始在计算sigma的分布的时候s的估计还是未知 

    if nargin < 4
       U = 1*randn(M, R) ;
       V = 1*randn(N, R) ;
    end
    LowRank =  U * V';
    Cov_U = repmat( eye(R,R), [1 1 M] );
    Cov_V = repmat( eye(R,R), [1 1 N] );
    if ~isfield(opts, 'tol')
        opts.tol = 1e-4;
    end
    if ~isfield(opts, 'minIter')
        opts.minIter = 20;
    end

    for iter = 1 : opts.maxIter 
        %% Gamma Phi Tao
         
        [opts,Gamma, Phi, E_inv_sigma, E_s] = AsymLaplaceDP(Y, U, V, Cov_U, Cov_V,Gamma,Phi,E_inv_sigma,E_s,Mask,opts);

         
        %% Lambda（yita）

        a1new    = (2*a1 + M + N)*ones(R,1);
        b1new    = (2*b1*ones(R,1)+ diag(V'*V) + diag(sum(Cov_V,3)) + diag(U'*U) + diag(sum(Cov_U,3)));
        lambda_inv = a1new./b1new; % R*1
        lambda     = b1new./a1new;
        LowRank =  U * V';

         if (opts.Prune&&mod(iter,5)==0)
            a      = lambda/sum(lambda);
           [a_sort, index] = sort(a);
            sum_temp     =0;
            t            =R+1;
            
            
           for i =1:R
               t = t-1;
               sum_temp  = sum_temp + a_sort(t);
               if sum_temp >= 0.99
                   break;
               end
           end
%            disp(lambda);
           index  =  index(t:R);
           U      = U(:, index);
           V      = V(:, index);
           lambda_inv = lambda_inv(index);
           R      = size(U,2);
           opts.r = R;
           Cov_U= Cov_U(index, index, :);
           Cov_V= Cov_V(index, index, :);
           disp([num2str(R) ' rank left']);
        end  
        
        Lambda_inv = diag(lambda_inv);
         

           %%  U
             sumnew = reshape(reshape(Phi,M*N,opts.K)*((E_inv_sigma(:,2).*E_s)), M,N);
             Ysum = Y.*sumnew;
             sumnew2 = reshape(reshape(Phi,M*N,opts.K)*E_inv_sigma(:,1),M,N);

             t =reshape(Cov_V,R*R,N);
            for i = 1 : M
                Cov_U(:,:,i) = inv( tao_s*(1-tao_s)*(reshape( t* sumnew(i,:)', R, R ) + bsxfun(@times, V', sumnew(i,:))*V) + Lambda_inv );
                U(i,:) = ((tao_s*(1-tao_s)*Ysum(i,:)- (0.5-tao_s)*sumnew2(i,:))*V) * Cov_U(:,:,i);
            end
            clear t;
           %%  V
            p      =reshape(Cov_U,R*R,M);
            for j = 1 : N
                    Cov_V(:,:,j) = inv( tao_s*(1-tao_s)*(reshape( p* sumnew(:,j), R, R ) + bsxfun(@times, U', sumnew(:,j)')*U) + Lambda_inv );
                    V(j,:) = ((tao_s*(1-tao_s)*Ysum(:,j)'- (0.5-tao_s)*sumnew2(:,j)')*U) * Cov_V(:,:,j);
            end
            clear p;
%              fprintf('U,V updated!\n')
            if norm(U*V' - LowRank , 'fro') / norm(LowRank, 'fro') < opts.tol && iter > opts.minIter
                    break;
            end
      
      
    end
    
    LowRank =  U * V';
    disp([num2str(size(Phi,3)) ' cluster left']);
       
end
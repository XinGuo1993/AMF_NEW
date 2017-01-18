function [U V LowRank Phi opts] = SkewNormVarInference(Y,Mask, opts, U, V)
    
    
    lamda_s     = opts.lamda_s;
    R           = opts.r;
    a0          = opts.a0;
    b0          = opts.b0;
    a1          = opts.a1;
    b1          = opts.b1;
    alpha0      = opts.alpha0;
    K           = opts.K;
    [M N]       = size(Y);
    Gamma    = ones(K,2);
    Gamma(:,2)=alpha0*ones(K,1);
    E_inv_sigma      = a0*ones(K,2); %待修改
    E_inv_sigma(:,2) = b0*ones(K,1);  % 待修改
    Phi      = 1/K *ones(M,N,K);
    E_s      = zeros(K,1);   % 对E（s）初始化，因为一开始在计算sigma的分布的时候s的估计还是未知 

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
         disp(['iter number: ' num2str(iter) ' of ' num2str(opts.maxIter)]);
        %% Gamma Phi Tao
         
        [opts, Gamma, Phi, E_inv_sigma,E_s] = SkewNormDP(Y, U, V, Cov_U, Cov_V,Gamma,Phi,E_inv_sigma,E_s, Mask,opts);

         
          %% Lambda (yita)

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
             sumnew = reshape(reshape(Phi,M*N,opts.K)*E_inv_sigma(:,2), M,N);
             Ysum = Y.*sumnew;
             sumnew2 = reshape(reshape(Phi,M*N,opts.K)*((E_inv_sigma(:,1).* E_s)),M,N);

             t =reshape(Cov_V,R*R,N);
            for i = 1 : M
                Cov_U(:,:,i) = inv( (1+lamda_s^2)*(reshape( t* sumnew(i,:)', R, R ) + bsxfun(@times, V', sumnew(i,:))*V) + Lambda_inv );
                U(i,:) = (((1+lamda_s^2)*Ysum(i,:)-lamda_s*sqrt(1+lamda_s^2)*sumnew2(i,:))*V) * Cov_U(:,:,i);
            end
            clear t;
           %%  V
            p      =reshape(Cov_U,R*R,M);
            for j = 1 : N
                    Cov_V(:,:,j) = inv((1+lamda_s^2)* reshape( p* sumnew(:,j), R, R ) + bsxfun(@times, U', sumnew(:,j)')*U + Lambda_inv );
                    V(j,:) = (((1+lamda_s^2)*Ysum(:,j)'- lamda_s*sqrt(1+lamda_s^2)*sumnew2(:,j)')*U) * Cov_V(:,:,j);
            end
            clear p;

            if norm(U*V' - LowRank , 'fro') / norm(LowRank, 'fro') < opts.tol && iter > opts.minIter
                    break;
            end
           %% error 计算 
            temp0    = (U*V').*Y;
            N_1 = reshape(V',R,1,N);
            N_2 = reshape(V',1,R,N);
            temp_N=  bsxfun(@times,N_1,N_2 )+Cov_V ;

            M_1 = reshape(U',R,1,M);
            M_2 = reshape(U',1,R,M);
            temp_M=  bsxfun(@times,M_1,M_2 )+Cov_U ;

            MatU =reshape(temp_M, R*R,M);
            MatV =reshape(temp_N, R*R,N);
            err2= Y.^2+MatU'*MatV-2*temp0;   %  size M*N
            guoxin_err2 = sum(sum(err2))
            clear temp0 temp_M temp_N;
            clear MatU;
            clear MatV;
    end
    
    LowRank =  U * V';
    disp([num2str(size(Phi,3)) ' cluster left']);
       
end
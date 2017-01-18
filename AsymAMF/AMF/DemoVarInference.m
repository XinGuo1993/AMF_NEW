function [U V LowRank Phi opts] = DemoVarInference(Y,Mask, opts, U, V)
    
    
    
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
    Tao      = a0*ones(K,2);
    Tao(:,2) = b0*ones(K,1);
    Phi      = 1/K *ones(M,N,K);

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
         
        [opts, Gamma, Phi, Tao] = DemoDP(Y, U, V, Cov_U, Cov_V,Gamma,Phi,Tao,Mask,opts);
         E_sigma     = (Tao(:,1)./Tao(:,2))';
         
          %% Lambda 

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

             sumnew = reshape(reshape(Phi,M*N,opts.K)*E_sigma', M,N);
             Ysum = Y.*sumnew;

             t =reshape(Cov_V,R*R,N);
            for i = 1 : M
                Cov_U(:,:,i) = inv( reshape( t* sumnew(i,:)', R, R ) + bsxfun(@times, V', sumnew(i,:))*V + Lambda_inv );
                U(i,:) = (Ysum(i,:)*V) * Cov_U(:,:,i);
            end
            clear t;
           %%  V
            p      =reshape(Cov_U,R*R,M);
            for j = 1 : N
                    Cov_V(:,:,j) = inv( reshape( p* sumnew(:,j), R, R ) + bsxfun(@times, U', sumnew(:,j)')*U + Lambda_inv );
                    V(j,:) = (Ysum(:,j)'*U) * Cov_V(:,:,j);
            end
            clear p;

            if norm(U*V' - LowRank , 'fro') / norm(LowRank, 'fro') < opts.tol && iter > opts.minIter
                    break;
            end
      
      
    end
    
    LowRank =  U * V';
    disp([num2str(size(Phi,3)) ' cluster left']);
       
end
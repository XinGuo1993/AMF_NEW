function Isum = quadToInf(fun,a,dx0,tol,method)
if nargin < 2 ,a=0 ;end
if nargin < 3 ,dx0=0.5 ;end
if nargin < 4 ,tol = 5e-4 ;end
if nargin < 5 ,method = 1 ;end
j=0;dx = dx0;Isum = 0;x1 = a; maxint = 35;
%fprintf('\n     j       dx      x2      I_j     Isum\n');
while j<maxint
    x2 = x1 + dx;
    switch method
        case 1, I = quad(fun,x1,x2);
        case 2, I = quadl(fun,x1,x2);
        otherwise, error(sprintf('method = %d not allowed',method));
    end    
    Isum = Isum + I;    
    %fprintf('%4d %8.1f %8.1f %12.8f %12.8f\n',j,dx,x2,I,Isum);
    if j>5 & abs(I/Isum) < tol,break; end
    j = j+1;x1 = x2;dx = 2*dx;
end
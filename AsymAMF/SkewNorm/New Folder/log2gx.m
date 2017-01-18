clear;
clc;
%% 确定积分界
A = 500;
B = 0;
C = 500;
eps = 1e-15;
mu0 = (B+sqrt(B^2 +8*A*C))/(4*C);
sigma0 = (B+sqrt(B^2 +8*A*C))/sqrt(4*C*(8*A*C+B^2+B*sqrt(B^2 +8*A*C) ));
low_bound = -mu0/(sqrt(2)*sigma0);
upper_bound = 1000;
% f=@(x)A*log(x)+B*x-C*x^2;
f=@(x)exp( (B-2*C*mu0)*sqrt(2)*sigma0*x - 2*C*sigma0^2 *x.^2 + A*log(1+ sqrt(2)*sigma0/mu0 * x));
L = low_bound + 0.001;
U = upper_bound-0.001;
Lower = f(L);
Upper = f(U);
%% 有进一步优化的必要
while Lower>eps
    L = 0.5*(low_bound+L)
    Lower = f(L);
end
while Upper>eps
    U = 2*U
    Upper = f(U);
end
%%
L
U
Lower
Upper

Idx = (round(U-L)+10)*10;
step = (U-L)/Idx;
result = 0;
for i = 1:step-1
    result = result + step/6*(f(L+step*(i-1)) + 4*f(L+step*(i-0.5)) + f(L+step*i));
end
result

J_plus(A,B,C)




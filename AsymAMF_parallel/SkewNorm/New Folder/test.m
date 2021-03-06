clear;
A = 627;
B = 0;
C = 14802;

mu0 = (B+sqrt(B^2 +8*A*C))/(4*C);
sigma0 = (B+sqrt(B^2 +8*A*C))/sqrt(4*C*(8*A*C+B^2+B*sqrt(B^2 +8*A*C) ));
low_bound = -mu0/(sqrt(2)*sigma0);
%quadl('f(x,0.6076,0.3389,1,2,3)',1,2)
%quad(['f(x,',num2str(mu0),',',num2str(sigma0),',',num2str(A),',',num2str(B),',',num2str(C),')'],-1,1)
I0 = quadToInf(['f(x,',num2str(mu0),',',num2str(sigma0),',',num2str(A),',',num2str(B),',',num2str(C),')'],low_bound);
if imag(I0) <0.01
    I = log(sqrt(2)*sigma0) + log(real(I0)) + A*log(mu0)+B*mu0-C*mu0^2;
end

guoxin = quadToInf(['f2(x,',num2str(A),',',num2str(B),',',num2str(C),')']);

% Çó½âmu ºÍsiama
A = 1;
B = 2;
C = 3;
mu1 = 4;
mu2 = 5;
sigma1 = 1/sqrt(1/mu1^2*(1/log(mu1)+1/log(mu1)^2+A) +2*C);
sigma2  = 1/sqrt(1/mu2^2*(1/log(mu2)+1/log(mu2)^2+A) +2*C);

upper_bound1 =  1 - mu1/(sqrt(2)*sigma1);
low_bound1 = - mu1/(sqrt(2)*sigma1);
low_bound2 = 1 - mu2/(sqrt(2)*sigma2);

I1 = quadl(['f2_1(x,',num2str(mu1),',',num2str(sigma1),',',num2str(A),',',num2str(B),',',num2str(C),')'],low_bound1,upper_bound1);
I2 = quadToInf(['f2_2(x,',num2str(mu2),',',num2str(sigma2),',',num2str(A),',',num2str(B),',',num2str(C),')'],low_bound2);

I11 = log(sqrt(2)*sigma1) + log(real(I1)) + A*log(mu1)+B*mu1-C*mu1^2;
I22 = log(sqrt(2)*sigma2) + log(real(I2)) + A*log(mu2)+B*mu2-C*mu2^2;
%% 确定mu0
f=@(x)1./(x.*log(x))+A./x+B-2*C.*x;
a = 0.001;
b = 1-0.001;
c = 1+0.001;
d = 10;
threshold = 0.0001;
%确定 初始点，由于函数单调，因而用二分法求零点最准确（直接使用fzero函数出现了无法收链到大于一解的情况）
fa=f(a);
fb=f(b);
fc=f(c);
fd=f(d);
while fa<0
    a = a/2;
    fa=f(a);
end
while fb>0
    b = (1+b)*0.5;
    fb=f(b);
end
while fc<0
    c = (1+c)*0.5;
    fc=f(c);
end
while fd>0
    d = 2*d;
    fd=f(d);
end
while min(abs(fa),abs(fb))>threshold
    tem = (a+b)*0.5;
    f_tem = f(tem);
    if f_tem>0
        a = tem;
        fa = f_tem;
    else
        b = tem;
        fb = f_tem; 
    end
    
end

while min(abs(fc),abs(fd))>threshold
    tem = (c+d)*0.5;
    f_tem = f(tem);
    if f_tem>0
        c = tem;
        fc = f_tem;
    else
        d = tem;
        fd = f_tem;
    end
end

if abs(fa) < abs(fb)
    tem1 = a;
else
    tem1 = b;
end
if abs(fc) < abs(fd)
    tem2 = c;
else
    tem2 = d;
end
% fprintf('erfen ending\n');

%%
mu1 = tem1;
mu2 = tem2;
sigma1 = 1/sqrt(1/mu1^2*(1/log(mu1)+1/log(mu1)^2+A) +2*C);
sigma2 = 1/sqrt(1/mu2^2*(1/log(mu2)+1/log(mu2)^2+A) +2*C);

upper_bound1 =  (1-1e-7 - mu1)/(sqrt(2)*sigma1);
low_bound1 = (1e-7 - mu1)/(sqrt(2)*sigma1);
low_bound2 = (1+1e-7 - mu2)/(sqrt(2)*sigma2);

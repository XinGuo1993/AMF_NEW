function result = J_log(A,B,C,J_sigma)
%UNTITLED8 此处显示有关此函数的摘要
% A x指数， B exp中x指数， C exp中x^2指数
%%%%% 注意输出的是log值，而不是原函数
%   此处显示详细说明
flag = 0;
if A>0.25
    flag =1;
    %% 确定mu0
    f=@(x)1./(x.*log(x))+A./x+B-2*C.*x;
    a = 1e-8;
    b = 1-1e-8;
    c = 1+1e-8;
    d = 10;
    threshold = 1e-8;
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

    fa_abs = abs(fa);
    fb_abs = abs(fb);
    fa2 = fa_abs+10;
    fb2 = fb_abs+10;
    while (min(fa_abs,fb_abs)>threshold) && (abs(fa2-fa_abs)+abs(fb2-fb_abs)>1e-8)
        fa2 = fa_abs;
        fb2 = fb_abs;
        tem = (a+b)*0.5;
        f_tem = f(tem);
        if f_tem>0
            a = tem;
            fa = f_tem;
            fa_abs = abs(fa);
        else
            b = tem;
            fb = f_tem; 
            fb_abs = abs(fb);
        end
    end

    fc_abs = abs(fc);
    fd_abs = abs(fd);
    fc2 = fc_abs+10;
    fd2 = fd_abs+10;
    while (min(fc_abs,fd_abs)>threshold) && (abs(fc2-fc_abs)+abs(fd2-fd_abs)>1e-8)
        fc2 = fc_abs;
        fd2 = fd_abs;
        tem = (c+d)*0.5;
        f_tem = f(tem);
        
        if f_tem>0
            c = tem;
            fc = f_tem;
            fc_abs = abs(fc);
        else
            d = tem;
            fd = f_tem;
            fd_abs = abs(fd);
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

    %%
    mu1 = tem1;
    mu2 = tem2;
    sigma1 = 1/sqrt(1/mu1^2*(1/log(mu1)+1/log(mu1)^2+A) +2*C);
    sigma2 = 1/sqrt(1/mu2^2*(1/log(mu2)+1/log(mu2)^2+A) +2*C);
    
 
  if isreal(sigma1)&&isreal(sigma1)
        upper_bound1 =  (1 - mu1)/(sqrt(2)*sigma1);
        low_bound1 = (0 - mu1)/(sqrt(2)*sigma1);
        low_bound2 = (1 - mu2)/(sqrt(2)*sigma2);

        I1 = quadl(['f_log(x,',num2str(mu1),',',num2str(sigma1),',',num2str(A),',',num2str(B),',',num2str(C),')'],low_bound1,upper_bound1);
        I2 = quadToInf(['f_log(x,',num2str(mu2),',',num2str(sigma2),',',num2str(A),',',num2str(B),',',num2str(C),')'],low_bound2);

        if imag(I1) <0.01
            I11 = log(sqrt(2)*sigma1) + log(real(I1)) + log(-log(mu1))+A*log(mu1)+B*mu1-C*mu1^2;
        end
        
        if imag(I2) <0.01
            I22 = log(sqrt(2)*sigma2) + log(real(I2)) + log(log(mu2))+A*log(mu2)+B*mu2-C*mu2^2;
            result = -exp(I11-J_sigma)+exp(I22-J_sigma);
             
        elseif abs(mu2-1)<1e-8
            I22=0;
            result = -exp(I11-J_sigma);
        else
            I22=quadl(['f_log2(x,',num2str(A),',',num2str(B),',',num2str(C),')'],1,mu2);
            result = -exp(I11-J_sigma);
        end 
  else
      flag =0;
  end

end

if flag ==0
    I22=quadToInf(['f_log2(x,',num2str(A),',',num2str(B),',',num2str(C),')']);
    result = I22/exp(J_sigma);
end

end


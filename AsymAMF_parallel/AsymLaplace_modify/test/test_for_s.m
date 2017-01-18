%测试收敛s的参数
clear;
clc;

beta_s = 1;
err = 10;
k=0;
alpha(1)=100000;
gamma(1)=100000;
save('guoxin.mat','alpha');
save('guoxin.mat','gamma','-append');

for alpha_s = 0.002:0.00001:0.004
    for gamma_s = -0.940: 0.00001:-0.90;
        sq_ab = sqrt(alpha_s*beta_s);
        E_s = sqrt(beta_s/alpha_s)*besselk(gamma_s+1,sq_ab)/besselk(gamma_s,sq_ab);
%         E_inv_s = sqrt(alpha_s/beta_s)*besselk(gamma_s+1,sq_ab)/besselk(gamma_s,sq_ab) - 2*gamma_s/beta_s;
        E_log_s = 0.5*log(beta_s/alpha_s) + diff_bessel(gamma_s,sq_ab);
        err = (E_s-4)^2+E_log_s^2;
        if(err<0.0000001)
            k=k+1;
            alpha_s
            gamma_s
            alpha(k)=alpha_s;
            gamma(k)=gamma_s;
            save('guoxin.mat','alpha','-append');
            save('guoxin.mat','gamma','-append');
        end
    end
end
clear
load('guoxin.mat')

function f_log = f_log(x,mu0,sigma0,A,B,C)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
% function for J_log
%      f_log = log(mu0+sqrt(2)*sigma0.*x)./log(mu0).*(1+sqrt(2)*sigma0/mu0.*x).^A .*exp(sqrt(2)*B*sigma0.*x - 2*sigma0^2*C.*x.^2 -2*sqrt(2)*mu0*sigma0*C.*x); 
      f_log = exp( log(log(mu0+sqrt(2)*sigma0.*x)/log(mu0)) + (B-2*C*mu0)*sqrt(2)*sigma0*x - 2*C*sigma0^2 *x.^2 + A*log(1+ sqrt(2)*sigma0/mu0 * x)  );
end



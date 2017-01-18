% test
clear;
clc;
mu = 1;
Z = 100;
BES = besselk(mu+1,Z)./besselk(mu,Z)
result = div_bessel(mu,Z)

% BES2 = exp(diff_bessel(mu,Z));
% plot(mu,BES)
% hold on;
% plot(mu,BES2)









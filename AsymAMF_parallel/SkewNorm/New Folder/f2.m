function f2 = f2(x,A,B,C)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
% function for J_log
    f2 =  exp(B*x - C*x.^2 + A*log(x) + log(log(x)));
end



function f_log2 = f_log2(x,A,B,C)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    f_log2 =  log(x).*exp(B*x - C*x.^2 + A*log(x));
end
